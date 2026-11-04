# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Genome-backed QC and bounded local rescue of an existing insect OR catalogue.

Does not alter a genome or the input catalogue. Family support, coding-model
state and evidence confidence are separate. All accepted models must round-trip
from genomic coordinates to CDS and protein; pseudogene translations retain *.
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import hashlib
import json
import logging
import math
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote

from .audit_rescue_endpoints import endpoint_options, oriented_bases, transcript, write_table
from .rescue_flagged_or import overlap_fraction, read_fai, read_table
from .summarize_rescue_exonerate import parse_output
from orsearcher.fasta import FastaRecord, read_fasta, reverse_complement, translate, write_fasta
from orsearcher.models import parse_attributes, parse_miniprot_gff, attach_reference_lengths
from orsearcher.sharding import Segment

LOG = logging.getLogger("catalogue_review")
CANONICAL = {"GT-AG", "GC-AG", "AT-AC"}


def command(args, output=None, timeout=120):
    if output:
        with Path(output).open("w") as stream:
            result = subprocess.run(args, stdout=stream, stderr=subprocess.PIPE, text=True, timeout=timeout)
    else:
        result = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    if output:
        Path(str(output) + ".stderr.log").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"exit={result.returncode}: {' '.join(args)}\n{result.stderr[-2000:]}")
    return result.stdout or ""


def annotations(path, wanted):
    out = {cid: [] for cid in wanted}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            parent = parse_attributes(fields[8]).get("Parent", "")
            if parent in out:
                out[parent].append((int(fields[3]), int(fields[4]), int(fields[7]) if fields[7].isdigit() else 0))
    if any(not values for values in out.values()):
        raise ValueError("Some selected candidates lack CDS coordinates")
    return out


def splice_sites(sequence, exons, strand):
    ordered = sorted(exons, reverse=strand == "-")
    sites = []
    for left, right in zip(ordered, ordered[1:]):
        if strand == "+":
            donor, acceptor = sequence[left[1]:left[1]+2], sequence[right[0]-3:right[0]-1]
        else:
            donor = reverse_complement(sequence[left[0]-3:left[0]-1])
            acceptor = reverse_complement(sequence[right[1]:right[1]+2])
        sites.append(f"{donor}-{acceptor}")
    return sites


def build_model(cid, source, seq, offset, exons, strand, phase=0, frameshifts=0,
                ref=None, qstart=0, qend=0, ref_length=0, identity=0, terminal=None):
    exons = sorted(exons, reverse=strand == "-")
    if phase:
        if strand == "+": exons[0] = (exons[0][0] + phase, exons[0][1])
        else: exons[0] = (exons[0][0], exons[0][1] - phase)
    if any(s < 1 or e > len(seq) or s > e for s, e in exons):
        raise ValueError("Invalid CDS coordinates")
    cds, positions = transcript(seq, exons, strand)
    protein = translate(cds).removesuffix("*")
    sites = splice_sites(seq, exons, strand)
    stops = []
    for index, aa in enumerate(protein):
        if aa != "*": continue
        local = positions[index*3:index*3+3]
        distance = min(min(p-(s-1), e-1-p) for p in local for s,e in exons if s-1 <= p <= e-1)
        stops.append(dict(codon=cds[index*3:index*3+3], aa_position=index+1,
                          coordinates=[p+offset+1 for p in local], edge_distance=distance,
                          ambiguous_context=sum(b not in "ACGT" for b in seq[max(0,min(local)-100):max(local)+101])))
    key = json.dumps([cid, exons, strand, source, ref, terminal], sort_keys=True)
    return dict(id=cid+".m"+hashlib.sha256(key.encode()).hexdigest()[:12], cid=cid, source=source,
                exons=[[s+offset,e+offset] for s,e in exons], strand=strand, cds=cds, protein=protein,
                length=len(protein), has_start=cds.startswith("ATG"), has_stop=translate(cds).endswith("*"),
                remainder=len(cds)%3, ambiguous=sum(b not in "ACGT" for b in cds), stops=stops,
                frameshifts=frameshifts, sites=sites, canonical=all(site in CANONICAL for site in sites),
                ref=ref or "", ref_cov=(qend-qstart+1)/ref_length if ref_length else 0,
                identity=identity, terminal=terminal or {}, overlap=0., bridge=False)


def terminal_variants(model, seq, offset, references):
    exons = [(s-offset,e-offset) for s,e in model["exons"]]
    cds, positions = transcript(seq, exons, model["strand"])
    starts, added = endpoint_options(seq, cds, positions, model["strand"], max_extension=450, max_trim=30)
    if added is None or model["frameshifts"] or model["stops"] or model["remainder"]:
        return []
    out=[]
    for shift in starts:
        if shift == 0 and added == 0: continue
        new = list(exons)
        if model["strand"] == "+":
            new[0]=(new[0][0]+shift,new[0][1]); new[-1]=(new[-1][0],new[-1][1]+added)
        else:
            new[0]=(new[0][0],new[0][1]-shift); new[-1]=(new[-1][0]-added,new[-1][1])
        if any(s>e for s,e in new): continue
        term=dict(start_adjustment_nt=shift, stop_extension_nt=added)
        if model["ref"]:
            metrics=model.get("alignment",{})
            term["n_difference"]=abs(-shift/3-(metrics.get("qstart",1)-1))
            term["c_difference"]=abs((added-3)/3-(metrics.get("qlen",0)-metrics.get("qend",0))) if added else 0
            if term["n_difference"]>30 or term["c_difference"]>30:
                continue
        try:
            new_model=build_model(model["cid"],model["source"]+"_terminal",seq,offset,new,model["strand"],
                                  frameshifts=model["frameshifts"],ref=model["ref"],terminal=term)
        except ValueError:
            continue
        new_model.update(ref_cov=model["ref_cov"],identity=model["identity"],alignment=model.get("alignment",{}))
        if model.get('alignment_tool'):new_model['alignment_tool']=model['alignment_tool']
        out.append(new_model)
    return out


def model_candidates(row, case, exons, protein, refs, lengths, args, original_model=None):
    cid=row["candidate_id"]
    start=max(1,int(row["start"])-args.flank_bp)
    end=min(lengths[row["sequence_id"]],int(row["end"])+args.flank_bp)
    case.mkdir(parents=True,exist_ok=True)
    fa=command([args.samtools,"faidx",str(args.genome),f"{row['sequence_id']}:{start}-{end}"])
    seq="".join(fa.splitlines()[1:]).upper(); offset=start-1
    if len(seq)!=end-start+1: raise ValueError("Region extraction length mismatch")
    write_fasta(case/"region.fna",[FastaRecord("region","",seq)])
    (case/"region.json").write_text(json.dumps(dict(offset=offset,end=end,contig=row["sequence_id"],strand=row["strand"])))
    ordered=sorted(exons,key=lambda item:item[0],reverse=row["strand"]=="-")
    orig=build_model(cid,"original",seq,offset,[(s-offset,e-offset) for s,e,_ in ordered],row["strand"],
                     phase=ordered[0][2],frameshifts=int(row["frameshifts"]))
    # Endpoint proposals must retain the original homolog endpoints BEFORE
    # terminal_variants runs. A post-hoc update silently permits arbitrary ATGs.
    if original_model is not None:
        for key in ['ref','ref_cov','identity','alignment','alignment_tool']:
            if key in original_model:orig[key]=original_model[key]
    orig.setdefault('alignment_tool','miniprot')
    orig["protein_match"]=orig["protein"].rstrip("*")==protein.rstrip("*")
    all_models=[orig,*terminal_variants(orig,seq,offset,refs)]
    query=list(refs.values())
    write_fasta(case/"references.faa",query)
    if query:
        command([args.miniprot,"-t",str(args.threads_per_job),"-k","5","-M","0","-L","15","-l","4","--gff","-G",str(args.max_intron),
                 "--outc","0.15","--outs","0.20","-N","12","--outn","3","-n","3","-p","0.2",
                 str(case/"region.fna"),str(case/"references.faa")],case/"miniprot.gff3",args.miniprot_timeout)
        segment=Segment("region",row["sequence_id"],"",offset,end,lengths[row["sequence_id"]],"region.fna",False)
        aligned=parse_miniprot_gff(case/"miniprot.gff3",{"region":segment})
        attach_reference_lengths(aligned,{r.id:len(r.sequence) for r in query})
        for item in aligned:
            if item.strand!=row["strand"]: continue
            features=sorted([f for f in item.features if f.feature_type=="CDS"],key=lambda f:f.start,reverse=item.strand=="-")
            if not features: continue
            new=build_model(cid,"miniprot",seq,offset,[(f.start,f.end) for f in features],item.strand,
                            int(features[0].phase) if features[0].phase.isdigit() else 0,
                            item.frameshifts,item.reference_id,item.reference_start,item.reference_end,
                            item.reference_length,item.identity)
            new["alignment"]=dict(qstart=item.reference_start,qend=item.reference_end,qlen=item.reference_length)
            new['alignment_tool']='miniprot'
            if overlap_fraction([(s,e) for s,e,_ in exons],new["exons"])<0.5: continue
            all_models.extend([new,*terminal_variants(new,seq,offset,refs)])
    for model in all_models:
        model["overlap"]=overlap_fraction([(s,e) for s,e,_ in exons],model["exons"])
    (case/"models.json").write_text(json.dumps(all_models))
    return all_models


def parse_hmm(path):
    out={}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"): continue
        f=line.split()
        out[f[0]]=dict(evalue=float(f[4]),score=float(f[5]),domain_score=float(f[8]))
    return out


def assess_models(models, args, stage, reference_records):
    stage.mkdir(parents=True,exist_ok=True)
    proteins=stage/"models.faa"
    write_fasta(proteins,[FastaRecord(m["id"],"",m["protein"].replace("*","X")) for m in models if m["protein"]])
    # GA is evidence, not a universal rejection criterion. Sub-GA results are
    # retained and may be supported by strong competing protein homology.
    command([args.hmmsearch,"--cpu",str(args.workers),"--noali","--tblout",str(stage/"or.tbl"),
             str(args.or_hmm),str(proteins)],stage/"or.out")
    or_hits=parse_hmm(stage/"or.tbl")
    gr_hits={}
    for index,hmm in enumerate(args.gr_hmm):
        path=stage/f"gr{index}.tbl"
        command([args.hmmsearch,"--cpu",str(args.workers),"--cut_ga","--noali","--tblout",str(path),
                 str(hmm),str(proteins)],stage/f"gr{index}.out")
        gr_hits.update(parse_hmm(path))
    ga=domain_ga=None
    for line in args.or_hmm.read_text().splitlines():
        if line.startswith("GA "):
            ga=float(line.split()[1]); domain_ga=float(line.split()[2].rstrip(";")); break
    hits={"or":defaultdict(list),"gr":defaultdict(list)}
    for family,db,limit in [("or",args.or_db,25),("gr",args.gr_db,1)]:
        path=stage/f"{family}.diamond.tsv"
        command([args.diamond,"blastp","--query",str(proteins),"--db",str(db),"--out",str(path),
                 "--outfmt","6","qseqid","sseqid","pident","length","qlen","slen","qstart","qend",
                 "sstart","send","evalue","bitscore","--max-target-seqs",str(limit),"--sensitive",
                 "--evalue","1e-3","--threads",str(args.workers*args.threads_per_job),"--quiet"])
        for line in path.read_text().splitlines():
            f=line.split("\t")
            hits[family][f[0]].append(dict(ref=f[1],identity=float(f[2]),length=int(f[3]),qlen=int(f[4]),
                                          slen=int(f[5]),qs=int(f[6]),qe=int(f[7]),ss=int(f[8]),se=int(f[9]),
                                          evalue=float(f[10]),bits=float(f[11])))
    for model in models:
        mid=model["id"]
        hm=or_hits.get(mid,{"evalue":1.,"score":0.,"domain_score":0.})
        oh=hits["or"].get(mid,[]); gh=hits["gr"].get(mid,[])
        ob=oh[0]["bits"] if oh else 0.; gb=gh[0]["bits"] if gh else 0.
        model.update(or_evalue=hm["evalue"],or_score=hm["score"],or_ga=hm["score"]>=(ga or 1e9) and hm["domain_score"]>=(domain_ga or 1e9),
                     gr_ga=mid in gr_hits,or_bits=ob,gr_bits=gb,
                     family="GR_likely" if gb-ob>=10 and gb else "OR" if oh and ob-gb>=10 else "OR_GR_uncertain")
        model["full_ref_support"]=[]
        model["reference_exact_matches"]=[]
        for hit in oh:
            ref=reference_records.get(hit["ref"])
            if ref and ref.sequence==model["protein"]:
                model["reference_exact_matches"].append(hit["ref"])
                # This supports compatibility with a known OR, but cannot be
                # counted as independent evidence for completeness confidence.
            # Divergent insect OR N termini have variable lengths. Relative
            # endpoint bounds avoid downgrading an intact gene solely because
            # a distant reference is longer or its initial residue is not M.
            subject_slack=max(20,math.ceil(hit["slen"]*.1))
            query_slack=max(20,math.ceil(model["length"]*.1))
            if (ref and hit["evalue"]<=1e-5 and hit["identity"]>=25
                    and (hit["se"]-hit["ss"]+1)/hit["slen"]>=0.8 and hit["ss"]<=subject_slack
                    and hit["slen"]-hit["se"]<=subject_slack and .75<=model["length"]/hit["slen"]<=1.25
                    and hit["qs"]<=query_slack and model["length"]-hit["qe"]<=query_slack):
                model["full_ref_support"].append(hit["ref"])
        if model["ref"] and model["ref"] in reference_records:
            ref=reference_records[model["ref"]]
            a=model.get("alignment",{})
            term=model.get("terminal",{})
            if (ref.sequence.startswith("M") and ref.sequence!=model["protein"] and model["ref_cov"]>=.8
                    and model["identity"]>=.25 and .75<=model["length"]/len(ref.sequence)<=1.25
                    and (term.get("n_difference",a.get("qstart",999)-1))<=30
                    and (term.get("c_difference",len(ref.sequence)-a.get("qend",0)))<=30):
                model["full_ref_support"].append(model["ref"])
        model["full_ref_support"]=sorted(set(model["full_ref_support"]))


def reviewed_results_path(run):
    """Accept a normal reusable review run or an older finalized cache layout."""
    for folder in [run, run/"reviewed_results"]:
        if (folder/"all_models.json").is_file() and (folder/"final/OR_catalogue.tsv").is_file():
            return folder
    raise ValueError(f"No review models/catalogue found under {run}")


def coding_intact(model):
    return not (model["frameshifts"] or model["stops"] or model["ambiguous"] or model["remainder"])


def complete_model(model):
    return (coding_intact(model) and model["has_start"] and model["has_stop"] and model["canonical"]
            and model["overlap"]>=.8 and not model["bridge"] and (model["family"]=="OR" or model.get('structure_or_eligible',False))
            and bool(model["full_ref_support"]) and 250<=model["length"]<=650
            and (model["or_ga"] or model["or_evalue"]<=1e-3 and model["or_bits"]>=100
                 or model.get('family_support')=='curated_homology_consensus' and len(model['full_ref_support'])>=2
                 or model.get('family_support')=='curated_full_homology_consensus'))


def model_rank(model):
    return (not model["bridge"],model.get("family")=="OR",model.get("or_ga",False),coding_intact(model),
            model["canonical"],model["overlap"],min(model["ref_cov"],1),model.get("or_bits",0),model["length"])


def decide(row, models):
    good=[m for m in models if complete_model(m)]
    if good:
        # Prefer an unchanged, independently supported original over arbitrary
        # length maximization. Otherwise rank endpoints and homolog evidence.
        good.sort(key=lambda m:(m["source"]=="original",len(m["full_ref_support"]),m["overlap"],
                                m.get("or_bits",0),-sum(m["terminal"].get(k,0) for k in ["n_difference","c_difference"])),reverse=True)
        chosen=good[0]
        concordant={m["ref"] for m in models if m["source"].startswith(("miniprot","exonerate"))
                    and m["cds"]==chosen["cds"] and coding_intact(m) and m["ref"]}
        sources={m["source"].split("_")[0] for m in models if m["cds"]==chosen["cds"] and coding_intact(m)}
        exact=set(chosen.get("reference_exact_matches",[]))
        independent=len((set(chosen["full_ref_support"])|concordant)-exact)
        confidence="reviewed_complete" if chosen["source"]=="original" and independent>=2 or "exonerate" in sources and "miniprot" in sources else "provisional_complete"
        return dict(category="complete",confidence=confidence,model=chosen,
                    reason="genomic_ATG_stop_intact_frame_canonical_splicing_external_full_reference_support",
                    stop_reference_count=0,stop_coordinates="",crosscheck_sources=",".join(sorted(sources)))
    trusted=[m for m in models if not m["bridge"] and m.get("family")=="OR" and m["canonical"]
             and m["overlap"]>=.75 and m["ref_cov"]>=.75 and m["ref"] and not m["ambiguous"]]
    stop_support=defaultdict(set); stop_models=defaultdict(list)
    for model in trusted:
        if model["frameshifts"]:
            continue  # frame-dependent downstream stops do not confirm a shared nonsense mutation
        for stop in model["stops"]:
            if stop["edge_distance"]>=30 and stop["ambiguous_context"]==0:
                key=",".join(map(str,stop["coordinates"]))
                stop_support[key].add(model["ref"]); stop_models[key].append(model)
    key=max(stop_support,key=lambda k:len(stop_support[k])) if stop_support else ""
    undisturbed=[m for m in trusted if coding_intact(m) and m["length"]>=250]
    if key and len(stop_support[key])>=3 and not undisturbed:
        chosen=max(stop_models[key],key=lambda m:(not m["frameshifts"],-len(m["stops"]),m["ref_cov"],m.get("or_bits",0)))
        tools={m["source"].split("_")[0] for m in stop_models[key]}
        return dict(category="pseudogene",confidence="genome_supported" if "exonerate" in tools else "multi_reference_supported",
                    model=chosen,reason="shared_exonic_stop_no_ambiguous_bases_no_credible_undisrupted_alternative",
                    stop_reference_count=len(stop_support[key]),stop_coordinates=key,crosscheck_sources=",".join(sorted(tools)))
    orig=next(m for m in models if m["source"]=="original")
    if row["original_category"]=="pseudogene" and (orig["stops"] or orig["frameshifts"]) and not orig["ambiguous"] and not undisturbed:
        return dict(category="pseudogene",confidence="provisional_pseudogene",model=orig,
                    reason="original_disruption_retained_but_insufficient_independent_structure_evidence",
                    stop_reference_count=len(stop_support[key]) if key else 0,stop_coordinates=key,crosscheck_sources="original")
    alternatives=[m for m in models if m["source"]!="original" and not m["bridge"] and m["family"]=="OR"
                  and coding_intact(m) and m["canonical"] and m["overlap"]>=.8]
    longer=[m for m in alternatives if m["length"]>=orig["length"]+40]
    chosen=max(longer,key=model_rank) if longer else orig
    return dict(category="partial",confidence="unresolved",model=chosen,
                reason="family_supported_but_terminal_or_coding_structure_not_resolved",
                stop_reference_count=len(stop_support[key]) if key else 0,stop_coordinates=key,crosscheck_sources="")


def mark_bridges(models, rows, original_exons):
    anchors=[r for r in rows if r["original_category"]=="complete" and int(r["protein_length"])>=250]
    byid={r["candidate_id"]:r for r in rows}
    for model in models:
        row=byid[model["cid"]]
        for other in anchors:
            oid=other["candidate_id"]
            if oid==model["cid"] or other["sequence_id"]!=row["sequence_id"] or other["strand"]!=row["strand"]: continue
            target=[(s,e) for s,e,_ in original_exons[oid]]
            if overlap_fraction(target,model["exons"])>=.5:
                model["bridge"]=True; model["bridge_to"]=oid; break


def second_check(row, case, models, refs, args):
    cid=row["candidate_id"]; output=case/"exonerate"; output.mkdir(exist_ok=True)
    sequence=next(read_fasta(case/"region.fna")).sequence
    meta=json.loads((case/"region.json").read_text()); offset=meta["offset"]
    viable=[m for m in models if not m["bridge"] and m["ref"] and m.get("family")=="OR"]
    selected=sorted(viable,key=model_rank,reverse=True)[:6]
    refids=list(dict.fromkeys(m["ref"] for m in selected))[:3]
    if not refids: refids=list(refs)[:3]
    if not refids: return [],dict(cid=cid,status="no_reference")
    target=next((m for m in selected if m["ref"]==refids[0]),models[0])
    left=max(offset+1,min(s for s,e in target["exons"])-5000)
    right=min(offset+len(sequence),max(e for s,e in target["exons"])+5000)
    if right-left+1>250000: return [],dict(cid=cid,status="window_too_large")
    new_offset=left-1; seq=sequence[left-offset-1:right-offset]
    write_fasta(output/"target.fna",[FastaRecord("region","",seq)])
    write_fasta(output/"query.faa",[refs[r] for r in refids])
    cmd=[args.exonerate,"--model","protein2genome","--querytype","protein","--targettype","dna",
         "--query",str(output/"query.faa"),"--target",str(output/"target.fna"),"--bestn","2","--score","100",
         "--maxintron",str(args.max_intron),"--refine","region","--refineboundary","300","--dpmemory","256",
         "--showalignment","no","--showtargetgff","yes","--showvulgar","yes",
         "--ryo","CDS_BEGIN\n%tcs\nCDS_END\n"]
    (output/"command.json").write_text(json.dumps(cmd))
    began=time.monotonic()
    try:
        command(cmd,output/"exonerate.out",args.exonerate_timeout)
    except subprocess.TimeoutExpired:
        return [],dict(cid=cid,status="timeout",seconds=round(time.monotonic()-began,2))
    except Exception as error:
        return [],dict(cid=cid,status="error",error=str(error))
    new_models=[]
    for alignment in parse_output(output/"exonerate.out"):
        if alignment["strand"]!=row["strand"] or not alignment["exons"]: continue
        ref=refs[alignment["reference_id"]]
        model=build_model(cid,"exonerate",seq,new_offset,alignment["exons"],row["strand"],
                          frameshifts=alignment["frameshifts"],ref=ref.id,qstart=alignment["qstart"]+1,
                          qend=alignment["qend"],ref_length=len(ref.sequence),identity=alignment.get("identity",0.))
        if model["cds"]!=alignment["reported_cds"].upper() and not model["frameshifts"]:
            LOG.warning("Exonerate coordinate/CDS mismatch: %s",cid); continue
        model["alignment"]=dict(qstart=alignment["qstart"]+1,qend=alignment["qend"],qlen=len(ref.sequence))
        model["overlap"]=overlap_fraction(models[0]["exons"],model["exons"])
        if model["overlap"]<.5: continue
        new_models.append(model)
        for variant in terminal_variants(model,seq,new_offset,refs):
            variant["overlap"]=overlap_fraction(models[0]["exons"],variant["exons"]); new_models.append(variant)
    return new_models,dict(cid=cid,status="completed",seconds=round(time.monotonic()-began,2),model_count=len(new_models))


def export_results(rows, byid, decisions, args, stats, originals):
    final=args.outdir/"final"; final.mkdir(exist_ok=True)
    records=[]; gff=["##gff-version 3"]
    for row in rows:
        cid=row["candidate_id"]; d=decisions[cid]; m=d["model"]
        # Do not merge biological loci just because proteins match. Possible
        # same-locus duplications are flagged below, not silently collapsed.
        rec=dict(candidate_id=cid,sequence_id=row["sequence_id"],strand=row["strand"],family=m["family"],
                 original_category=row["original_category"],model_category=d["category"],confidence=d["confidence"],
                 original_length_aa=int(row["protein_length"]),final_length_aa=m["length"],model_source=m["source"],
                 selected_model_id=m["id"],start=min(s for s,e in m["exons"]),end=max(e for s,e in m["exons"]),
                 exon_count=len(m["exons"]),has_start=int(m["has_start"]),has_stop=int(m["has_stop"]),
                 frameshifts=m["frameshifts"],internal_stops=len(m["stops"]),cds_remainder=m["remainder"],ambiguous_bases=m["ambiguous"],
                 canonical_splicing=int(m["canonical"]),or_hmm_ga=int(m["or_ga"]),or_hmm_evalue=m["or_evalue"],
                 gr_hmm_ga=int(m["gr_ga"]),or_bitscore=m["or_bits"],gr_bitscore=m["gr_bits"],
                 compatible_full_reference_count=len(m["full_ref_support"]),
                 independent_full_reference_count=len(set(m["full_ref_support"])-set(m["reference_exact_matches"])),
                 reference_exact_matches=",".join(m["reference_exact_matches"]),
                 selected_cds_overlap=round(m["overlap"],4),stop_reference_count=d["stop_reference_count"],
                 stop_coordinates=d["stop_coordinates"],crosscheck_sources=d["crosscheck_sources"],reason=d["reason"],
                 same_locus_overlap_warning="")
        for previous in records:
            pm=decisions[previous["candidate_id"]]["model"]
            if previous["sequence_id"]==row["sequence_id"] and pm["strand"]==m["strand"]:
                if overlap_fraction(pm["exons"],m["exons"])>=.8 and overlap_fraction(m["exons"],pm["exons"])>=.8:
                    rec["same_locus_overlap_warning"]=previous["candidate_id"]
        records.append(rec)
        attrs=f"ID={cid};model_id={m['id']};model_category={d['category']};confidence={d['confidence']};source_model={m['source']}"
        gff.append("\t".join(map(str,[row["sequence_id"],"OR_model_review","mRNA",rec["start"],rec["end"],".",m["strand"],".",attrs])))
        cumulative=0
        for i,(s,e) in enumerate(m["exons"],1):
            gff.append("\t".join(map(str,[row["sequence_id"],"OR_model_review","CDS",s,e,".",m["strand"],(3-cumulative%3)%3,f"ID={cid}.CDS{i};Parent={cid}"])))
            cumulative+=e-s+1
        # Coordinates are normalized to actual translated CDS, with no invented
        # bases or frame corrections. A terminal stop is included in these CDSs.
        original=next(x for x in byid[cid] if x["source"]=="original")
        rec["model_changed"]=int(m["cds"]!=original["cds"])
    write_table(final/"OR_catalogue.tsv",records,list(records[0]))
    write_table(final/"OR_reclassification.tsv",records,list(records[0]))
    (final/"OR_models.gff3").write_text("\n".join(gff)+"\n")
    for kind in ["complete","partial","pseudogene"]:
        members=[r for r in records if r["model_category"]==kind]
        write_table(final/f"OR_{kind}.tsv",members,list(records[0]))
        for seqtype,extension in [("protein","proteins.faa"),("cds","cds.fna")]:
            write_fasta(final/f"OR_{kind}.{extension}",[FastaRecord(r["candidate_id"],f"category={kind} confidence={r['confidence']} model={r['selected_model_id']}",decisions[r["candidate_id"]]["model"][seqtype]) for r in members])
    reviewed=[r for r in records if r["model_category"]=="complete" and r["confidence"]=="reviewed_complete"]
    write_fasta(final/"OR_complete_reviewed.cds.fna",[FastaRecord(r["candidate_id"],"",decisions[r["candidate_id"]]["model"]["cds"]) for r in reviewed])
    original_audit=[]; model_audit=[]; stop_audit=[]
    for row in rows:
        cid=row["candidate_id"]; orig=next(m for m in byid[cid] if m["source"]=="original")
        original_audit.append(dict(candidate_id=cid,original_category=row["original_category"],protein_matches_original=int(orig["protein_match"]),
                                   length_aa=orig["length"],has_start=int(orig["has_start"]),has_stop=int(orig["has_stop"]),internal_stops=len(orig["stops"]),
                                   frameshifts=orig["frameshifts"],cds_remainder=orig["remainder"],ambiguous_bases=orig["ambiguous"],splice_sites=",".join(orig["sites"])))
        for model in byid[cid]:
            model_audit.append(dict(candidate_id=cid,model_id=model["id"],source=model["source"],reference=model["ref"],
                                    alignment_tool=model.get('alignment_tool') or ('exonerate' if 'exonerate' in model['source'] else 'miniprot'),
                                    alignment_reference_origin=model.get('alignment_reference_origin','unknown'),
                                    alignment_reference_qualified=int(model.get('alignment_reference_qualified',False)),
                                    reference_sequence_key=model.get('reference_sequence_key',''),
                                    reference_coverage=round(model["ref_cov"],4),length_aa=model["length"],has_start=int(model["has_start"]),has_stop=int(model["has_stop"]),
                                    internal_stops=len(model["stops"]),frameshifts=model["frameshifts"],canonical_splicing=int(model["canonical"]),
                                    cds_overlap=round(model["overlap"],4),bridge=int(model["bridge"]),family=model["family"],or_ga=int(model["or_ga"]),
                                    full_reference_support=",".join(model["full_ref_support"]),accepted_complete=int(complete_model(model))))
            for stop in model["stops"]:
                stop_audit.append(dict(candidate_id=cid,model_id=model["id"],source=model["source"],reference=model["ref"],**stop))
    write_table(args.outdir/"original_model_qc.tsv",original_audit,list(original_audit[0]))
    write_table(args.outdir/"alternative_model_evidence.tsv",model_audit,list(model_audit[0]))
    write_table(args.outdir/"internal_stop_evidence.tsv",stop_audit,list(stop_audit[0]) if stop_audit else ["candidate_id","model_id"])
    transitions=Counter(f"{r['original_category']}->{r['model_category']}" for r in records)
    stats.update(candidate_count=len(records),original_counts=dict(Counter(r["original_category"] for r in records)),
                 final_counts=dict(Counter(r["model_category"] for r in records)),confidence_counts=dict(Counter(r["confidence"] for r in records)),
                 transitions=dict(transitions),changed_models=sum(r["model_changed"] for r in records),
                 overlapping_final_models=sum(bool(r["same_locus_overlap_warning"]) for r in records),
                 family_counts=dict(Counter(r["family"] for r in records)),
                 original_protein_mismatches=sum(not r["protein_matches_original"] for r in original_audit))
    stats.setdefault("limitations", ["No RNA/read validation", "Not a genome-wide sensitivity benchmark", "Records not equivalent to independently confirmed biological genes", "No automatic same-locus merging in this diagnostic run"])
    (final/"OR_summary.json").write_text(json.dumps(stats,indent=2)+"\n")
    return stats



from __future__ import annotations
import argparse,json,re,sys,shutil
from pathlib import Path
from . import __version__

class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,argparse.RawDescriptionHelpFormatter):pass

def parser():
    p=argparse.ArgumentParser(prog='orsearcher',description='Insect OR genome annotation, bounded model repair and auditable locus consolidation.')
    p.add_argument('--version',action='version',version='ORsearcher '+__version__)
    commands=p.add_subparsers(dest='command',required=True)
    commands.add_parser('check',help='Check binaries and R/ape, without installing or changing software')
    d=commands.add_parser('references',help='Show default reference availability and provenance')
    d.add_argument('--reference-dir',type=Path)
    b=commands.add_parser('prepare-benchmark',help='Extract OR/GR annotation for held-out evaluation only')
    b.add_argument('--annotation',type=Path,required=True);b.add_argument('--outdir',type=Path,required=True)
    benchmark_inputs=[b]
    b=commands.add_parser('benchmark',help='Compare exported models with held-out annotation')
    b.add_argument('--annotation',type=Path,required=True);b.add_argument('--prediction',type=Path,required=True);b.add_argument('--outdir',type=Path,required=True)
    benchmark_inputs.append(b)
    b.add_argument('--sample',help='Select one sample when --prediction is a batch output directory')
    b.add_argument('--prediction-scope',choices=['auto','final','legacy','candidates'],default='auto',help='Auto prefers final results; explicitly select legacy or unverified discovery candidates')
    b.add_argument('--seqid-map',type=Path,help='Explicit TSV: prediction_seqid, annotation_seqid; never strip accession versions implicitly')
    b.add_argument('--min-reference-coverage',type=float,default=.5,help='Primary locus-recovery threshold; coverage tiers 50/75/90%% are also reported')
    b.add_argument('--min-predicted-support',type=float,default=.3,help='Minimum fraction of predicted coding bases supported by one annotated transcript')
    for benchmark_input in benchmark_inputs:
        benchmark_input.add_argument('--family-map',type=Path,help='Comparator-only TSV: feature_id, family (OR/GR/OTHER), optional name; explicit labels override automatic annotation names')
        benchmark_input.add_argument('--or-name-pattern',action='append',default=[],help='Additional case-sensitive regex for OR comparator identifiers; repeatable')
        benchmark_input.add_argument('--gr-name-pattern',action='append',default=[],help='Additional case-sensitive regex for GR comparator identifiers; repeatable')
    r=commands.add_parser('run',formatter_class=HelpFormatter,help='Fresh complete run; samples in a manifest are processed sequentially')
    r.epilog = """Examples:
  orsearcher --genome genome.fa
  orsearcher --genome genome.fa --threads 12 --preset G2 --outdir result_g2
  orsearcher --genome genome.fa --replace-or custom_proteins --add-gr extra_proteins
  orsearcher --genome genome.fa --add-hmm OR extra_profile
  orsearcher --genome genome.fa --rounds 2 --next-round --preset G3

Commands: run, check, references, benchmark, prepare-benchmark, refine, polish.
Run `orsearcher COMMAND -h` for command-specific help.
"""
    essentials = r.add_argument_group('Genome, CPU budget and preset')
    inputs = essentials.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--genome', type=Path, help='Genome nucleotide FASTA (plain text or BGZF-compressed)')
    essentials.add_argument('--threads', type=int, default=36, help='Maximum CPU budget for one genome')
    essentials.add_argument('--preset', type=str.upper, choices=('G0','G1','G2','G3'), default='G0',
        help='G0 baseline, G1 relaxed discovery, G2 relaxed identity/reconstruction, G3 relaxed completeness support. Explicit options override the preset')
    inputs.add_argument('--manifest', type=Path, help='Alternative TSV with sample, genome and optional annotation columns')
    essentials.add_argument('--outdir', type=Path, help='Output directory (default: orsearcher_<sample>)')
    essentials.add_argument('--sample', help='Sample ID (default: inferred from the genome filename)')
    refs = r.add_argument_group('Replace or add reference files for this run')
    refs.add_argument('--replace-hmm', nargs=2, action='append', metavar=('OR|GR','FILE'), default=[],
        help='Replace bundled profile(s) for the given family with amino-acid HMMER3 FILE. Repeat for multiple files or families')
    refs.add_argument('--replace-or', type=Path, metavar='FASTA', help='Replace bundled OR proteins and infer template eligibility')
    refs.add_argument('--replace-gr', type=Path, metavar='FASTA', help='Replace bundled GR proteins')
    refs.add_argument('--add-hmm', nargs=2, action='append', metavar=('OR|GR','FILE'), default=[],
        help='Add amino-acid HMMER3 FILE to the given family. Repeat for multiple files or families')
    refs.add_argument('--add-or', type=Path, action='append', metavar='FASTA', default=[], help='Add OR proteins (repeatable)')
    refs.add_argument('--add-gr', type=Path, action='append', metavar='FASTA', default=[], help='Add GR proteins (repeatable)')
    run = r.add_argument_group('Run control and advanced inputs')
    run.add_argument('--annotation', type=Path, help='Optional existing genome annotation in GFF3 format')
    run.add_argument('--rounds', type=int, default=1, help='Number of full workflow runs, each with its own discovery passes')
    run.add_argument('--next-round', action='store_true', help='Start the next-round override block. Other explicit options inherit')
    run.add_argument('--round-seeds', choices=['final-intact','off'], default='final-intact', help='Use intact OR models from the previous round as search-only seeds. Manifests require off')
    run.add_argument('--search-seeds', type=Path, help='Extra search-only protein FASTA. Use none to clear inherited seeds')
    run.add_argument('--or-reference', type=Path, help='Advanced replacement preserving FASTA IDs for reproducing previous analyses')
    run.add_argument('--gr-reference', type=Path, help='Advanced replacement preserving FASTA IDs for reproducing previous analyses')
    run.add_argument('--or-hmm', type=Path, help='Advanced replacement for OR profiles')
    run.add_argument('--gr-hmm', type=Path, nargs='+', help='Advanced replacement for all GR profile files')
    run.add_argument('--reference-dir', type=Path, help='Use a compatible external default reference bundle')
    run.add_argument('--complete-reference', type=Path, help='Explicit OR template subset with matching IDs and sequences. Use none for automatic qualification')
    run.add_argument('--no-background', action='store_true', help='Disable other-protein competition')
    run.add_argument('--exclude-reference-pattern', action='append', default=[], help='Exclude reference headers matching this regex (repeatable)')
    run.add_argument('--exclude-reference-proteins', type=Path, help='Exclude proteins matching sequences in this FASTA exactly')
    run.add_argument('--log-level', choices=['debug','info','warning'], default='info', help='Console logging level')
    run.add_argument('--quiet', action='store_true', help='Show warnings and errors on the console and keep full file logs')
    run.add_argument('--jobs', type=int, help='Concurrent genome shards (default: up to 4 within the CPU budget)')
    run.add_argument('--workers', type=int, help='Concurrent local alignments (default: up to 18 within the CPU budget)')
    run.add_argument('--resume', action='store_true', help='Resume identical inputs, settings and software')
    discovery = r.add_argument_group('Genome search and local alignment')
    discovery.add_argument('--chunk-size',type=int,default=20_000_000,help='Maximum genome window size in bases before adding overlaps')
    discovery.add_argument('--overlap',type=int,default=500_000,help='Window overlap; not a guarantee for loci longer than this')
    discovery.add_argument('--shard-size',type=int,default=1_000_000_000,help='Pack windows and small contigs into a limited number of files')
    discovery.add_argument('--max-intron',type=int,default=200_000,help='Maximum intron length in bases for protein-to-genome alignment')
    discovery.add_argument('--max-window-bp',type=int,default=500_000,help='Maximum local reconstruction window in bases')
    discovery.add_argument('--flank-bp',type=int,default=100_000,help='Requested local flank on each side in bases, subject to the window limit')
    discovery.add_argument('--max-hits',type=int,default=2000,help='Maximum Miniprot secondary alignments per query and shard')
    discovery.add_argument('--exonerate-timeout',type=int,default=60,help='Exonerate time limit per local command in seconds')
    discovery.add_argument('--exonerate-word-length',type=int,choices=[3,4,5,6],default=3,help='Sensitive local protein seeding; 6 approximates the older heuristic')
    discovery.add_argument('--exonerate-legacy-rescue',choices=['auto','off'],default='auto',help='One bounded word-length-6 retry where external Miniprot/Exonerate structure does not agree')
    discovery.add_argument('--exonerate-legacy-budget',type=int,default=500,help='Deterministic per-sample candidate allowance for complementary reconstruction; deferred cases are logged')
    discovery.add_argument('--local-timeout',type=int,default=120,help='Miniprot time limit per local command in seconds')
    translated = r.add_argument_group('Translated search with TBLASTN')
    translated.add_argument('--translated-search',choices=['auto','off'],default='auto',help='Independent whole-genome TBLASTN and bounded local reconstruction')
    translated.add_argument('--translated-evalue',type=float,default=1e-3,help='Maximum TBLASTN E-value for translated-search hits')
    translated.add_argument('--translated-timeout',type=int,default=7200,help='Maximum TBLASTN seconds per packed shard')
    translated.add_argument('--translated-max-regions',type=int,default=5000,help='Reconstruction budget; deferred anchors are explicitly exported')
    translated.add_argument('--translated-flank-bp',type=int,default=15000,help='Requested flank on each side of a translated-search anchor in bases')
    translated.add_argument('--translated-references',type=int,default=6,help='Maximum distinct full reference proteins per translated anchor')
    recall = r.add_argument_group('Exon-set search with MetaEuk and GeneWise')
    recall.add_argument('--recall-search',choices=['auto','off'],default='auto',help='Independent MetaEuk raw exon sets and bounded GeneWise rescue')
    recall.add_argument('--recall-only',action='store_true',help='Export unfiltered discovery candidates; defer OR identity, full-length and pseudogene classification')
    recall.add_argument('--recall-evalue',type=float,default=0.1,help='Loose combined MetaEuk exon-set E-value, NOT an OR confirmation threshold')
    recall.add_argument('--recall-sensitivity',type=float,default=5.7,help='MetaEuk/MMseqs sensitivity; 7.5 is substantially more expensive')
    recall.add_argument('--recall-min-fragment-aa',type=int,default=30,help='Global translated ORF minimum; shorter exons remain addressable by local GeneWise')
    recall.add_argument('--recall-search-timeout',type=int,default=7200,help='MetaEuk search seconds per packed shard')
    recall.add_argument('--recall-max-regions',type=int,default=500,help='Maximum GeneWise anchor regions per sample; other candidates and deferred anchors retained')
    recall.add_argument('--recall-exon-sets',type=int,default=20,help='MetaEuk exon sets per reference/contig/strand; raw calls retained before overlap reduction')
    recall.add_argument('--recall-references',type=int,default=4,help='Distinct external references per GeneWise region')
    recall.add_argument('--recall-complementary-regions',type=int,default=100,help='Additional sensitive Miniprot/Exonerate regions after incomplete GeneWise; 0 disables')
    recall.add_argument('--recall-timeout',type=int,default=20,help='GeneWise seconds per reference/region')
    recall.add_argument('--recall-window-bp',type=int,default=60000,help='Maximum MetaEuk/GeneWise reconstruction window in bases')
    recall.add_argument('--recall-flank-bp',type=int,default=5000,help='Requested flank on each side of a MetaEuk/GeneWise anchor in bases')
    recall.add_argument('--wise-config',type=Path,help='GeneWise config directory; infer from its conda prefix if absent')
    structure = r.add_argument_group('Structure reconstruction and family evidence')
    structure.add_argument('--annotation-mode',choices=['unified','legacy'],default='unified',help='Unified structure-first annotation; legacy preserves the older final routing')
    structure.add_argument('--structure-max-loci',type=int,default=500,help='Actual structural-review allowance; deferred groups remain unresolved, not deleted')
    structure.add_argument('--structure-flank-bp',type=int,default=20000,help='Requested flank on each side of a locus under structural review in bases')
    structure.add_argument('--structure-references',type=int,default=8,help='Distinct audited external proteins; best whole-protein support included first')
    structure.add_argument('--structure-rescue-budget',type=int,default=64,help='Bounded extra alternative-reference checks for truncated/ambiguous groups; 0 disables')
    structure.add_argument('--structure-rescue-references',type=int,default=3,help='Maximum distinct audited proteins per extra check')
    structure.add_argument('--structure-rescue-timeout',type=int,default=20,help='Seconds per extra check; no implicit second full timeout')
    structure.add_argument('--background-reference',type=Path,help='Audited reviewed insect non-OR/non-GR FASTA, protein-only competition')
    structure.add_argument('--phylogeny-timeout',type=int,default=900,help='Bound each alignment/tree/QC command; failure leaves identity evidence unresolved')
    structure.add_argument('--phylogeny-max-aa',type=int,default=1500,help='Long models retained but omitted from OR/GR diagnostic alignment')
    polish = r.add_argument_group('Final model repair')
    polish.add_argument('--polish-repair-budget',type=int,default=256,help='Additional fixed loci to repair; 0 keeps consolidation but disables repair')
    polish.add_argument('--polish-references',type=int,default=2,help='Distinct reference proteins per repair locus')
    polish.add_argument('--polish-timeout',type=int,default=180,help='Maximum seconds for each repair alignment command')
    polish.add_argument('--polish-flank-bp',type=int,default=80000,help='Requested repair flank, bounded by max-window-bp')
    from .tuning import add_arguments
    add_arguments(r)
    f=commands.add_parser('refine',parents=[r],add_help=False,help='New structure/identity analysis from cached model JSONs, without re-searching genomes')
    f.add_argument('--candidate-models',type=Path,action='append',required=True,help='Annotation-blind model list JSON; repeat for stages; never a target GFF')
    f.add_argument('--evidence-cache',type=Path,help='Read-only compatible protein SQLite cache copied to the new output')
    q=commands.add_parser('polish',parents=[r],add_help=False,help='Consolidate and repair a previous reviewed sample into a NEW run')
    q.add_argument('--previous-sample',type=Path,required=True,help='Previous samples/<sample> directory containing structure/all_models.json')
    q.add_argument('--allow-unverified-previous',action='store_true',help='Explicitly accept missing/changed previous-genome metadata (not proof of equivalence; prefer refine)')
    return p

def main(argv=None):
    p=parser()
    argv=list(sys.argv[1:] if argv is None else argv)
    if not argv or (argv[0].startswith('-') and argv[0] != '--version'):
        argv = ['run', *(argv or ['--help'])]
    first=argv[:argv.index('--next-round')] if '--next-round' in argv else argv
    args=p.parse_args(first)
    if args.command in ('run','refine','polish'):
        from .presets import apply_preset, set_run_defaults
        apply_preset(args, first[1:], p._subparsers._group_actions[0].choices[args.command])
        set_run_defaults(args)
    if getattr(args,'rounds',1)>1 or '--next-round' in argv:
        try:
            from .workflow import dispatch
            return dispatch(argv,p,main)
        except (ValueError,OSError) as error:
            print('ORsearcher failed: '+str(error),file=sys.stderr);return 1
    from .engine import dependency_check,run,load_samples
    try:
        if args.command=='check':print(json.dumps(dependency_check(args),indent=2));return 0
        if args.command=='references':
            from .defaults import inventory
            print(json.dumps(inventory(args.reference_dir),indent=2));return 0
        if args.command in ['benchmark','prepare-benchmark']:
            from .benchmark import evaluate,prepare
            print(json.dumps(evaluate(args) if args.command=='benchmark' else str(prepare(args)),indent=2));return 0
        from .defaults import resolve_references
        from .tuning import validate
        args.infer_templates = args.complete_reference is not None and str(args.complete_reference).lower() == 'none'
        for name in ['search_seeds','complete_reference']:
            if str(getattr(args,name,None)).lower()=='none':setattr(args,name,None)
        validate(args)
        if args.rounds<1:raise ValueError('--rounds must be positive')
        if args.next_round:raise ValueError('--next-round requires multi-round run')
        if getattr(args,'search_seeds',None) and not args.search_seeds.is_file():raise ValueError('Missing --search-seeds FASTA')
        if args.command!='run' and (args.rounds!=1 or args.search_seeds):raise ValueError('Whole-pipeline rounds and search seeds require run')
        resolve_references(args)
        for field in ['threads','jobs','workers','chunk_size','shard_size','max_intron','max_window_bp','flank_bp','max_hits','exonerate_timeout','local_timeout','translated_timeout','translated_max_regions','translated_flank_bp','translated_references','exonerate_legacy_budget','recall_search_timeout','recall_max_regions','recall_exon_sets','recall_references','recall_timeout','recall_window_bp','recall_flank_bp','polish_references','polish_timeout','polish_flank_bp']:
            if getattr(args,field)<=0:p.error('--'+field.replace('_','-')+' must be positive')
        if args.jobs>args.threads or args.workers>args.threads:p.error('jobs/workers must not exceed threads')
        if not 0<=args.overlap<args.chunk_size or args.shard_size<args.chunk_size:p.error('Require 0 <= overlap < chunk-size <= shard-size')
        if args.overlap<args.max_intron:p.error('overlap should be at least max-intron; increase it for long gene spans')
        if not 0<args.translated_evalue<=10:p.error('translated-evalue must be in (0,10]')
        if args.translated_flank_bp*2>=args.max_window_bp:p.error('translated-flank-bp too large for max-window-bp')
        if not 0<args.recall_evalue<=100:p.error('recall-evalue must be in (0,100]')
        if not 1<=args.recall_sensitivity<=7.5:p.error('recall-sensitivity must be in [1,7.5]')
        if args.recall_min_fragment_aa<10:p.error('recall-min-fragment-aa must be >=10')
        if args.recall_complementary_regions<0:p.error('recall-complementary-regions must be >=0')
        for field in ['structure_max_loci','structure_flank_bp','structure_references','phylogeny_timeout','phylogeny_max_aa']:
            if getattr(args,field)<=0:p.error('--'+field.replace('_','-')+' must be positive')
        if args.structure_rescue_budget<0:p.error('--structure-rescue-budget must be >=0')
        if args.polish_repair_budget<0:p.error('--polish-repair-budget must be >=0')
        if args.structure_rescue_references<=0 or args.structure_rescue_timeout<=0:p.error('structure-rescue references/timeout must be positive')
        if args.command=='refine' and (args.manifest or args.recall_only):p.error('refine requires one --genome and cannot use --recall-only')
        if args.command=='polish' and (args.manifest or args.recall_only or args.annotation):p.error('polish requires one --genome; no manifest, annotation or recall-only')
        if args.command in ['run','refine'] and args.annotation_mode!='unified':p.error('The release workflow requires --annotation-mode unified; legacy results can be imported with polish')
        if args.recall_flank_bp*2>=args.recall_window_bp:p.error('recall-flank-bp too large for recall-window-bp')
        if args.recall_only and args.recall_search=='off':p.error('recall-only requires recall-search auto')
        if args.recall_search!='off' and args.wise_config is None:
            binary=shutil.which('genewise')
            if binary:args.wise_config=Path(binary).resolve().parents[1]/'share/wise2/wisecfg'
        samples=load_samples(args)
        if not samples or len({s.get('sample') for s in samples})!=len(samples):p.error('Manifest is empty or contains duplicate sample IDs')
        for row in samples:
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',row.get('sample','')) or row['sample'] in ['.','..']:p.error('sample must be a safe nonempty ASCII identifier')
            if not row.get('genome'):p.error('Manifest requires a genome path for each sample')
        for pattern in args.exclude_reference_pattern:re.compile(pattern)
        for path in [args.or_reference,args.gr_reference,args.or_hmm,*args.gr_hmm,*([args.complete_reference] if args.complete_reference else []),*([args.background_reference] if args.background_reference else []),*([args.exclude_reference_proteins] if args.exclude_reference_proteins else []),*[Path(s['genome']) for s in samples],*[Path(s['annotation']) for s in samples if s.get('annotation')]]:
            if not path.is_file() or not path.stat().st_size:p.error('Missing/empty input: '+str(path))
        if args.command=='polish':
            args.previous_sample=args.previous_sample.resolve()
            for name in ['structure/all_models.json','annotation_qc/identity_audit.tsv','nonredundant/OR_catalogue.tsv']:
                if not (args.previous_sample/name).is_file():p.error('Previous sample lacks '+name)
            target=args.outdir.resolve()
            if target==args.previous_sample or args.previous_sample in target.parents or target in args.previous_sample.parents:p.error('polish output must be separate from the previous sample tree')
            from .provenance import verify_previous_genome
            args.previous_genome_provenance=verify_previous_genome(args.previous_sample,args.genome,args.allow_unverified_previous)
        return run(args)
    except Exception as error:
        print(f'ORsearcher failed: {error}',file=sys.stderr);return 1


# Computational workflow and reporting limits

This document describes the 2.3 release-candidate workflow. It is a software
protocol, not a claim of benchmarked insect-wide sensitivity or specificity.

## Candidate discovery

The input genome is streamed into overlapping windows and packed shards,
reducing boundary losses without creating one file for every short contig.
Coordinates are restored to the original assembly. Two protein-guided search
passes use external OR references and a bounded set of first-pass/optional
annotation-derived seeds. A second pass still runs if no adaptive seeds pass.

Independent TBLASTN and MetaEuk exon-set searches can add discovery anchors;
bounded GeneWise/Miniprot/Exonerate reconstruction supplies candidate models.
HSPs are localization evidence, not directly exported CDS or pseudogene
damage evidence. Broad discovery records survive subsequent filtering.

## Reference qualification and family identity

Reference proteins are audited for domain support and reconstruction
eligibility. Automatic gap-fill references can be masked outside supported OR
domains for discovery; original proteins remain available for provenance.
Fragment references can localize a locus without proving a full-length model.

OR identity combines protein homology, OR/GR profile evidence, curated/distinct
reference support, diagnostic tree evidence and optional other-protein
competition. A single profile hit, length threshold or tree-side assignment
does not establish identity. Strong conflicting evidence remains unresolved
rather than being silently forced into OR or GR. Reference-tree votes require
alignment/reference QC and are not experimental validation.

## Structural review

Local models are checked for start/stop codons, reading frame, internal stops,
ambiguous sequence, splice structure, repeated-domain/fusion warnings and
reference-end compatibility. Independent tools or distinct reference proteins
provide different kinds of support; repeated parameters/aliases of the same
tool and sequence are not independent biological evidence. All reconstruction
uses assembly bases; reference-guided sequence grafting is prohibited.

## Practical polishing in 2.3

After structural/family review, fixed coding anchors consolidate compatible
alternatives at the same location. Coding overlap, frame and splice structure
are used instead of genomic-span proximity or peptide identity alone. A
bridge-like model must not transitively merge adjacent genes. Competing paths
are held with a recorded reason instead of counted as independent main genes.

Incomplete models and models whose complementary check previously timed out
can receive one bounded reconstruction attempt using selected qualified
references and a less expensive local alignment strategy. Repair budgets,
reference count, flank and timeout are explicit. Exhausting a budget or time
limit leaves incomplete evidence; it does not prove absence or pseudogenization.

Final representatives undergo OR/GR/background evidence checks and a
genome-to-CDS-to-protein round-trip before export. Complete practical models do
not require absolute two-tool unanimity. This is an explicit usability choice,
not a calibrated probability of correctness. Assembly-supported disrupted
candidates retain their caveats and links to intact overlapping alternatives.

## Minimum reporting for a study

Record the exact software version, genome accession/checksum, reference
manifest and overrides, native-tool versions, genetic-code assumption,
resource settings and stage budgets. Report independent loci separately from
alternative models. Include the counts of incomplete, unresolved and deferred
cases rather than reporting only the complete set.

If target-species proteins were in the search references, their recovery is
not a held-out accuracy estimate. Cross-order validation requires appropriate
independent truth/evidence and taxon-aware benchmarking. No read/RNA-based
damage confirmation or functional assay is part of this workflow.

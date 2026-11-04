# Output contract

## Primary model sets

Each sample's `results/` directory has three model classes. Each class has a
protein FASTA, genome-derived CDS FASTA and GFF3 with matched model identifiers.
Empty classes still have empty FASTA files and a valid GFF3 header.

| Class | Interpretation | Important exclusions |
| --- | --- | --- |
| `complete` | Practical intact model with start/stop, acceptable splice/frame structure, full-length-qualified homology and OR identity support | Not experimental full-length/function confirmation; not necessarily reviewed by two independent aligners |
| `partial` | Credible OR sequence but full-length model requirements are unmet | Not automatically a pseudogene; may be recoverable with better evidence or assembly |
| `pseudogene` | Retained assembly-supported disruptive model/candidate from structural review | Not inferred from timeout, shortness, missing termini or failed alignment alone; assembly error remains possible |

FASTA descriptions retain classification and source-model information. GFF3
coordinates are 1-based inclusive and refer to the original genome, not shard
coordinates. CDS phases follow transcriptional order. Protein/CDS sequences
are reconstructed from the assembly rather than copied from references.

## Summary and provenance

`OR_catalogue.tsv` is the primary human-readable summary. Treat its header as
the authoritative schema for the installed version. Information includes
model ID, contig/strand/coordinates, class and confidence, protein length,
source-model/group IDs, reference support and coverage, start/stop and damage
indicators, and alternative/overlap disposition where applicable.

The catalogue additionally records OR/GR homology scores, profile support,
other-protein competition, canonical splicing, repair provenance and whether
an earlier strict complete review exists. `independent_locus=False` links
damaged alternatives to intact or previously selected damaged anchors; those
rows must not increase the independent gene count.

In standalone `polish`, a retained pseudogene's `damage_evidence` is explicitly
`inherited_assembly_supported_not_revalidated_by_polish`. Protein identity is
recomputed under active references, but the original damage conclusion is
not a new independent validation. Use a fresh `refine`/`run` for a structural
reassessment, and reads/RNA or experiments for biological confirmation.

`OR_models.gff3`, when present, is the combined coordinate catalogue; use the
class-specific GFF3 files for direct category-specific workflows.

`summary.json` records counts and completion warnings. `export_validation.json`
records coordinate-to-CDS and CDS-to-protein round-trip checks. A successful
round-trip establishes internal consistency, not biological model accuracy.

## Audit outputs are part of the result

Polishing preserves model assignments, competing coding-overlap paths,
identity-held models, selected/deferred repair work and command outcomes.
Earlier `recall_candidates/`, `structure/`, `annotation_qc/`,
`identity_unresolved/`, `excluded_non_OR/` and `unreviewed_candidates/`
directories retain upstream evidence when produced. Their records are not
additional independent genes in the final count.

Do not add category counts blindly if a disrupted alternative is linked to an
intact representative at the same locus. Use the independence/overlap
disposition in the catalogue/audit. A practical main-set choice does not mean
every held alternative is a false gene or GR.

## Logging and completion

At batch level, `batch.log` contains concise operational information;
`debug.log` preserves detailed program diagnostics. Each sample has its own
run log, stage status/summary files and local tool command/output files.
Input manifests and checksums make reference/configuration choices traceable.

`completed` means the configured workflow finished. `completed_with_warnings`
means outputs were written but limitations such as deferred work, failed
reconstruction or unresolved competition remain. Neither status means all OR
genes in the species have been found. `error`/`failed` requires inspecting the
reported stage and tool diagnostics.

For downstream selection-pressure analysis, start from complete models,
establish orthology, and inspect protein-guided codon alignments before testing.
Do not use category names as substitutes for those checks.

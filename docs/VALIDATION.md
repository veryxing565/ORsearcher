# Release and scoring validation

## 2.3.0rc2 — scorer compatibility, 2026-10-29

- Native Ubuntu/Python 3.13: **283 regression tests passed**, including R/ape.
- Final sample-directory and direct `results/` inputs both recovered the two
  expected models in the existing synthetic fixture. Its comparator was
  self-derived solely for I/O checking, not independent biological validation.
- Real whole-annotation parsing passed for Drosophila melanogaster, Bombyx
  mori and Dendroctonus ponderosae. This tests names/GFF compatibility, not
  three new genome-prediction runs or complete curated OR/GR inventories.
- A native test exposed legitimate multipart trans-splicing of an unrelated
  Drosophila gene. The reader was corrected and regression tests added; the
  failed pre-correction QA directory was retained rather than rewritten.
- An existing candidate10 standard Drosophila result was rescored read-only.
  It matched 60/60 automatically recognized coding-OR comparator loci at 50%,
  59 at 75%, 58 at 90%, and 14 strict complete coordinate/phase models. These
  are results of that **older cached standard run**, not accuracy gains from
  rc2, not a held-out test, and not proof the annotation lists every true OR.
- Source payload comparison confirmed only the scorer, CLI and version file
  changed under `src/`; predictor code and reference bytes were unchanged.
- No new whole-genome prediction/calibration runs were started. This supports
  entering small multi-order development calibration, not freezing for a
  publication-level accuracy comparison. See [benchmarking](BENCHMARK.md).

The following is the historical rc1 integration validation, not a claim that
its whole-genome smoke was rerun for rc2.

## 2.3.0rc1 — integration validation

Validated on 2026-10-28. These are software integration checks, **not** an
estimate of gene-annotation sensitivity, specificity or functional accuracy.

## Native Linux checks

- Existing isolated analysis environment, Python 3.11, with the executable
  dependencies reported by `orsearcher check`.
- **236 regression tests passed**, including native R/ape tree-QC cases.
- The existing deployed candidate10 CLI/engine/refinement source fingerprints
  matched the development baseline before integration; old code/results were
  not overwritten.

## End-to-end synthetic fixture

`tools/make_smoke_genome.py` plants a reference-derived OR with a canonical
synthetic intron on each strand, plus a random-background contig. Because the
planted sequence is an exact reference, recovery is deliberately **not** an
independent discovery or cross-species benchmark.

The final full `run` used 12 threads, two shard jobs and six local workers,
with small declared candidate/reconstruction budgets. It completed in
**97.85 seconds** and produced:

| Check | Observed result |
| --- | --- |
| Independent complete models | 2: one plus-strand, one minus-strand |
| Partial / pseudogene records | 0 / 0 |
| Required three-category protein/CDS/GFF3 files | All present, including empty classes |
| Genome/GFF/CDS/protein mismatches | 0 |
| Complete coding-format failures | 0 |
| Competing alternatives | Retained in the audit, not counted as independent genes |

The status was `completed_with_warnings`: the intentionally small GeneWise
budget deferred anchors, some anchors remained unresolved, and coding-overlap
alternatives were held. These warnings now propagate to the sample/batch
summary instead of being hidden by successful final export.

## Previously timed-out repair scheduling

A **copy** of the synthetic reviewed sample was used; only that copy's old
crosscheck statuses were set to `timeout`. This is a simulated prior-failure
fixture, not an observed real failure and not evidence that a biological
model was previously wrong.

Standalone `polish` then selected two fixed loci. Both native repair tasks
completed, producing 12 alternative alignment models. Final output remained
two independent complete loci, with zero coordinate/sequence mismatches.
Successful alignment records did not inflate gene counts. Original test
outputs and all earlier real-genome analyses remained untouched.

## Distribution checks

- Standard wheel and source-distribution builds succeeded.
- Clean-target installation and imports/CLI checks succeeded.
- Public artifacts include reference provenance/checksums and exclude the
  restricted derived `GR.faa`; local reference preparation reproduces its
  pinned checksum from the official source archive.
- The public source ZIP omits caches, analysis outputs, genomes and private
  reference-preparation scratch files.

## What has not been established

This release has **not** had a fresh 11-GB real-genome run or a multi-order
biological benchmark. The practical length/homology rules are not calibrated
across insects. Assembly-supported disruptions are not read/RNA-validated
pseudogenes. The release-candidate label is intentional.

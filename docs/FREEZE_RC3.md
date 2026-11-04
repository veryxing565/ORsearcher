# 2.3.0rc3 candidate — freeze contract, 2026-10-31

This is a software release candidate, not a frozen release. Scope: insect genome OR discovery, reconstruction, OR/GR/background discrimination, bounded repair, three output categories, batch operation and reproducible audit. No paper preparation, new reference payloads or unrelated biological analyses.

## Required gates

1. Native regression suite passes. Genome/GFF/CDS/protein and phase checks have zero mismatches for all eight final outputs and both large Locusta candidate pools. No candidate may be skipped to obtain a speedup.
2. Fresh eight-run regression uses unchanged genomes, reference pools, prediction thresholds and budgets. One genome at a time; 36 compute threads maximum, 4 shards, 18 single-thread reconstruction workers; CPU affinity restricted to 36 logical CPUs. All prediction exits must be zero; warnings and missing scores remain explicit.
3. For each of six scorable arms, 50/75/90% one-to-one recovery and strict complete structure counts must not decrease versus rc2. GR overlaps, duplicate/fusion indicators must not increase without individual adjudication. A count increase alone is not evidence of improved accuracy.
4. Coordinate/phase identity is reported separately from complete-model identity against a non-partial reference. Retain the historical strict metric unchanged. Partial models expose the unmet completeness requirements.
5. Indexed extraction must match samtools exactly on randomized wrapped/CRLF and real genome intervals; eight final output validations and the two million-scale candidate validations must pass. Real-interval extraction median performance target: at least 2x over three repeat timings; whole-pipeline time must be reported separately. Full Locusta runtime and memory must not exceed the rc2 baseline by more than 10% without explanation and review.
6. Record source/input hashes, actual dependency binaries and versions, environment export, commands, elapsed time, memory, timeout/unresolved counts. Package must install offline using available build dependencies and pass CLI smoke checks; do not distribute restricted GR data.

Engineering candidate acceptance requires all applicable automated gates plus manual review of unresolved warnings. General insect annotation quality remains blocked until the poor held-out-species recovery and missing Locusta comparator are addressed. Do not lower standards or relabel uncertain models merely to pass a gate. Any frozen release must explicitly state the supported scope and remaining reference dependence; no general accuracy claim follows from this calibration.

## Changes

- Indexed plain-FASTA reads and bounded 8 MiB interval cache replace per-model samtools spawning in export validation. Compressed/missing-index inputs retain samtools fallback. All models retain CDS, protein, phase and completeness checks; progress every 10,000 models.
- Repair planning examines all retained alternatives when finding qualified templates instead of silently stopping after eight. Existing per-locus reference limits, window bounds, global repair budget and identity thresholds remain intact.
- Final catalogue adds `completeness_blockers`; benchmark adds coordinate/phase agreement counts without changing strict structure or recall definitions.

No whole-genome improvement is established until the background evaluation finishes. No automatic promotion to a final release.

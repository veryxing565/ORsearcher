# Command-line reference

Run `orsearcher -h` for the current installed options. The tables below list
G0 defaults. Lengths are in bases or amino acids as specified. Coverage
thresholds use fractions from 0 to 1, and amino-acid identity uses percent.

## Genome, CPU budget and preset

| Option | Meaning |
| --- | --- |
| `--genome` | Genome nucleotide FASTA (plain text or BGZF-compressed) |
| `--threads` | Maximum CPU budget for one genome (default: 36) |
| `--preset` | G0 baseline, G1 relaxed discovery, G2 relaxed identity/reconstruction, G3 relaxed completeness support. Explicit options override the preset (default: G0) Choices: G0, G1, G2, G3. |
| `--manifest` | Alternative TSV with sample, genome and optional annotation columns |
| `--outdir` | Output directory (default: `orsearcher_<sample>`) |
| `--sample` | Sample ID (default: inferred from the genome filename) |

## Replace or add reference files for this run

| Option | Meaning |
| --- | --- |
| `--replace-hmm OR FILE`, `--replace-hmm GR FILE` | Replace bundled profile(s) for the given family with amino-acid HMMER3 FILE. Repeat for multiple files or families |
| `--replace-or` | Replace bundled OR proteins and infer template eligibility |
| `--replace-gr` | Replace bundled GR proteins |
| `--add-hmm OR FILE`, `--add-hmm GR FILE` | Add amino-acid HMMER3 FILE to the given family. Repeat for multiple files or families |
| `--add-or` | Add OR proteins (repeatable) |
| `--add-gr` | Add GR proteins (repeatable) |

## Run control and advanced inputs

| Option | Meaning |
| --- | --- |
| `--annotation` | Optional existing genome annotation in GFF3 format |
| `--rounds` | Number of full workflow runs, each with its own discovery passes (default: 1) |
| `--next-round` | Start the next-round override block. Other explicit options inherit (default: False) |
| `--round-seeds` | Use intact OR models from the previous round as search-only seeds. Manifests require off (default: final-intact) Choices: final-intact, off. |
| `--search-seeds` | Extra search-only protein FASTA. Use none to clear inherited seeds |
| `--or-reference` | Advanced replacement preserving FASTA IDs for reproducing previous analyses |
| `--gr-reference` | Advanced replacement preserving FASTA IDs for reproducing previous analyses |
| `--or-hmm` | Advanced replacement for OR profiles |
| `--gr-hmm` | Advanced replacement for all GR profile files |
| `--reference-dir` | Use a compatible external default reference bundle |
| `--complete-reference` | Explicit OR template subset with matching IDs and sequences. Use none for automatic qualification |
| `--no-background` | Disable other-protein competition (default: False) |
| `--exclude-reference-pattern` | Exclude reference headers matching this regex (repeatable) |
| `--exclude-reference-proteins` | Exclude proteins matching sequences in this FASTA exactly |
| `--log-level` | Console logging level (default: info) Choices: debug, info, warning. |
| `--quiet` | Show warnings and errors on the console and keep full file logs (default: False) |
| `--jobs` | Concurrent genome shards (default: up to 4 within the CPU budget) |
| `--workers` | Concurrent local alignments (default: up to 18 within the CPU budget) |
| `--resume` | Resume identical inputs, settings and software (default: False) |

## Genome search and local alignment

| Option | Meaning |
| --- | --- |
| `--chunk-size` | Maximum genome window size in bases before adding overlaps (default: 20000000) |
| `--overlap` | Window overlap; not a guarantee for loci longer than this (default: 500000) |
| `--shard-size` | Pack windows and small contigs into a limited number of files (default: 1000000000) |
| `--max-intron` | Maximum intron length in bases for protein-to-genome alignment (default: 200000) |
| `--max-window-bp` | Maximum local reconstruction window in bases (default: 500000) |
| `--flank-bp` | Requested local flank on each side in bases, subject to the window limit (default: 100000) |
| `--max-hits` | Maximum Miniprot secondary alignments per query and shard (default: 2000) |
| `--exonerate-timeout` | Exonerate time limit per local command in seconds (default: 60) |
| `--exonerate-word-length` | Sensitive local protein seeding; 6 approximates the older heuristic (default: 3) Choices: 3, 4, 5, 6. |
| `--exonerate-legacy-rescue` | One bounded word-length-6 retry where external Miniprot/Exonerate structure does not agree (default: auto) Choices: auto, off. |
| `--exonerate-legacy-budget` | Deterministic per-sample candidate allowance for complementary reconstruction; deferred cases are logged (default: 500) |
| `--local-timeout` | Miniprot time limit per local command in seconds (default: 120) |

## Translated search with TBLASTN

| Option | Meaning |
| --- | --- |
| `--translated-search` | Independent whole-genome TBLASTN and bounded local reconstruction (default: auto) Choices: auto, off. |
| `--translated-evalue` | Maximum TBLASTN E-value for translated-search hits (default: 0.001) |
| `--translated-timeout` | Maximum TBLASTN seconds per packed shard (default: 7200) |
| `--translated-max-regions` | Reconstruction budget; deferred anchors are explicitly exported (default: 5000) |
| `--translated-flank-bp` | Requested flank on each side of a translated-search anchor in bases (default: 15000) |
| `--translated-references` | Maximum distinct full reference proteins per translated anchor (default: 6) |

## Exon-set search with MetaEuk and GeneWise

| Option | Meaning |
| --- | --- |
| `--recall-search` | Independent MetaEuk raw exon sets and bounded GeneWise rescue (default: auto) Choices: auto, off. |
| `--recall-only` | Export unfiltered discovery candidates; defer OR identity, full-length and pseudogene classification (default: False) |
| `--recall-evalue` | Loose combined MetaEuk exon-set E-value, NOT an OR confirmation threshold (default: 0.1) |
| `--recall-sensitivity` | MetaEuk/MMseqs sensitivity; 7.5 is substantially more expensive (default: 5.7) |
| `--recall-min-fragment-aa` | Global translated ORF minimum; shorter exons remain addressable by local GeneWise (default: 30) |
| `--recall-search-timeout` | MetaEuk search seconds per packed shard (default: 7200) |
| `--recall-max-regions` | Maximum GeneWise anchor regions per sample; other candidates and deferred anchors retained (default: 500) |
| `--recall-exon-sets` | MetaEuk exon sets per reference/contig/strand; raw calls retained before overlap reduction (default: 20) |
| `--recall-references` | Distinct external references per GeneWise region (default: 4) |
| `--recall-complementary-regions` | Additional sensitive Miniprot/Exonerate regions after incomplete GeneWise; 0 disables (default: 100) |
| `--recall-timeout` | GeneWise seconds per reference/region (default: 20) |
| `--recall-window-bp` | Maximum MetaEuk/GeneWise reconstruction window in bases (default: 60000) |
| `--recall-flank-bp` | Requested flank on each side of a MetaEuk/GeneWise anchor in bases (default: 5000) |
| `--wise-config` | GeneWise config directory; infer from its conda prefix if absent |

## Structure reconstruction and family evidence

| Option | Meaning |
| --- | --- |
| `--annotation-mode` | Unified structure-first annotation; legacy preserves the older final routing (default: unified) Choices: unified, legacy. |
| `--structure-max-loci` | Actual structural-review allowance; deferred groups remain unresolved, not deleted (default: 500) |
| `--structure-flank-bp` | Requested flank on each side of a locus under structural review in bases (default: 20000) |
| `--structure-references` | Distinct audited external proteins; best whole-protein support included first (default: 8) |
| `--structure-rescue-budget` | Bounded extra alternative-reference checks for truncated/ambiguous groups; 0 disables (default: 64) |
| `--structure-rescue-references` | Maximum distinct audited proteins per extra check (default: 3) |
| `--structure-rescue-timeout` | Seconds per extra check; no implicit second full timeout (default: 20) |
| `--background-reference` | Audited reviewed insect non-OR/non-GR FASTA, protein-only competition |
| `--phylogeny-timeout` | Bound each alignment/tree/QC command; failure leaves identity evidence unresolved (default: 900) |
| `--phylogeny-max-aa` | Long models retained but omitted from OR/GR diagnostic alignment (default: 1500) |

## Final model repair

| Option | Meaning |
| --- | --- |
| `--polish-repair-budget` | Additional fixed loci to repair; 0 keeps consolidation but disables repair (default: 256) |
| `--polish-references` | Distinct reference proteins per repair locus (default: 2) |
| `--polish-timeout` | Maximum seconds for each repair alignment command (default: 180) |
| `--polish-flank-bp` | Requested repair flank, bounded by max-window-bp (default: 80000) |

## Discovery and adaptive seed thresholds

| Option | Meaning |
| --- | --- |
| `--adaptive-min-aa` | Minimum adaptive search-seed protein length (default: 250) |
| `--adaptive-max-aa` | Maximum adaptive search-seed protein length (default: 650) |
| `--adaptive-reference-coverage` | Minimum adaptive-seed reference coverage fraction (default: 0.7) |
| `--adaptive-min-identity` | Minimum adaptive-seed identity percent (default: 20.0) |
| `--first-seed-hits` | Minimum Miniprot seed hits in the first discovery pass (default: 6) |
| `--second-seed-hits` | Minimum Miniprot seed hits in the relaxed discovery pass (default: 3) |
| `--first-output-coverage` | First-pass Miniprot --outc coverage fraction (default: 0.25) |
| `--second-output-coverage` | Second-pass Miniprot --outc coverage fraction (default: 0.15) |
| `--first-alignment-ratio` | First-pass Miniprot -p secondary alignment ratio (default: 0.3) |
| `--second-alignment-ratio` | Second-pass Miniprot -p secondary alignment ratio (default: 0.2) |
| `--alignment-output-score` | Miniprot --outs relative output score (default: 0.2) |
| `--candidate-min-aa` | Minimum protein length admitted to structural reconstruction (default: 80) |
| `--candidate-profile-evalue` | Maximum OR-profile E-value for structural admission (default: 0.1) |
| `--candidate-bits` | Minimum homolog bitscore for structural admission (default: 35.0) |
| `--candidate-evalue` | Maximum homolog E-value for structural admission (default: 0.001) |
| `--candidate-query-coverage` | Minimum query coverage for structural admission (default: 0.35) |
| `--adaptive-require-full-reference`, `--no-adaptive-require-full-reference` | Require qualified full-reference support for internal adaptive search seeds (default: True) |

## Family identity and protein evidence thresholds

| Option | Meaning |
| --- | --- |
| `--evidence-or-targets` | Maximum retained external OR DIAMOND targets per protein (default: 50) |
| `--evidence-gr-targets` | Maximum retained external GR DIAMOND targets per protein (default: 5) |
| `--evidence-evalue` | DIAMOND protein-evidence search E-value; must encompass identity thresholds (default: 0.001) |
| `--evidence-domain-overlap` | Minimum reference OR-domain overlap for non-curated homology evidence (default: 0.25) |
| `--identity-evalue` | Maximum E-value for credible OR and GR homologs (same criteria) (default: 1e-05) |
| `--identity-bits` | Minimum credible OR and GR homolog bitscore (default: 50.0) |
| `--identity-query-coverage` | Minimum credible OR and GR query coverage (default: 0.5) |
| `--identity-min-identity` | Minimum credible OR and GR amino-acid identity percent (default: 25.0) |
| `--distant-evalue` | Maximum E-value for distant OR and GR homologs (default: 0.001) |
| `--distant-bits` | Minimum distant OR and GR homolog bitscore (default: 35.0) |
| `--distant-query-coverage` | Minimum distant OR and GR query coverage (default: 0.35) |
| `--distant-min-identity` | Minimum distant OR and GR identity percent (default: 15.0) |
| `--profile-evalue` | Maximum profile domain E-value for OR/GR support (default: 0.001) |
| `--profile-coverage` | Minimum OR/GR profile coverage fraction (default: 0.25) |
| `--family-score-margin` | Required OR-minus-GR or GR-minus-OR bitscore margin (default: 10.0) |
| `--distant-score-margin` | Required OR-minus-GR margin for the distant-OR path (default: 15.0) |
| `--background-score-margin` | Required background-minus-best-family bitscore margin (default: 30.0) |
| `--identity-use-tree`, `--no-identity-use-tree` | Retain QC-approved diagnostic-tree votes during final identity assessment; default: disabled to preserve rc3 behavior. A tree alone does not establish OR identity. (default: False) |

## Template selection and model assessment thresholds

| Option | Meaning |
| --- | --- |
| `--complete-min-aa` | Minimum protein length for complete coding models (default: 250) |
| `--complete-max-aa` | Maximum protein length for complete coding models (default: 650) |
| `--partial-min-aa` | Minimum accepted partial-model protein length (default: 150) |
| `--partial-max-aa` | Maximum accepted partial-model protein length (default: 650) |
| `--model-hit-evalue` | Maximum homology E-value for final complete/partial model assessment (default: 1e-20) |
| `--model-hit-bits` | Minimum homology bitscore for final model assessment (default: 100.0) |
| `--model-hit-identity` | Minimum homology identity percent for final model assessment (default: 30.0) |
| `--model-hit-coverage` | Minimum homology query coverage for final model assessment (default: 0.7) |
| `--partial-hit-identity` | Minimum homolog identity percent for accepted partial models (default: 35.0) |
| `--complete-hit-bits` | Minimum whole-homolog bitscore for completeness support (default: 120.0) |
| `--complete-query-coverage` | Minimum query coverage for whole-homolog completeness support (default: 0.8) |
| `--complete-reference-coverage` | Minimum subject coverage for whole-homolog completeness support (default: 0.85) |
| `--complete-terminal-slack-aa` | Minimum allowed terminal slack in whole-homolog alignment, aa (default: 25.0) |
| `--complete-terminal-slack-fraction` | Allowed terminal slack as fraction of reference length (default: 0.1) |
| `--complete-evidence` | Completeness evidence: qualified-template (rc3 default), external-homology (whole active external OR homolog), or coding-and-identity (intact coding plus OR identity; limited full-length confidence) (default: qualified-template) Choices: qualified-template, external-homology, coding-and-identity. |
| `--reconstruction-reference-policy` | Local reconstruction references: audited (default) or all active external OR proteins; self-generated search seeds never qualify (default: audited) Choices: audited, active-external. |

## Repair admission thresholds

| Option | Meaning |
| --- | --- |
| `--repair-hit-bits` | Minimum external-reference bitscore to plan a polish repair (default: 80.0) |
| `--repair-hit-evalue` | Maximum external-reference E-value to plan a polish repair (default: 1e-10) |
| `--repair-hit-identity` | Minimum external-reference identity percent to plan a polish repair (default: 30.0) |

## Preset differences

All unlisted settings use G0 values in every preset. An explicitly supplied
option overrides the selected preset. Values below refer to one full workflow run.

| Parameter | G0 | G1 | G2 | G3 |
| --- | --- | --- | --- | --- |
| `--first-seed-hits` | 6 | 5 | 5 | 5 |
| `--second-seed-hits` | 3 | 2 | 2 | 2 |
| `--first-output-coverage` | 0.25 | 0.2 | 0.2 | 0.2 |
| `--second-output-coverage` | 0.15 | 0.1 | 0.1 | 0.1 |
| `--first-alignment-ratio` | 0.3 | 0.25 | 0.25 | 0.25 |
| `--second-alignment-ratio` | 0.2 | 0.15 | 0.15 | 0.15 |
| `--alignment-output-score` | 0.2 | 0.15 | 0.15 | 0.15 |
| `--candidate-min-aa` | 80 | 60 | 60 | 60 |
| `--candidate-bits` | 35.0 | 30.0 | 30.0 | 30.0 |
| `--candidate-query-coverage` | 0.35 | 0.25 | 0.25 | 0.25 |
| `--adaptive-require-full-reference` | enabled | enabled | disabled | disabled |
| `--adaptive-reference-coverage` | 0.7 | 0.7 | 0.65 | 0.65 |
| `--evidence-or-targets` | 50 | 50 | 100 | 100 |
| `--evidence-gr-targets` | 5 | 5 | 20 | 20 |
| `--identity-bits` | 50.0 | 50.0 | 45.0 | 45.0 |
| `--identity-min-identity` | 25.0 | 25.0 | 20.0 | 20.0 |
| `--identity-query-coverage` | 0.5 | 0.5 | 0.4 | 0.4 |
| `--distant-bits` | 35.0 | 35.0 | 30.0 | 30.0 |
| `--distant-query-coverage` | 0.35 | 0.35 | 0.3 | 0.3 |
| `--identity-use-tree` | disabled | disabled | enabled | enabled |
| `--reconstruction-reference-policy` | audited | audited | active-external | active-external |
| `--structure-max-loci` | 500 | 500 | 1000 | 1000 |
| `--structure-references` | 8 | 8 | 10 | 10 |
| `--structure-rescue-budget` | 64 | 64 | 128 | 128 |
| `--recall-max-regions` | 500 | 500 | 1000 | 1000 |
| `--recall-complementary-regions` | 100 | 100 | 200 | 200 |
| `--translated-max-regions` | 5000 | 5000 | 7500 | 7500 |
| `--polish-repair-budget` | 256 | 256 | 512 | 512 |
| `--repair-hit-bits` | 80.0 | 80.0 | 60.0 | 60.0 |
| `--repair-hit-evalue` | 1e-10 | 1e-10 | 1e-08 | 1e-08 |
| `--repair-hit-identity` | 30.0 | 30.0 | 25.0 | 25.0 |
| `--model-hit-evalue` | 1e-20 | 1e-20 | 1e-12 | 1e-12 |
| `--model-hit-bits` | 100.0 | 100.0 | 80.0 | 80.0 |
| `--model-hit-identity` | 30.0 | 30.0 | 25.0 | 25.0 |
| `--model-hit-coverage` | 0.7 | 0.7 | 0.65 | 0.65 |
| `--partial-hit-identity` | 35.0 | 35.0 | 30.0 | 30.0 |
| `--complete-evidence` | qualified-template | qualified-template | qualified-template | external-homology |
| `--complete-hit-bits` | 120.0 | 120.0 | 120.0 | 100.0 |
| `--complete-query-coverage` | 0.8 | 0.8 | 0.8 | 0.75 |
| `--complete-reference-coverage` | 0.85 | 0.85 | 0.85 | 0.8 |
| `--complete-terminal-slack-aa` | 25.0 | 25.0 | 25.0 | 35.0 |
| `--complete-terminal-slack-fraction` | 0.1 | 0.1 | 0.1 | 0.15 |

## Other commands

`orsearcher check` verifies analysis tools. `orsearcher --version` prints the
software version. The remaining commands accept the options below. `refine`
and `polish` also use the relevant input and tuning options listed above.

### `orsearcher references`

| Option | Meaning |
| --- | --- |

### `orsearcher prepare-benchmark`

| Option | Meaning |
| --- | --- |
| `--family-map` | Comparator-only TSV: feature_id, family (OR/GR/OTHER), optional name; explicit labels override automatic annotation names |
| `--or-name-pattern` | Additional case-sensitive regex for OR comparator identifiers; repeatable |
| `--gr-name-pattern` | Additional case-sensitive regex for GR comparator identifiers; repeatable |

### `orsearcher benchmark`

| Option | Meaning |
| --- | --- |
| `--sample` | Select one sample when --prediction is a batch output directory |
| `--prediction-scope` | Auto prefers final results; explicitly select legacy or unverified discovery candidates Choices: auto, final, legacy, candidates. |
| `--seqid-map` | Explicit TSV: prediction_seqid, annotation_seqid; never strip accession versions implicitly |
| `--min-reference-coverage` | Primary locus-recovery threshold; coverage tiers 50/75/90% are also reported |
| `--min-predicted-support` | Minimum fraction of predicted coding bases supported by one annotated transcript |
| `--family-map` | Comparator-only TSV: feature_id, family (OR/GR/OTHER), optional name; explicit labels override automatic annotation names |
| `--or-name-pattern` | Additional case-sensitive regex for OR comparator identifiers; repeatable |
| `--gr-name-pattern` | Additional case-sensitive regex for GR comparator identifiers; repeatable |

### `orsearcher refine`

| Option | Meaning |
| --- | --- |
| `--candidate-models` | Annotation-blind model list JSON; repeat for stages; never a target GFF |
| `--evidence-cache` | Read-only compatible protein SQLite cache copied to the new output |

### `orsearcher polish`

| Option | Meaning |
| --- | --- |
| `--previous-sample` | Previous samples/<sample> directory containing structure/all_models.json |
| `--allow-unverified-previous` | Explicitly accept missing/changed previous-genome metadata (not proof of equivalence; prefer refine) |

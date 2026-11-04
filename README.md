# ORsearcher: A Specialized Tool for Odorant Receptor Prediction

```text
   ___  ____                      _
  / _ \|  _ \ ___  ___  __ _ _ __| |__   ___ _ __
 | | | | |_) / __|/ _ \/ _` | '__| '_ \ / _ \ '__|
 | |_| |  _ <\__ \  __/ (_| | |  | | | |  __/ |
  \___/|_| \_\___/\___|\__,_|_|  |_| |_|\___|_|
```

**Version 2.3.0** · Insect odorant receptor annotation from a genome assembly

ORsearcher identifies candidate odorant receptor (OR) loci, reconstructs their protein-coding gene models, and evaluates supporting evidence from OR, gustatory receptor (GR), and other protein sequences. The tool outputs complete gene models, partial gene models, and candidate pseudogenes, together with the evidence supporting each classification. ORsearcher includes a built-in reference database containing 3,309 OR protein sequences from 46 species, along with corresponding GR and background protein databases, allowing users to initiate searches using only a target genome. When the built-in OR reference sequences are evolutionarily distant from the target species, users can supplement the reference database with additional sequences through configurable parameters to improve search performance.

[Installation](#installation) · [Quick start](#quick-start) · [Outputs](#outputs) ·
[Reference file formats](#reference-file-formats) · [Parameters](#parameters)

## Installation

Use Linux or WSL2 with Python 3.10 or later. From the source directory, create
the environment and install the package.

```bash
mamba env create -f environment.yml
conda activate orsearcher23
python -m pip install .
orsearcher --version
orsearcher check
```

The package includes 3,309 curated OR proteins from 46 species across 14
orders, 1,530 UniProtKB GR proteins, 7,739 background proteins and three Pfam
profiles. `orsearcher references` shows the installed files and their
provenance. All default prediction inputs are available offline after installation.

## Quick start

### 1. Annotate one genome

```bash
orsearcher --genome insect.fna
```

This runs G0 with a 36-thread CPU budget. For `insect.fna`, the sample is
`insect` and the default output directory is `orsearcher_insect/`. Use a
fresh `--outdir` for another analysis. Use `--threads` on smaller machines.

```bash
orsearcher --genome insect.fna --threads 12 --outdir insect_g0
```

The `orsearcher run --genome insect.fna` form is equivalent. Use plain
nucleotide FASTA or BGZF-compressed FASTA. `samtools faidx` may create an index
beside the input genome, so use a writable genome directory or an existing index.

The default uses the full curated OR panel. For a held-out test, exclude
the target species from the OR panel and its template subset with an explicit
reference-header pattern, for example `--exclude-reference-pattern 'species=Drosophila_melanogaster\b'`.
For study reproduction, also select the study's matching GR and background inputs.

### 2. Choose a parameter preset

| Preset | Settings |
| --- | --- |
| `G0` | Baseline search, reconstruction, identity and completeness criteria. |
| `G1` | Relaxes discovery and candidate admission. |
| `G2` | Adds more permissive homology and reconstruction settings, larger local budgets and eligible tree evidence. |
| `G3` | Adds broader completeness support from active external OR homologs and relaxes the whole-alignment cutoffs. |

```bash
orsearcher --genome insect.fna --preset G1 --outdir insect_g1
orsearcher --genome insect.fna --preset G2 --outdir insect_g2
orsearcher --genome insect.fna --preset G3 --outdir insect_g3
```

All four presets run the same workflow once. Explicit options override a
preset regardless of their position on the command line. The four settings
were compared in fruit fly with the study reference panel. Preset choice and
reference coverage can affect which loci are recovered and how models are classified.

### 3. Replace or add reference inputs

| Input | Replace the family input | Add to the active family input |
| --- | --- | --- |
| OR proteins | `--replace-or FILE` | `--add-or FILE` |
| GR proteins | `--replace-gr FILE` | `--add-gr FILE` |
| Profiles | `--replace-hmm OR FILE` or `--replace-hmm GR FILE` | `--add-hmm OR FILE` or `--add-hmm GR FILE` |

```bash
orsearcher --genome insect.fna --replace-or my_ORs --outdir insect_custom_or
orsearcher --genome insect.fna --preset G2 --add-or close_ORs --add-gr close_GRs --outdir insect_augmented
orsearcher --genome insect.fna --replace-hmm GR my_GR_profiles --add-hmm OR extra_OR_profile --outdir insect_profiles
```

Replacement selects the baseline for that role. Additions extend the baseline
selected for the run. Repeat `--add-or`, `--add-gr` and either HMM option to
provide multiple files. Changes apply to the current run and are recorded in
`input_manifest.json` and `shared/reference_inputs/selection.json`.

### 4. Adjust a parameter or run another round

```bash
orsearcher --genome insect.fna --preset G2 --structure-max-loci 750 --outdir insect_custom
orsearcher --genome insect.fna --preset G0 --rounds 2 --outdir insect_rounds \
  --next-round --preset G3 --add-or extra_ORs
```

Each round uses the same workflow. Explicit options carry over unless
overridden in a `--next-round` block. A new preset sets the baseline for
unspecified options in that round. Intact models from the preceding round can
seed discovery. Use `--round-seeds off` to use external references only.
See the [CLI reference](docs/CLI.md) for multi-sample and cached workflows.

## Outputs

| Path under the output directory | Contents |
| --- | --- |
| `samples/<sample>/results/OR_catalogue.tsv` | Main catalogue with coordinates, class, evidence and completeness blockers. |
| `samples/<sample>/results/OR_models.gff3` | Coordinates and CDS structure for exported models. |
| `samples/<sample>/results/OR_{complete,partial,pseudogene}.*` | Separate GFF3, CDS FASTA and protein FASTA files for each class. |
| `samples/<sample>/results/summary.json` | Counts of independent loci and exported models, classification and review status. |
| `samples/<sample>/results/export_validation.json` | Checks of exported model coordinates and sequences. |
| `samples/<sample>/recall_candidates/` | Discovery candidates, including models that did not qualify for final export. |
| `samples/<sample>/phylogeny/` | OR/GR tree and alignment when generated. |
| `batch.log`, `samples/<sample>/run.log`, `debug.log` | English run progress, sample progress and detailed diagnostics. |
| `effective_parameters.json`, `input_manifest.json` | Resolved preset, settings, inputs and reference checksums. |

`batch_summary.json` records run status. Multiple complete rounds write the
same layout under `rounds/round01/`, `rounds/round02/` and so on. See
[output details](docs/OUTPUTS.md) for intermediate files and interpretation.

## Reference file formats

**OR and GR proteins.** Supply standard protein FASTA. The header needs only
an identifier after `>`. Descriptions, accession and species information are
optional, and their order is unrestricted. ORsearcher records the original
headers and gives supplied records internal identifiers. It accepts wrapped
or unwrapped sequences, lower case and a single terminal `*`. Use amino-acid
sequences of at least 50 residues without gaps or internal stops. See the
[FASTA example](examples/reference_proteins.fasta). The example illustrates
format only.

Submit both full-length proteins and fragments in the same OR file. ORsearcher
checks protein length and profile support to decide which supplied sequences
can support completeness. A fragment can still help find a candidate locus.

**HMMs.** Supply amino-acid profiles in HMMER3 text format. A file can contain
multiple profiles ending in `//`. Specify the role (`OR` or `GR`) before each
filename. Use unique `NAME` values for distinct profiles in one family.
ORsearcher applies each profile's own gathering cutoffs where present and uses
E-value and coverage evidence for profiles without gathering cutoffs. File
extensions are unrestricted for both HMMs and protein FASTA.

## Parameters

Run `orsearcher -h` for the option list. The [CLI reference](docs/CLI.md)
describes every option by module and gives the exact G0–G3 settings.

| Module | Example options |
| --- | --- |
| Input and scheduling | `--genome`, `--threads`, `--preset`, `--jobs`, `--workers` |
| Genome search | `--chunk-size`, `--overlap`, `--first-seed-hits`, `--candidate-min-aa` |
| TBLASTN and MetaEuk searches | `--translated-max-regions`, `--recall-max-regions`, `--recall-sensitivity` |
| Structure reconstruction | `--structure-max-loci`, `--structure-references`, `--reconstruction-reference-policy` |
| Identity and completeness | `--identity-bits`, `--family-score-margin`, `--identity-use-tree`, `--complete-evidence` |
| Repair and export | `--polish-repair-budget`, `--polish-timeout`, `--resume` |

See [reference provenance](docs/DATA_PROVENANCE.md), [methods](docs/METHODS.md)
and [release verification](docs/RELEASE_VALIDATION.md). Code is available
under the [MIT License](LICENSE). The bundled data have separate terms.

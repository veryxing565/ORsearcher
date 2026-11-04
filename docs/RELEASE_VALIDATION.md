# ORsearcher 2.3.0 release verification

The release was verified on x86-64 Linux under WSL2 using the supplied
`environment.yml`. All 303 regression tests passed with the native tools
installed. The 72 parameters in the prior G0–G3 comparison table matched
the resolved settings for each preset. Explicit option precedence was also
checked.

## Installed analysis

The final bundled OR FASTA and its template subset match the frozen full
study panel byte for byte. The final wheel was installed in a fresh Python
environment, all 303 regression tests passed, and a G0 run with the genome
as its sole option recovered both planted OR loci as complete models with
zero coordinate-to-sequence mismatches. The bundled GR, background and
profile files match those tested in the previous wheel. The final run
retained audit warnings for unresolved GeneWise anchors, anchors without
qualified GeneWise references, and competing overlapping paths.

### Earlier interface verification

The wheel was installed and run from outside the source directory. A small
synthetic genome contained two planted spliced OR models on opposite strands
and a negative sequence. The G0 command supplied only the genome:

```bash
python tools/make_smoke_genome.py --out smoke.fna
orsearcher --genome smoke.fna
```

The interface checks below used the earlier UniProtKB OR panel with the
same executable and the same GR, background and profile files.

| Preset | Reference settings exercised | Independent loci | Exported complete models |
| --- | --- | ---: | ---: |
| G0 | All bundled inputs, 36-thread default | 2 | 2 |
| G1 | OR and GR additions and an added OR profile | 2 | 2 |
| G2 | Replacement OR and GR proteins and replacement GR profile | 2 | 2 |
| G3 | All six replacement and addition options plus an explicit threshold override | 2 | 2 |

All four analyses ran the two Miniprot passes, TBLASTN, MetaEuk/GeneWise,
structural review, OR/GR phylogeny and final export. Each matched the assembly
for the exported coordinates and sequences. G0 and G1 finished with recorded
warnings for unresolved GeneWise anchors. Their final outputs and audit
records were available. G2 and G3 completed without those warnings.

Custom input checks exercised extensionless FASTA and HMM files, a gzip FASTA
with no filename suffix, duplicate FASTA IDs, and a profile lacking gathering
cutoffs. A supplied OR fragment was retained for search and excluded from
template eligibility. The regression tests also check malformed HMM and DNA
input errors.

An additional installed run completed two rounds with a G0-to-G3 preset
change and inherited replacement inputs. The G3 run with both replacements
and additions resumed with the same input files and settings.

The wheel, source distribution and GitHub source ZIP were checked for
all seven reference files, their recorded checksums and public data provenance.

## Scope

The planted proteins were represented in the reference inputs. These tests
check installation, option behavior, pipeline execution and exports. The
study's genome-scale benchmarks used their separately recorded inputs and
exclusions. The bundled UniProtKB GR panel has sparse orthopteran coverage.
Exact inputs and data terms are documented in the
[bundle manifest](../src/orsearcher/data/defaults/manifest.json).

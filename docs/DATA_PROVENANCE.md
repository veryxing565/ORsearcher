# Bundled reference data

The public 2.3.0 package includes a self-contained reference bundle. The
[manifest](../src/orsearcher/data/defaults/manifest.json) records the exact
files, counts, source queries, accessions, selection rules and SHA-256 hashes.
`orsearcher references` reports the bundle selected at runtime.

## Protein panels

The OR panel is the curated OR reference database v1.0.0 compiled from
published insect OR datasets, iOrME and deposited protein translations. It
contains 3,309 search proteins from 46 species in 14 orders, including 915
partial-rescue entries and 38 Orco entries. The bundled
[sequence source table](../src/orsearcher/data/defaults/OR_sources.tsv) records
the original source, available accessions, species, curation pool and template
eligibility of each entry. The maintainers include this panel with permission
to redistribute it in ORsearcher. Individual source datasets retain their terms.

`complete_OR.faa` contains the original 1,373 templates. These entries have
an author complete or curator nonfragment/nonpseudogene label, are 250–650
amino acids long and passed the collection's residue and curation checks.
Template eligibility records the level of support available for reconstruction.

The GR panel contains 1,530 proteins selected from UniProtKB release 2026_03.
Candidate entries were selected by receptor name in eight insect orders and
screened against the GR Pfam profiles. The selection removes identical
protein sequences, limits contribution by taxon and order, and retains the
UniProt accession and taxonomy ID in each header. Source queries and counts
are in the manifest.

`background.faa` contains previously filtered reviewed insect UniProtKB
proteins used for other-family competition. It is not an exhaustive list of
non-OR/non-GR proteins. The GR and background files are modified derivatives of UniProtKB
data, distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Credit the [UniProt Consortium](https://www.uniprot.org/) and refer to the
[UniProt data terms](https://www.uniprot.org/api-documentation/support-data).

## Profile models

The bundle uses Pfam PF02949.26 for OR evidence and PF08395.18 and PF06151.19
for GR evidence. The files are the pinned downloads from InterPro. Pfam data
are distributed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
See the [InterPro data notice](https://proteinswebteam.github.io/interpro-blog/2023/03/15/InterPro-93.0-new-features-and-updates/).

## Reproducing earlier analyses

The study's held-out tests used sample-specific subsets of the bundled OR
panel and its templates. The study used a different GR library and
sample-specific exclusions on GR and background proteins. The experimental
GR derivative was derived from the supplementary
data of [Mitchell et al. (2026)](https://doi.org/10.1038/s41467-026-73718-5) and
is kept separate from the public release. Use the original approved inputs,
checksums and command settings to reproduce an earlier run. The
[reference guide](REFERENCES.md) explains how to select an external bundle.

The MIT software license covers project code. Data retain the terms listed
above and in the manifest.

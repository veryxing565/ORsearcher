# Contributing

Please keep algorithmic changes separate
from dataset changes, and include regression tests for both positive and
negative cases.

## Development

Use Python 3.10 or newer on Linux. Install the native analysis tools with the
provided `environment.yml` when testing real analyses. The pure-Python suite
does not require those tools; native R tests are skipped when R is absent.

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python -m pip install build
python -m build
```

## Reporting a problem

Include the software version, operating system, command line, dependency
report, relevant concise log excerpt, and the failure status file. A small
synthetic reproducer is preferred. Do not attach passwords, SSH keys, cloud
tokens, private reference collections, or unpublished genomes without consent.

Distinguish execution failures from scientific disagreement: a local alignment
timeout is not evidence of a pseudogene; several models at one locus are not
necessarily several genes. Preserve raw model provenance and held alternatives
when changing consolidation or classification rules.

## Reference data and releases

New reference datasets need source URLs/accessions, collection date, sequence
identifiers, completeness policy, checksums, and license/redistribution status.
Do not silently label a public supplement or database download as MIT. A
release must pass tests, build/install outside the checkout, include all
required runtime data, and contain no machine-specific paths or credentials.
Full-genome benchmarks must identify their target genomes and distinguish
same-species-reference performance from held-out evaluation.

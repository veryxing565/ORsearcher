# Reference selection

The installed defaults provide all reference roles for a new prediction.
Use `orsearcher references` to inspect the active bundle. The
[README](../README.md#reference-file-formats) describes input formats.

| Role | Bundled file | Public option |
| --- | --- | --- |
| OR proteins | `OR.faa` | `--replace-or`, `--add-or` |
| GR proteins | `GR.faa` | `--replace-gr`, `--add-gr` |
| OR profiles | `OR.hmm` | `--replace-hmm OR`, `--add-hmm OR` |
| GR profiles | `GR_1.hmm`, `GR_2.hmm` | `--replace-hmm GR`, `--add-hmm GR` |
| OR template subset | `complete_OR.faa` | Inferred for replacement OR proteins. Audited additions can extend the subset. |
| Other protein competition | `background.faa` | `--background-reference` |

Use `--reference-dir DIR` or `ORSEARCHER_REFERENCE_DIR` for a compatible
external bundle with the same filenames. The command-line directory takes
priority over the environment setting. A replacement overrides the selected
bundle for that role. If both replacement and addition are supplied, the
addition extends the replacement. HMM replacement is scoped to the family
specified on the command line. Multiple model files can be used for each
family.

For legacy reference identifiers, `--or-reference`, `--gr-reference`,
`--or-hmm` and `--gr-hmm` continue to accept explicit input paths.
`--complete-reference` pairs with `--or-reference` when an independently
prepared subset is available. Subset entries must match the active OR
reference by ID and sequence. Use `--complete-reference none` to select
templates through the runtime audit.

See [data provenance and terms](DATA_PROVENANCE.md) for the public files.

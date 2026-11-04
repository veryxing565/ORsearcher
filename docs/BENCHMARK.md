# Annotation benchmarking (scorer schema 2)

The scorer compares coordinates with a held-out GFF3 annotation. It does not
run genome prediction, supply target annotation to the predictor, establish
functional OR identity, or prove that every unannotated candidate is false.
Never compare old/new summary numbers without noting the scorer schema and
rescoring both with the same inputs and rules.

## Select the final outputs

```bash
orsearcher benchmark --annotation curated_same_assembly.gff3 \
  --prediction run_A/samples/insect_A --outdir score_A
```

The default priority is `results/` (polished final models), directly supplied
`OR_catalogue.tsv` + `OR_models.gff3`, legacy `nonredundant/`, then a directly
supplied candidate pool. A half-present preferred pair is an error, not a
reason to fall back silently. A directly supplied `results/` folder also works.
Use `--prediction-scope final` to require final-model filenames. Candidate
scopes are always labelled unverified, with zero complete-model recoveries.

For a batch root use `--sample insect_A`; omission works only if exactly one
sample directory exists. Explicit `--prediction-scope legacy` or `candidates`
allows analysis of intermediate evidence without silently replacing the final
result. The output directory must be separate and empty; old scores are not
overwritten. Examples below use a new output directory on every call.

## Cross-insect reference selection

Automatic symbols include `Or67d`, `OR1`, `Orco`, `LmigOR1`, `CkiaOr123` and GR
equivalents. Species prefixes without separators require a capitalized family
boundary. Gene/transcript/CDS product labels such as `odorant receptor 1` also
work when the gene is named `LOC...`. Binding/associated/accessory proteins
are not selected from their product phrase. This is deliberately conservative:
it cannot identify every arbitrary gene ID or guarantee that the source
annotation correctly calls every OR/GR.

For a curated list, provide `--family-map families.tsv` to `benchmark` or
`prepare-benchmark`. This maps annotation IDs/aliases, **not protein similarity**:

```text
feature_id	family	name
gene-LOC123	OR	SpeciesOR1
XP_012345.1	GR	SpeciesGR1
gene-LOC999	OTHER	not_a_receptor
```

Keys may match exact feature `ID`, `gene`, `Name`, `gene_name`, `locus_tag`,
`protein_id`, `transcript_id`, `Alias` or full `Dbxref` tokens. A matched CDS or
transcript assigns its parent gene. Explicit labels override automatic names
for that gene; conflicting explicit families, absent keys or keys matching
multiple genes are errors. `OTHER` excludes a misleading automatic label.
Orco is included as OR by default; explicitly map it to OTHER if the declared
benchmark excludes Orco. Additional naming schemes may use repeatable
`--or-name-pattern '^custom_OR'` / `--gr-name-pattern '^custom_GR'` regexes;
audit the selected inventory before running comparisons.

`prepare-benchmark --annotation annotation.gff3 --family-map families.tsv
--outdir comparator` extracts the selected comparator and inventory. Keep it
outside prediction inputs. GFF3 or `.gff3.gz` is supported, not GTF. Parent
chains, shared CDS parents, CDS directly under genes and transcript isoforms
are handled; a transcript must resolve to one gene. Unsupported or ambiguous
hierarchies fail rather than silently choosing a parent.

Declared multipart/trans-spliced non-receptor genes do not block whole-file
reading (for example Drosophila `mod(mdg4)`). A selected receptor split across
different contigs/strands remains an explicit unsupported-comparator error;
the scorer does not invent an ordinary cis-spliced structure for it.

## Coordinate identifiers and input safety

Annotation and predictions must describe the **same assembly and coordinate
system**. Accessions including their versions must agree. An optional
`--seqid-map contigs.tsv` performs explicit name substitution only:

```text
prediction_seqid	annotation_seqid
chr1	NC_012345.1
```

It is not a liftover and does not prove the sequences are the same. The program
does not guess by stripping accession versions or `chr` prefixes. Completely
disjoint annotation/prediction contig IDs fail; predictions on valid contigs
without any annotated receptor remain unmatched, not invalid.

Duplicate candidate IDs, missing CDS, orphan CDS parents, inconsistent
contigs/strands/bounds and invalid coordinates fail. Repeated identical CDS
rows are deduplicated; interval unions prevent double-counting coding bases.
Overlapping non-identical CDS remain flagged and cannot count as exact intact
structure. Such structures may be legitimate slippage annotations or errors;
they are not silently repaired by the scorer.

## What the numbers mean

- Primary locus recovery: at least 50% of one annotated OR transcript's CDS,
  with at least 30% of predicted coding bases supported by that transcript.
  `--min-reference-coverage` / `--min-predicted-support` can change the primary
  rule; record those values and use the same rule for all tools.
- Additional 50%, 75% and 90% recovery tiers use the declared predicted-support
  cutoff. Coverage uses CDS interval unions, not genomic span or protein length.
- A deterministic maximum-cardinality one-to-one matching allows each model
  to claim at most one gene and each gene at most one model. Isoforms do not
  inflate gene counts. Tiers are solved independently, so assignments can vary
  in ambiguous cases. Potential multi-gene models and duplicate predictions
  are flagged; their complete overlap evidence is retained.
- `independent_locus=False` alternatives are retained in diagnostics but not
  primary recovery. Older files without that field are treated as exported
  representatives; one-to-one matching still prevents duplicate locus counts.
- Strict structure recovery requires a `complete` prediction with identical
  CDS boundaries and phases to a comparator transcript, known phases, and no
  comparator partial/pseudogene flags. Nonzero leading CDS phase also prevents
  treating the comparator as intact. The separate within-3-bp endpoint flag
  is only diagnostic. **Neither coordinate agreement nor the complete label
  proves biological full length or function.**
- OR genes without an annotated CDS are listed as unscorable, not silently
  counted as missed coding loci. A zero coding-OR denominator is an error.
- All qualifying GR overlaps are reported, even when an OR match scores
  better. An overlap is an investigation flag, not proof of GR identity. If no
  coding GR comparator is recognized, the summary explicitly warns that zero
  overlaps cannot demonstrate GR specificity.
- Unmatched predictions do not automatically supply a false-positive rate,
  precision or FDR. Those require an independently justified truth/negative set.

Outputs: `benchmark_summary.json` (schema, thresholds, input SHA-256 and counts),
`OR_locus_recovery.tsv` (one row per scoreable OR locus),
`candidate_annotation_comparison.tsv` (best diagnostic match plus assigned locus),
`all_annotation_overlaps.tsv` (all transcript-level overlaps),
`missed_OR_annotation.tsv` and `unscorable_reference_loci.tsv`.

## 中文：现在能否冻结并开始论文比较？

结论：评分修正后可以进入**少量开发物种校准**，不宜直接宣布软件已达到昆虫纲通用准确率，也不宜立刻冻结开展论文的正式主比较。评分器修好是“尺子能用”，不是“准确率已经测好”。

建议先用 3 个物种、至少 3 个目，例如果蝇、赤拟谷盗、家蚕或蜜蜂；实际以已有基因组、同版本注释和可追溯 OR/GR 目录为准。果蝇已经反复参与开发，只能作为开发检查，不再冒充独立测试物种。先核查最终位点找回、75%/90% 覆盖、严格结构一致、GR 重叠及耗时/内存/超时。只处理暴露的明显错误，不为个别物种无限调参。

冻结前仍需明确一件重要的事：程序当前默认 OR 库是 **390 条、4 种直翅目**；2026-10-28 新整理的昆虫库 v1 是 **3309 条、46 物种、14 目**，尚未替换这个默认。新库中的 primary 记录不都等于全长模板，partial/Orco 也不能不经检查就用作全长修补依据。应固定“搜索全集”和“合格全长模板子集”各自的文件及校验值；不得把新数据库整理完成误写成新库已经经过完整流程验证。

通过开发检查后，再冻结代码、OR/GR/HMM/背景库、完整模板、参数、依赖版本和输入校验值。正式测试另留未参与调参的物种，不根据正式结果反复调参再继续称它们为独立测试。常规跨物种测试去掉目标物种参考及换名后的完全相同蛋白；另做留出整个目的远缘压力测试，分开报告。各工具采用相同原始基因组、可比参考范围与明确资源/时间预算。目标注释只用于评分，不进入预测。

研究使用与公开分发是两回事：新 v1 库逐来源再分发许可尚未审计完成，现有 GR 衍生库不随公开包分发。可以先开展研究校准，但公开 GitHub/投稿复现材料需保留署名、明确许可，必要时用下载与本地准备脚本替代直接打包。

For GFF3 attribute semantics see [NCBI's official GFF3 guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

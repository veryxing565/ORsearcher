# 2.3.0

- Add a single-genome entry point with a 36-thread default and G0–G3 presets.
- Bundle the 3,309-protein curated OR panel, a UniProtKB GR panel, Pfam profiles and background proteins.
- Support run-local reference replacement and augmentation with standard protein FASTA or HMMER3 profiles.
- Record reference input provenance, resolved settings and English progress logs.

# 2.3.0rc4 (2026-11-01)

CLI-first tuning, explicit completeness policies, search-only seeds and repeated normal-pipeline rounds with inherited command-line overrides. Default final assessment replay preserves the prior complete/partial categories. rc3 remains the fixed baseline for the historical competitor benchmark.

# 2.3.0rc2 — 2026-10-29

- Fix annotation benchmarking to prioritize polished `results/` outputs; accept direct result/sample/batch paths, require explicit sample selection for ambiguous batches, and expose explicit final/legacy/candidate scopes.
- Recognize conservative cross-insect OR/GR symbols and product annotations. Add explicit comparator family maps, identifier patterns, gzip GFF3 and user-specified contig aliases; never silently strip accession versions.
- Separate coding-evaluable loci from references without CDS. Validate file pairs, IDs, Parent chains, coordinates/strands and catalogue/GFF bounds; do not silently score malformed input as biological misses.
- Report 50/75/90% coding coverage and strict coordinate/phase agreement separately. Use one-to-one maximum-cardinality matching, retain all overlaps, flag potential fusions/duplicates and exclude non-independent alternatives from primary recall.
- Keep discovery-only candidates distinct from complete models, partial comparator models distinct from intact truth, and GR overlaps distinct from proven false positives. Record scorer schema, parameters and input SHA-256 hashes.
- No search/classification thresholds, predictor algorithms or reference payloads changed. This fixes the measuring tool; it is not evidence of improved biological accuracy.

# 2.3.0rc1 — 2026-10-28

- Integrate genome discovery, structural/identity review, non-transitive coding-overlap consolidation and bounded local repair into one portable package.
- Add standalone `polish` for old reviewed samples, without overwriting old evidence or re-running whole-genome discovery.
- Export complete, partial and assembly-supported pseudogene candidates as protein/CDS/GFF3, one summary catalogue and retained lossless model audits.
- Recompute OR/GR/background competition when replacing reference datasets; explicit options replace defaults instead of appending.
- Bound repair windows and time budgets; cache only successful repairs. Do not classify a timeout as a pseudogene.
- Distinguish overlapping disrupted alternatives from independent loci, including pseudogene-versus-pseudogene overlaps.
- Add portable defaults, provenance/checksums, concise operational logs, debug logs, mamba environment, package builds and CI.
- Preserve the existing locust-focused OR and insect GR reference panel for private use. Exclude the restricted derived GR payload from public packages and provide an explicitly acknowledged local preparation helper.
- This is a release candidate: the software integration tests are not evidence of calibrated accuracy or exhaustive performance across insect orders.

# 2.2.0.dev4 — 2026-10-13

- 独立 candidate10。只修改结构代表选择和审计，不改发现参数及身份阈值；保留所有旧源码和结果。
- 结构比较不再经由 distant_hits 的查询覆盖门槛或 50% 查询覆盖过滤，同一个显著局部 HSP 不因长查询分母变大而消失。按共同参考覆盖、bits/identity 比较，不混排不同参考的原生 identity。
- 共同证据相同后参考编码起止与有效原始锚点来源，不按蛋白长度最大化；融合、非标准剪接、损伤、重复域和完整结构判据保留。
- 新增逐个片段替代的选择审计，查询比例只作诊断并明确不进入排序，不因代表保留自动提升身份或完整类别。
- 新增相同已完成结构组上的选择/新背景/新树对照；源损伤与坐标证据复用、输入哈希和范围明确。它不是新的局部比对或全基因组运行，另做完整后处理的新局部重建对照。
- 新增 15 项回归测试，共 191 项，覆盖相同 HSP/查询长度、位置与参考别名、原生指标混排、较短但更强证据的路径及损伤/融合约束。

# 2.2.0.dev3 — 2026-10-13

- 独立 candidate09；保留 candidate07/08 的源码、运行目录和发现证据。
- Profile 阴性 OR 新增“强、近全长 curated 同源 + 不同蛋白序列的显著辅助同源 + OR/GR 分数优势 + 无 GR profile”的路径；单条强命中仍不能确认，真实 GR/背景冲突仍待定。减少这类模型对一次诊断树正票的依赖，并非宣称所有树判定均已稳定。
- 明确区分结构复核准入与最终 OR 身份；弱显著 curated 候选可进入修补，但不能据此进入最终 OR 集。结构选择和审计不再跳过待定身份组。
- 片段代表按共同参考蛋白的同源覆盖和分数比较，不混排不同参考的原生 identity/ref_cov，也不按最长蛋白选择。非标准剪接、损伤、融合和重复域不凭长度升类。
- 默认最多 64 组、每组 20 秒、最多 3 条不同合格新参考的 Exonerate 补证；使用种子长度 6，不额外再开 legacy 重试。超时/延期/缺少新参考均记录，仍是同一个真实比对器，不增加独立工具票。
- Complementary legacy 种子额度优先用于未解决结构及近缘损伤冲突，不再仅按原始 OR bits 分配；空间预算轮访和明确延期保留。
- 主 Exonerate 段错误重试只用剩余时限。额外重建的移码事件参与后续损伤审计，不修改组装碱基。
- 新增缓存病例诊断工具；对照脚本接受基线版本标签，完整类别精确结构增减与任意类别精确结构增减分别记录。
- 原生联调发现新比对结果在补证前尚无融合标记；先初始化空间融合证据再选择补证，失败验证目录保留。缓存诊断无新树，不交付 reviewed_complete 标签。
- 新增 29 项回归测试，共 176 项。程序效果需要原生联调、真实蛋白负控和两组果蝇对照完成后评估，不提前宣称恢复了全部缺失 OR。

# 2.2.0.dev2 — 2026-10-12

- 独立 candidate08，不改写 candidate07 源码、结果或证据缓存。
- 树分歧但无真实 GR 竞争不再否决多渠道 OR；远缘路径要求强 domain、显著审计同源/参考共识和对 GR 的优势，同等放宽 GR 竞争计分以避免不对称。真实竞争冲突仍保留。
- 合格原生比对进入修补但不确认身份；预算轮访固定空间锚点，端点兼容优先于小幅蛋白分数收益。
- 同框、无损伤、剪接兼容、近包含完整代表片段技术别名；不合并损伤片段、不同框或串联基因。输出明确模型记录不是确认基因数。
- 共享载荷结构模型存储，原坐标/来源/比对器独立，读取兼容旧版；减少重复大 JSON。身份审计在最终身份引起的结构置信度调整后写出。
- 新增 exact-sequence-disjoint 文献 OR/GR/reviewed 其他蛋白身份检验；新增对 candidate07 同口径接受位点/结构增减比较。
- 新增 16 项回归测试，合计 147 项。原树否决单测按新的明确策略更新，真实 GR 冲突用例仍保留；尚不宣称昆虫纲泛化、最终特异性或全长恢复率已提高。

# 2.2.0.dev1 — 2026-10-07

- 独立 candidate07，保持现有发现层，增加默认统一的结构优先/身份后判定；新 `refine` 命令可只读发现缓存并在新目录重做注释。
- 完整外观模型也复核；结构预算明确，弱同源可进入修补但不自动成为 OR；全蛋白兼容支持优先于长度，参考额度按不同蛋白序列计，允许审计合格的长参考参与修补。
- OR/GR profile 细节与 GR 同源列表留存；profile、竞争同源、经 QC 的参考树与 reviewed 非 OR 背景分渠道综合，冲突/缺证据保留身份待定；排除不擦除原始候选。
- 背景仅作蛋白级竞争；官方 UniProt 下载、阳性去污染、OR/GR profile 审计、来源和哈希独立记录，不增加全基因组参考搜索负担。
- 参考树新增最近三个参考一致性、参考分支错误率和占位 QC；树与对齐命令有超时，失败不作自动身份否定。
- 结构窗口上限防止越界；缺失/超时局部输出不使邻居复核崩溃。假基因仍需相同组装的可重复损伤证据，不修改碱基、不把终止/移码提示直接确认为假基因。
- 新增 20 项身份、冲突、参考序列去重、修补优先级、可运行额度、缺失输出和 CLI 回归测试，另加 2 项服务器原生 R 树 QC 测试；合计 131 项（无 R 的本地环境跳过原生两项）。旧参考额度测试更新为不同蛋白序列的额度语义。
- 输入处置映射与所有未复核候选单独保存，避免只输出少量身份待定模型而隐藏低优先级/延期范围。
- 新果蝇两组后处理基准分别检查候选保留、最终接受召回、完整剪接/端点匹配和错误排除；不是完整全基因组重跑。

# 2.1.1.dev4 — 2026-10-07

- 修复长锚点吞并内部短 HSP 的问题：双向跨度重叠去重，保留来源与参考，不传递合并位置。
- 长 MetaEuk 集合按外显子边界生成可运行局部锚点，保留完整原始模型；超过全模型上限的调用仍可提供局部线索，不作为 CDS 导出。
- 窗口/参考预检后再分配预算；不可运行项单列，不占有效重建额度。HSP 与 MetaEuk 通道轮流调度，通道内排名/空间分轮，不将不同原生分数全局混排。
- 分块不能继承父模型任意位置的 HMM 支持；家族证据只辅助优先级，不删除未过滤候选池。
- GeneWise 不充分时可追加有界 Miniprot/Exonerate，默认最多 100 个真正需要的区域，延期和错误明确记录。
- 增加只读全局搜索缓存、全新局部重建的比较脚本，验证证据配置、记录输入来源并保留全部旧候选；不得冒充新的全基因组测试。
- 互补重建严格使用预算检查通过的同一窗口，避免超长 HSP 链在 GeneWise 已回退到局部区域后再次占满窗口或失败。
- 新增 22 项调度、分块、来源、额度与补充重建回归测试，合计 109 项。

# 2.1.1.dev3 — 2026-10-06

- 以召回优先为当前目标：新增 `--recall-only`，候选池不做 OR/GR 删除或完整/假基因定论，所有记录为 candidate/unverified。
- 新增独立 MetaEuk 原始外显子集合搜索，不先贪心移除重叠调用；外部参考 ID 安全重编码并保留映射，坐标恢复到原始 contig。
- GeneWise/BLOSUM45 有界局部重建，不把短模型当成“已经覆盖该基因”，重访弱单参考锚点；保留所有合法结构替代，不修改基因组碱基。
- MetaEuk 集合长度不冒充真实参考覆盖率；GeneWise 多基因输出不将总体查询覆盖分配给每个子模型。
- 独立未过滤候选 TSV/GFF/FASTA、坐标/序列验证、错误和延期状态；与旧最终集和旧全模型候选池分别比较，目标 GFF 只用于参考排除/事后评估。
- 两个集中环境新增包，无其他软件升级。高灵敏度全局参数开销实测偏高，采用可调预算默认，不宣称已恢复全部 OR。
- 新增 17 项候选范围、解析、坐标、去重、来源、超时清理和比较语义回归测试；合计 87 项。固定锚点只合并参考来源，不传递扩展基因座；自生成参考不能遮蔽新搜索。
- 实测兼容 MetaEuk 7 的未指定 CDS phase：仅在每个外显子长度为三的倍数时按零处理；非零 phase/非整密码子拒绝静默转换。Linux 超时清理 MetaEuk 自身进程组，避免遗留 MMseqs 子进程。
- 全模式摘要传播补漏阶段的失败/延期警告；测试分别记录 MetaEuk、GeneWise 补回的基因名称，不将替代模型数当作召回提高。

# 2.1.1.dev2 — 2026-10-06

- 标记 external/self_generated/target_annotation/unknown 比对参考；自生成与目标注释查询不参与独立工具验证、损伤参考计数和高一致性排名。
- 合格外部蛋白的完全匹配用于结构兼容性选择，不冒充独立验证；跨工具参考支持按工具/参考蛋白序列去重，不按终端变体数量投票。
- 去掉模型选择中的长度优先，参考覆盖与外部 OR 证据先于端点推测；近等强内部剪接冲突归入结构未解决。
- 全长蛋白支持只按有/无参与排名，不按命中参考数量投票；修正短模型因数据库参考数量更多而覆盖更好端点模型的问题。
- 外部 Miniprot/Exonerate substantial 结构不一致时，预算内补一次旧长度 6 种子 Exonerate，限制 30 秒并保留两套参数路径；段错误仅在剩余时限内关闭 refinement 重试一次，两个参数方案仍只算一个比对器。
- 独立重建参考按不同蛋白序列选择；互补重建失败/延期进入警告，不覆盖成功的主要路径。
- 新增同框、剪接兼容、近包含的无损伤片段技术别名规则，不拼接或传递合并不同基因。
- 局部结构复核常规使用三条不同参考，首先纳入最强全长同源参考，避免其被高一致性的短局部比对排除。
- 增加真实案例缓存回放和 19 项回归测试，以及 model_selection.tsv/替代模型来源审计；全流程结果完成前不宣称综合性能改善。

# 2.1.1.dev1 — 2026-10-06

- 修复结构复核中 original_terminal 生成之前丢失参考端点和真实比对器来源的问题，避免无同源端点约束的延伸。
- 同参考、同链、查询坐标单调的多 HSP 链参与参考排序与局部窗口规划；查询重置不串联串联拷贝，HSP 链仍不直接输出 CDS。
- 默认局部 Exonerate 蛋白种子长度由 6 改为 3，局部 Miniprot 使用 k=5；有完整编码端点的高覆盖模型才可跳过补充重建。
- 每个补搜锚点最多六条不同的完整参考；所有候选仍通过共享 OR/GR 与 HMM 证据和结构检验，不放宽完整/假基因身份规则。
- 补搜池化保留不同真实比对器、参考和剪接替代模型，结构复核不再只收到单个代表。
- 仅同阅读框、同剪接、近乎相同 CDS 和端点的无损伤片段可记为完整模型技术别名，记录 alias_rule；不按蛋白相似性合并串联基因。
- 区分未解锚点与比对命令失败；窗口过大、无参考等状态写入独立案例记录。
- 果蝇正常参考和严格留出全流程将分别评估。尚不宣称昆虫纲通用性能或远缘全长恢复成功。

# 2.1.0 — 2026-10-02

- Independent directory; 2.0 sources and runs are preserved.
- Audit every OR reference with its profile and provenance. Mask automatic gap-fill discovery queries outside significant OR domain coordinates; original sequences remain unchanged.
- Non-domain matches to automatic references do not establish OR identity. Retain low-priority evidence separately. Sub-GA and curated multi-reference support remain reviewable.
- Full-length reference support requires reconstruction-qualified references. Adaptive seeds also require full-reference support.
- Original/global and local Miniprot are one tool. Reviewed complete requires both Miniprot and Exonerate, identical splice junctions, endpoints within one codon, and non-identical-target reference support.
- Reproducible pseudogene damage requires both actual aligners and distinct reference sequences, not merely two source names or duplicate IDs.
- Plausible longer same-locus coding models prevent unqualified complete classification of shorter models.
- Splice-compatible alternative CDS models share a non-transitive locus representative; originals remain in final/. Residual overlapping representatives lose reviewed-complete confidence.
- Empty local reference sets yield no_reference instead of crashing. Exonerate segmentation faults get one bounded retry without refinement; unsuccessful output is never accepted.
- Maximum second-aligner window default is 500 kb. Oversized/timeout cases remain explicitly unresolved, not pseudogenes.
- R tree diagnostics use iterative postorder counts instead of recursive descendant vectors. Plot failure is a warning. Diagnostic failure does not discard structure results and is explicitly reported.
- Held-out annotation preparation/comparison commands report locus recovery, splice/end consistency, GR overlap and unmatched predictions. Unmatched does not automatically mean false positive.
- This is not a validated universal insect classifier or experimental confirmation of OR/pseudogene function.
# Final test-driven corrections

- Diagnostic trees use the prepared/excluded reference FASTAs, not raw inputs.
- Full-length support is counted by unique reference protein sequence; exact
  target matches are recognized even if DIAMOND's top-hit limit omits an ID.
- Normal-reference and strict target-species-excluded fruit-fly tests are
  reported separately, including sensitivity failures rather than only runs
  with successful exit status.

# 2.1.1.dev0 (2026-10-03)

- 定位果蝇严格留出漏检到发现阶段，增加独立全基因组 TBLASTN 锚点发现及有界局部敏感 Miniprot/Exonerate 重建。
- HSP 不直接作为基因模型、不直接证明损伤；重建模型统一做 OR/GR 竞争、完整性和融合复核。
- 原始坐标池化与边界去重，非传递 CDS 空间锚点；留出参考排除贯穿补搜。
- 保留真实比对器来源，Exonerate 发现模型和终端变体不能冒充 Miniprot 独立验证。
- 增加搜索时限、重建预算、延期/失败导出，以及注释只用于事后评估的分阶段召回诊断。
- 仍为开发版：短局部窗口、远缘基因重建、候选筛选和多昆虫目泛化尚需实测。


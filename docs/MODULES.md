# 模块地图

2026-11-11 初版。只做归类，**不移动任何文件**：脚本路径写死在 slurm 脚本、`docs/review_deps.txt`、SPEC 和已完成作业的来源记录里（例如 V0 的 `v0_result.json` 按路径记录 248 的 sha256）。

归类依据是每个脚本文件头 docstring 的第一段，以及 README、`docs/PROJECT_SUMMARY.md` 和此前的结论记录，并非逐行读码。状态一栏只写有出处的结论。编号表示开发先后，不表示模块边界。

| 模块 | 一句话 | 脚本 | 状态 |
|---|---|---|---|
| M0 基础设施 | 共用库、环境、审查循环 | `scripts/lib/`、`lib_insert.py`、`tools/`、`env/`、`docs/` | 在用 |
| M1 结构发现（Layer 1–4） | 不用转座酶模型，从组装比较中发现插入等位、重建、定极性 | 20–99 中的发现、基准与审计脚本，`_retired/` 中的 Layer 3 | Layer 2 锁定，Layer 3 关闭 |
| M2 语料与物种扩展 | 选物种、选基因组、建 panel、下载、建语料 | 90、91、95、98、102、140–143、180、191、202 | 语料 `corpus_v2`：40,541 个基因组 |
| M3 事件与插入数据库 | 跨 panel 去重成事件，推断空位点，导出 v1 记录格式 | 100、101、103–110、197–201、203–218 | Insertion DB v1 已导出 |
| M4 蛋白 / CDS 与 bag | 插入序列的 CDS、Pfam、按最长 ORF 聚成 bag，扩展 bag 位点，蛋白上下文导出 | 111–113、132–135、150–163、190、192–196、220–232 | bag 规则已确认（最长 ORF，不改） |
| M5 全语料拷贝扫描（isref） | 用 bag 代表序列 blastn 扫全语料，得到拷贝与位点 | 233–242、248、249 | **进行中**（SPEC 项 0 完成，项 1 待做） |
| M6 靶位与背景分析 | 元件插在什么基因里，跨物种占据多少独立背景 | 170–172、181–186 | 分析性；183、184 的计数已撤回 |

---

## M0 基础设施

| 文件 | 作用 |
|---|---|
| `scripts/lib/fasta.py`、`paf.py`、`tsd.py`、`util.py` | FASTA 随机访问、PAF 读取、TSD 检测（stage 30/31 共用）、union-find 与区间合并 |
| `scripts/lib_insert.py` | 从 locus 行取插入序列的唯一入口（修正 `lcp_bp` 坐标系 bug，2026-09-24） |
| `scripts/local_record_objects.py` | 本地审计：候选区间组（非已验证事件） |
| `tools/codex-review`、`docs/SPEC.md`、`docs/review_deps.txt`、`docs/reviews/`、`docs/DECISIONS.md` | 审查循环（CLAUDE.md） |
| `env/` | conda 环境定义 |

## M1 结构发现（Layer 1–4）

架构：Layer 1 发现（Arm A / Arm B）→ Layer 2 等位重建 → Layer 3 极性 → Layer 4 插入前靶位 → Layer 5 注释（最后才做）。

| 脚本 | 层 | 作用 | 状态 |
|---|---|---|---|
| `_retired/10_skani_cluster.py` | L1 | skani ANI 聚类（已由 M2 的 panel 取代） | 退役 |
| `20_armA_multicopy.py` | L1 | Arm A：基因组内 ≥3 拷贝的重复 | 参数锁定（300 ≤ L < 5000，min_copies 3）；minimap2 `-f 0` 修复 |
| `30_armB_graph.py` | L1 | Arm B：ANI 簇上 minigraph 图泡 → 空/满等位 | 在用；>5 kb 泡另存不丢 |
| `_retired/30b_armB_pangraph.py`、`_retired/31_armB_refine.py` | L1 | PanGraph 引擎、核苷酸级精修 | 退役 |
| `32_nested_search.py` | L1 | 在 >5 kb 泡内搜嵌套的 <5 kb 插入 | 已运行（E. coli、K. pneu）；审阅后**结果不可用**，见 `docs/STAGE32_REVIEW.md` |
| `_retired/40_merge_catalog.py` | L1 | 两臂合并 | 退役（由 86 取代） |
| `50_annotate.py` | L5 | 注释，最后做，按独立证据轴报告 | 已运行 |
| `60`–`64` | 基准 | 60 Arm A vs ISEScan 列联表；61 读段金标准分层（从未在真实数据上运行）；62 A∪B 召回；63 Arm A 漏检分类；64 重复惩罚诊断 | 结论见 README |
| `70_allele_reconstructor.py` | L2 | 锚点放置、精确与容错分解 | **锁定 2026-08-28** |
| `71`–`73` | L2 | 71 与 Arm B 对照；72 不一致审计；73 冻结回归（8 个门） | 73 是锁定门 |
| `_retired/80`–`85`、`88`、`89` | L3 | 外群选择、等位观察、敏感性、局部谱系极性、ingroup 精炼、汇总、锚点可比对性、回归 | **关闭 2026-08-29**；P1 层 |
| `_retired/86_pre_insertion_target_catalog.py` | L4 | 只取 P1 的插入前靶位目录 | 被 86_catalogue 取代 |
| `_retired/92_structural_polarity.py` | L3 | 仅凭 TSD 定极性（两基因组） | 退役 |
| `86_catalogue.py` | 产品 | 每 locus / 每事件一行，加序列 | 在用 |
| `87_novelty_positional.py` | 注释 | 按位置判断插入是否与已知 IS 重叠 | 在用 |
| `96`、`97` | 分析 | 同一锚定位点上的两个独立插入，及两个零假设 | 已关闭：1 例干净、1 例边缘 |
| `99_mobility_join.py` | 证据 | Layer 2 插入与 Arm A 移动性证据的连接（LINE M） | 已运行 |
| `120`–`123` | 审计 | 复现并修复 `lcp_bp` 坐标系问题（86/100 切长等位）与反向互补 | 已修复，见 `lib_insert.py` |
| `124_catalogue_summary.py` | 汇总 | 各物种 catalogue 目录的汇总数字（result_note 引用的数） | 在用 |
| `125`、`129`–`131` | 审计 | 多等位：70 只比最长与最短，中间等位的信息与新位点 | 已运行 |
| `126`–`128` | 审计 | 跨批次"新事件"是新位点还是同一位点换了键 | 已运行 |

## M2 语料与物种扩展

| 脚本 | 作用 |
|---|---|
| `90_cluster_density_report.py` | 以采样拓扑（而非基因组数）衡量数据库容量 |
| `91_seed_panel_sampler.py` | 以种子为中心、两两可比的 panel |
| `95_species_census.py` | 各物种可用基因组普查（按可用数而非总数排序） |
| `98_species_ani_calibration.py` | 按物种自身分布定 ingroup ANI 下限 |
| `102_rarefaction.py` | 增加基因组是否还有收益（回溯稀疏化） |
| `140`–`143` | SPIRE MAG：可用量、可挖掘物种、下载列表、ANI 分布形状 |
| `180_queue_species_download.py` | 单物种下载队列（只取 Complete / Chromosome） |
| `191_make_corpus.py` | 拼接语料并记录 contig 长度 |
| `202_build_genome_manifest.py` | 已验证的 accession → 基因组路径映射 |

## M3 事件与插入数据库

| 脚本 | 作用 |
|---|---|
| `100_event_dedup.py` | 把各 panel 的 locus 合并为不同的事件（键须与 junction 无关） |
| `101_layer12_validity.py` | 无冻结基线物种上的 Layer 1/2 有效性 |
| `103`–`105`、`107` | LINE M：M1 门诊断、Arm A 缺口、物种内存在矩阵、基因组内拷贝结构 |
| `106_element_families.py` | 跨物种元件家族（按测得的阈值聚类） |
| `108`、`109` | "≥2 个不同上下文"的组成匹配零假设与定向零假设 |
| `110_build_database.py` | 每个（空靶位，插入）对一行，全部物种 |
| `197`–`199` | 插入位置表（空基因组坐标系）、单点查询、SQLite 索引 |
| `200_build_insertions_all.py`、`201_validate_positions.py` | 把序列、空位点、位置合到一个目录；在 FNA 中核对位置 |
| `203`–`212` | carrier 推断空位点：推断、合成测试、报告、断点等价修订、交叉表、审计 A、最终验收、保守子集 |
| `213`–`218` | v1 记录格式导出：前置条件、源清单、记录、验收、验收器的负面对照、dataset.json |

## M4 蛋白 / CDS 与 bag

| 脚本 | 作用 |
|---|---|
| `111`–`113` | 第一代 bag：CDS 聚类、Channel A/B bag、nc 重切检查 |
| `132`–`135` | bag 内位点独立性、ORF 数分布、最长 ORF 的脆弱性、bag 语料完整性 |
| `150`–`154` | 拷贝转为位点的潜力；切除重建空位点（151 循环论证已撤回，由 152 取代）；正确性预测 |
| `160`–`163` | 从拷贝切出 bag 位点、回贴核对、bags_v5、逐字段比较 |
| `190`、`192`–`196` | tblastn 蛋白上下文 [2 kb][CDS][2 kb]：抽取、复核、按基因组、全 bag 一次搜索、分拆、JSON |
| `220`–`226` | CDS 试点：选样、Prodigal 与六框两条路线、对照、Pfam、比较、上下文挽救、报告 |
| `227`–`229` | 主导 CDS（最长 ORF，已确认）聚类成 bag、展开成 occurrence、六框片段过滤 |
| `230`–`232` | 单位点 bag 是否像转座酶、为何只出现一次、按结构（而非域名）判断是否为 IS |

## M5 全语料拷贝扫描（isref）

受 `docs/SPEC.md` 管理，进度与遗留状态见 `docs/legacy_shard_inventory.md`。

| 脚本 / 作业 | 作用 |
|---|---|
| `233_isref_copies.py`、`235_copy_flanks.py` | 早期单查询版：扫描结果转拷贝，重建拷贝周围的空等位 |
| `234_fna_to_gbk.py` | FNA 转 GenBank（供 ISCompare） |
| `236`–`238` | 查询集选择、跨 bag 等价记录、可疑 bag 分组 |
| `239_chain_copies.py` | 多查询 HSP 链接成拷贝 |
| `240_verify_shards.py` | 退役（SPEC r2） |
| `241_audit_shards.py` | 附录 A 的探索性统计（不是验收项） |
| `242_collapse_sites.py` | SPEC 项 3：拷贝合并为位点（待按 SPEC 重写） |
| `248_trusted_manifest.py` | SPEC 项 0：V0 可信清单（已完成，job 27788985） |
| `249_io_probe.py` | Lustre I/O 探测（诊断用） |
| 待写：243、244、246 | SPEC 项 1 的 outfmt 7 记账与发布，项 2 的遗留对照 |
| `slurm/isref_v1/`、`slurm/isref_all/` | 作业脚本；`C_sweep.sh`、`C_retry.sh` 为遗留扫描 |

## M6 靶位与背景分析

| 脚本 | 作用 |
|---|---|
| `170_preinsertion_gene_context.py` | 靶位在插入前是不是基因 |
| `171`、`172` | IS 元件插在哪些基因里；元件类型 × 靶基因类型 |
| `181_cross_species_expansion.py` | 候选元件跨全部物种的背景 |
| `182_background_breadth.py` | 每个元件家族占据的独立靶背景数 |
| `183`–`186` | 挖掘方法比较（计数被 184 撤回）、四级背景（被 185 撤回）、背景审计、记录级四对象 |

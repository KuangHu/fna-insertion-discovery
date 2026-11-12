# 基因组语料（M2）

2026-11-12。规范见 `docs/SPEC.md` 项 M2-C；代码 `scripts/250_build_corpus.py`。`$W` = `/global/scratch/users/kh36969/fna_ins_discovery`。

## 语料表是什么

`$W/isref_v1/corpus_v2.tsv`：40,541 行（`assembly`、`species`、`fna_path`），sha256 `b33819c772c522cc066892589ac97c7845614f49579ff250ec8b55ae585ab28e`。M5 全语料扫描（SPEC 项 0–3）的 shard 划分、V0 可信清单和遗留扫描都以它为准，**不改动**。

它是 2026-11-09 用一次性脚本生成的，那个脚本已经不在了。现在 `250 build` 按以下规则重建它，与原文件逐字节相同（2026-11-12 实测）：

- 来源是 `config/corpus_sources.tsv` 中的 19 个目录，每个目录的文件数固定写在配置里；目录少了或多了文件都会失败，而不是悄悄产生另一份语料。
- 取每个目录顶层的全部 `*.fna`，不进入子目录（`_plasmid_only/` 中是有意搁置的只含质粒的"完整基因组"，S. enterica 55、P. aeruginosa 7、L. monocytogenes 2、B. pertussis 1）。
- 物种标签取目录的标签，路径不解析符号链接（E. coli 的 13,027 个文件是指向另一个项目暂存区的符号链接），按 accession 排序。

## 重复：GCA 与 GCF

E. coli 的来源目录同时含有 GenBank（GCA）与 RefSeq（GCF）两份记录：7,723 个 GCA，5,304 个 GCF，每个 GCF 都与某个 GCA 同一 9 位编号。全量比较（`M2_corpus.sh`，见下）：5,298 对序列完全相同；5 对是同一组装的不同版本（语料中的 GCF 是 `.1`，GCA 是 `.2`，NCBI 汇总中与 GCA `.2` 配对的是 GCF `.3`）；1 对同版本但内容不同（GCF_002925525.1 只含 GCA_002925525.1 7 条记录中的 5 条，少了 46,701 与 18,617 bp 两条；NCBI 汇总也标为 `different`）。其他物种没有 GCF。所以语料中实际是约 35,237 个不同组装，不是 40,541。

M1 的 panel 不受影响：23,664 个 panel 清单中，没有一个同时含一对中的两个成员。

用户决定（2026-11-12）：语料表不改，由 `250 pairs` 产生一份映射 `corpus_duplicates.tsv`，M5 按分离株计数时据此把 `collapse = yes` 的成员并入代表（GCA 优先，其次新版本）。见 SPEC 项 3 的 S-0。按此规则 `collapse = yes` 5,303 个（相同 5,298 + 不同版本 5），合并后 **35,238** 个；GCF_002925525.1 按规则为 `differs`、不合并（它与 GCA 是同一分离株、RefSeq 少收两条记录；是否例外合并由用户决定）。

`250 pairs` 读取语料中的每一个基因组（读入字节数须等于文件大小），在 `corpus_content.tsv` 中逐个登记路径、字节数、文件 sha256 与序列内容摘要。它另外列出**不同编号、序列完全相同**的组装（`corpus_identical_other_id.tsv`），只报告、不合并：克隆性物种（B. pertussis、T. pallidum、M. pneumoniae 等）的不同分离株可以序列相同；也有同一份数据重复提交的情况（如 E. coli K-12 MG1655 的 GCA_000005845.2 与 GCA_004358965.1）。合不合并由用户另行决定。全量结果：125 组、344 个组装（含 GCF 成员），涉及 275 个不同编号；若全部合并，再少 150 个编号。没有跨物种的组。E. coli 33 组 145 个，K. pneumoniae 20 组，L. monocytogenes 13 组，M. tuberculosis 11 组，C. trachomatis 10 组，其余 ≤ 8 组。最大一组是 E. coli GCA_030345415.1–GCA_030345635.1 的 12 个连续编号（各一条 4,687,476 bp 染色体，序列相同）。

## 缺口

`250 gaps` 对照 2026-11-12 下载的 GenBank 细菌汇总（`$W/ncbi_meta/assembly_summary_genbank_bacteria.20261112.txt`，sha256 在同名 `.sha256` 文件），以 2026-09-01（物种普查的日期）为截止日：

| | 数目 |
|---|---:|
| 截止日前已发布、语料中没有 | 1,540（E. coli 1,243，M. tuberculosis 213，E. faecalis 35，其余 ≤ 12） |
| 截止日后发布 | 3,342 |
| 语料中有、汇总的可用集中已没有 | 6（S. aureus 4，E. faecalis 2） |
| 同编号但语料中的版本较旧 | 40 |

E. coli 的 1,243 个即早先已知的"Step 1 增量"（约 1,238 个），从未下载。其余 297 个是普查之后才进入汇总、发布日期早于截止日的组装。

`gaps` 同时写出下载队列 `corpus_download_queue.tsv`（只含截止日前缺失、有 NCBI 路径的行）与搁置清单 `corpus_plasmid_only.tsv`（65 行；同编号的多个版本都列出、都计数，实际数据中没有）。

截止日前缺失的 1,540 个已于 2026-11-12 下载到暂存目录 `$W/corpus_staging/20261112/<species>/`（`tools/download_staging.sh`，在登录节点运行，用户批准），每个都通过 NCBI md5 与 `gunzip -t` 校验，0 个失败，共 7.1 GB。审查第 3 轮之后工具收紧（md5 清单单独下载并检查 curl 退出码，URL 须有主机名），用新工具重新核对，1,540 个全部通过，未重新下载（`download.recheck3.log`）；第 4 轮后 `.ok` 改为逐字段核对，再次复核全部通过（`download.recheck4.log`）。**它们不在语料中**；并入语料需要用户另行决定，并修订 SPEC 项 0–3（40,541 写在其中）。截止日后的 3,342 个尚未下载。

## 全量运行

作业 27807675（2026-11-12，n0012.cf1，提交 `8aba32f`，26 分钟），结果在 `$W/corpus_m2/runs/27807675/`：`build` 与 `corpus_v2.tsv` 逐字节相同；`gaps` 与登录节点上的 dev 运行逐字节相同；`pairs` 读取 40,541 个基因组共 179.7 GB，`corpus_content.tsv` 逐个登记。三个子目录都是 `status: complete`，产物 sha256 与结果 JSON 一致。运行记录中 `git_dirty: true`，原因是仓库中有与 M2 无关的未跟踪文件（M5 的 241、242、V1、V2 脚本）；作业脚本已核对 250、作业脚本本身与配置与提交逐字节相同。M5 消费映射时按 SPEC S-0 核对 `pairs.json` 的 `corpus_sha256` 与产物 sha256。

## 运行

全量运行用 `slurm/corpus/M2_corpus.sh`（build、gaps、pairs 依次运行，输入以 sha256 固定）。登录节点上只跑测试 `tests/test_build_corpus.py`，以及用小样本语料表（几十个基因组，sha256 照样固定）运行的 `pairs`。

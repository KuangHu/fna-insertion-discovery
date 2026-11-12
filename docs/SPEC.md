# SPEC：全语料扫描 —— V2 shard 验证 / 241 遗留 shard 审计 / 242 位点合并

状态：修订版 r6。第 1 轮审查见运行 r01（FAIL，16 条）；运行 r02 因 API 断线 ERROR，不计轮数；第 2 轮见运行 r03（FAIL，7 条）；第 3 轮见运行 r04（FAIL，5 条）。第 3 轮后达到 3 轮上限，用户 2026-11-11 决定并入第 3 轮的 5 条意见、再审 1 轮（第 4 轮，运行 r07，FAIL，2 条）。用户随后特许第 5 轮，范围限定见项 P。用户审定前不提交，也不进入项 1–3 的审查循环。

本文件每一条验收标准都写成可检查的条款。审查者按"当前项"审对应一节，外加第 0 节与"流程"。条款末尾的 `[R1-n]` 表示该条处理了第 1 轮（运行 r01）的第 n 条意见，`[R2-n]` 表示第 2 轮（运行 r03）的第 n 条，`[R3-n]` 表示第 3 轮（运行 r04）的第 n 条，`[R4-n]` 表示第 4 轮（运行 r07）的第 n 条，`[R5-1]` 表示第 5 轮（运行 r08）的唯一一条；`[用户 n]` 表示用户 2026-11-11 第一次决定的第 n 点，`[用户二 n]` 表示第二次决定的第 n 点。

**审查材料。** 每一项的必需材料在 `docs/review_deps.txt` 中按项声明，由 `codex-review`（v4，仓库副本 `tools/codex-review`）逐个检查：存在、非空、文本、不超过 200,000 字节、未被 `.gitignore` 忽略。声明的文件不论已提交、已修改还是未跟踪，都全文发送；以 `DIFF <旧> <新>` 声明的一对文件，两侧都同样检查，但只发送 `diff -u`。任一不满足，wrapper 在调用审查者之前直接输出 `VERDICT: ERROR`。审查者若发现被审内容引用了清单之外、未提供的文件，也判 ERROR，不得猜测其内容。

## 修订记录

- r6.3（用户 2026-11-11 第十次决定）：新增**项 M1-B：Arm B 结构等位发现**（M1 模块，见 `docs/MODULES.md`），用 gpt-6-astra 审查后按意见修改。只新增此项，项 0–3 不变。
- r6.2（用户 2026-11-11 第九次决定：先修改再重跑，未送审）：V0-7 新增"逐次记录"。起因：job 27788052 的失败记录中出错基因组列表为空（父进程达到上限即退出循环，出错 worker 的结果未收回）。放宽或删除：无。
- r6.1（项 0 第 5 轮之前；用户 2026-11-11 第八次决定第 2、3 点）：只改项 0。新增 V0-7（读取基因组时对瞬时 I/O 错误有限重读），V0-6 补充相应测试，项 0 写明第 5 轮的范围。起因：共享 scratch 多次在读取中途报 `Input/output error` / `Cannot send after transport endpoint shutdown`（`docs/lrc_io_incident.md`）。放宽或删除：无。重读不放宽任何内容检查：每次重读都是完整的重新读取，原有的全部检查照常执行，另加两次完整读取一致的要求。
- r6（本版；处理第 4 轮的 2 条意见，全部接受；第 5 轮由用户特许，范围限定）：
  - S-5a COMPONENT (ii)：不再用跨度净差。改为分别约束 h 内两个方向的缺口列总数：`length − query 跨度`（subject 多出、query 一侧为缺口的列数）与 `length − subject 跨度`（另一方向），各自 ≤ 50 bp。插入与删除不能再互相抵消。r5 修订记录中"不会多吸收"的说法改为下文 S-5a 中有依据的陈述。S-17 新增抵消反例 (a3c) [R4-1]。
  - V2-6 新增 (g) 坐标合法性：查询坐标在 `[1, qlen]` 内，subject 坐标经 min/max 处理负链后在 `[1, slen]` 内，`slen` 等于可信清单中该 contig 的长度，`length` 不小于两侧跨度，失配加缺口列不超过 `length`。新增负面测试与一条合法负链行的正面对照 [R4-2]。
  - 连带改动：V0-1 与 V2-4(c)(d) 的基因组可信清单新增每条 contig 的长度，**V0 清单须重建**；A-D1a 新增遗留一侧的坐标检查（以遗留行自报的 `qlen`、`slen` 为界），`slen` 与可信长度不符与 `qlen` 一样单独计数、不判畸形；S-3 新增 `s0 ≥ 0` 与 `s1 ≤ contig 长度` [R4-2]。
  - 事实核对：遗留与 V2 的输出格式都含 `qlen`、`slen` 两列（`C_sweep.sh:67`、`V2_verified.sh:125`），两处的长度过滤都保留全部 14 列（`C_sweep.sh:72`、`V2_verified.sh:147`）。
  - 第 5 轮（运行 r08）唯一一条意见：S-17 (a3) 仍用旧的"长度差"写法，与 (ii)、(a3c) 冲突；已按审查者的建议改写 (a3) [R5-1]。此项修改未再送审（用户限定不开第 6 轮），由用户定稿时确认。
  - 放宽或删除：无。
- r5（处理第 3 轮的 5 条意见，全部接受；第 4 轮由用户特许）：
  - G-9 区分**数据产物**与**结果 JSON**：三方摘要相等只适用于数据产物；结果 JSON 的路径、字节数、sha256 只由提交记录绑定，不要求它包含自身的摘要。V2-9、V2-12、S-16 相应改写 [R3-1]。
  - S-5a E2：去掉"块内查询区间并集"的写法。改为按 subject 顺序检查整条链的查询单调性，并用实际接触 L1、R0 的两条边界 HSP（必须唯一）核对查询间隙；有歧义即不可评估、不升级。S-17 加入审查者给出的反例 [R3-2]。
  - S-5a COMPONENT (ii)：去掉无法唯一定义的坐标投射，改为只用 subject 区间与端点即可机械检查的条件（外层链中单条 HSP 连续覆盖 b，且该 HSP 内没有大的长度差）。不可评估时 b 保持候选 [R3-3]。
  - A-D：新增 A-D1a 逐行格式与身份校验，先于比较；新增资格类判词 `LEGACY_MALFORMED` 及计数；V2 一侧校验失败判 `NO_V2_TWIN`（原因 `V2_ROWS_INVALID`）；idx 在 A-D 内定义；禁止跳过异常行后继续比较。A-D8 加入未知 id、残缺末行测试 [R3-4]。
  - G-3、S-9：方向共识按 bag 分别给出；位点涉及多个已知方向的 bag 时，位点方向取 `UNKNOWN`；`MIXED` 只用于同一 bag 内部。S-17 加入跨 bag 方向相反的测试 [R3-5]。
  - 项 P 的 P-8 改为对照上一轮审查 `r04.md` 与差异 `r04.SPEC.md → SPEC.md`（理由见 P-8；早先各轮已由其后一轮核对）。审查输入随之不再整份包含 r01/r03/r04 的 SPEC 快照与 r01/r03 审查文本。
  - 放宽或删除：无。COMPONENT 判定 (ii) 被替换而不是删除：新条件同样是吸收 b 的必要条件，缺少时 b 保持候选，不会因此多吸收。（r6 更正：后半句没有依据，第 4 轮已指出，见 r6 记录与 S-5a。）
- r4（处理第 2 轮的 7 条意见，全部接受）：
  - G-9 增加"发布前字节一致"与"消费时三方摘要相等"（结果 JSON 内的产物摘要 = 提交记录中的摘要 = 当前文件重算的摘要），写明字段映射 [R2-1]。
  - 新增**项 0：V0 可信清单**，并改写流程：项 0 审查通过后，经用户批准运行 V0、提交可信配置，然后才进入项 1 的完整审查，解除循环依赖 [R2-2]。
  - A-D3 规定判定顺序与唯一边界 b（第一个差异查询），前缀模式要求 b 处遗留为 V2 的真子集；相应修正 A-D8 [R2-3]。
  - 超长记录解析后成为"被打断锚"，参与单元构建，多条查询解析出的同一外层拷贝只计一次；S-6 写明固定的处理顺序 [R2-4]。
  - E2 绑定到外层记录自身的 HSP 链（用 239 的 `chain()` 重建），写明方向、共线与坐标公式（含负链）[R2-5]。
  - E1：覆盖 ≥ 0.50 只表示"未满足 E1"，候选保持候选；"组成部分"改为须有正面证据才能吸收 [R2-6]。
  - S-1、S-16 恢复"未消费 shard 及原因"，原因取值固定 [R2-7]。
  - 放宽或删除：无。
- r3（在 r02 ERROR 之后，按用户第二次决定修订）：
  - A-S 移出所有验收要求，改为附录 A（可选的探索性分析或历史记录）；项 2 的验收内容只剩 A-D [用户二 4]。
  - A-D 的判词全部改为观察（`*_OBSERVED`），新增 A-D4"观察与结论必须分开"；新增 `LEGACY_OUTPUT_DELETED`、`LEGACY_OUTPUT_ABSENT`；限制依据新的 64 shard 清单 `docs/legacy_shard_inventory.md` 更新：仍能直接对照的短 shard 是 58 与 63 [用户二 4, 用户二 6]。
  - "不同 bag + 包含度 ≥0.90"降为嵌套**候选**。新增 S-5a：候选要升级为单独计数的嵌套位点，必须同时有元件不同（E1）与插入边界（E2）的证据；超长记录在同样证据下，可以解析为被打断的外层拷贝 [用户二 5]。
  - 放宽或删除：A-S 不再是验收要求（用户决定）。其他条款都没有放宽。
- r2：
  - 第 1 轮 16 条全部接受，逐条落实于标注 `[R1-n]` 的条款。
  - 用户决定：路线 B（全部 64 个 shard 最终都须 `VERIFIED_V2`）；位点合并原则（防止串联桥接、成对重叠不作生物学证据、片段不另成位点、保留独立支持的嵌套插入、无法解析者标为不确定）；保留链方向与查询方向；240 退役；审查基础设施 v2。
  - 结构性改动：
    - 新增 G-9 统一发布协议，以提交记录取代 done marker；
    - 241 拆为 A-D（遗留与 V2 的确定性直接比较）和 A-S（非门控的统计迹象筛查），A-S 的参考集只取 `VERIFIED_V2`；
    - 242 由单一的双向 0.50 规则改为锚 / 支持 / 嵌套 / 不确定结构的分级模型，阈值 0.80 / 0.90。
  - 放宽或删除了什么、为什么：
    - 删除 A-6 的概率上界 10⁻⁵ 与族错误率 6.4×10⁻⁴：它们不成立 [R1-6]；
    - 删除"失败时正式目录不得留下文件"：节点崩溃时无法兑现 [R1-5]，改为更强的"未提交即不可见且不可消费"；
    - 241 不再为生产门控：路线 B 下生产只消费 `VERIFIED_V2`；
    - 242 不再使用 0.50：它无法区分"同一拷贝"与"嵌套"，见 S-5。

---

## 项 P：SPEC pre-review

审查命令：`codex-review "SPEC pre-review" review-hard`

**本项的审查对象是 SPEC 本身，不是代码。** 现有代码不满足某条标准，不构成本项 FAIL 的理由，代码在项 1–3 中审查。但代码与其他材料可以用来核对 SPEC 中的事实陈述：默认值、参数来源、文件格式、行号。

**判定标准**
- **P-1 正确性**：每条标准在技术上成立；SPEC 所述事实与材料一致。
- **P-2 一致性**：条款之间、项与项之间无矛盾，尤其是文件命名、命名空间、提交记录、门控条件在 V0 / V2 / 241 / 242 之间能否对上。
- **P-3 可检查性**：每条标准都能由测试，或对输出、源码的机械核对来判定。不能检查的标准须指出。
- **P-4 统计假设**：指出没有依据或可能不成立的假设，并说明假设不成立时误判的方向。
- **P-5 位点合并**：指出会虚增或丢失位点的规则，包括片段、长度变体、串联、嵌套、异链、跨原点。
- **P-6 阈值来源**：每个数值阈值是否按 G-8 标注，所写来源是否真的支持该数值。
- **P-7 完整性判定**：是否存在任何路径，能在证据不足时把 shard 或位点标为完成或已验证。
- **P-8 不得削弱**：对照最近两轮已完成的审查（`r04.md`、`r07.md`）及第 3 轮被审 SPEC 到本版的差异（`diff -u r04.SPEC.md SPEC.md`，由 wrapper 生成）检查：
  - 这两轮每条被接受的意见是否确已处理（第 4 轮被审的 SPEC r5 与 r4 的差别，以及 r6 与 r5 的差别，都包含在这份差异中，由修订记录区分）；
  - 差异中任何验收标准被放宽、删除或改为不可检查，而修订记录没有给出理由的，判 FAIL。

  更早的轮次已由其后一轮核对过（r03 核对 r01 的意见与 r01→r03 的改动，r04 核对 r03 的意见与 r03→r04 的改动）；r04 之后对任何条款的改动，包括撤回早先的修正，都必然出现在这份差异中。因此本轮只发送这一份差异，不再整份发送历次 SPEC 快照与更早的审查文本 [用户 2026-11-11 第四次决定]。

**第 5 轮的审查范围** [用户 2026-11-11 第五次决定]。第 5 轮只审 r6 的两处修改及其连带改动：S-5a COMPONENT (ii) 与 S-17 (a3c)；V2-6(g) 及其测试；V0-1、V2-4(c)(d)、A-D1a、S-3 中由此引起的改动；以及相应的修订记录。P-1 至 P-8 照常适用于这些条款，P-8 以差异 `r04.SPEC.md → SPEC.md` 与第 4 轮审查 `r07.md` 为依据。审查者若发现与上述改动无关的问题，仍请列出，但在条目开头标 `[范围外]`；这类问题由实现者记入 `docs/DECISIONS.md` 标为 deferred，不开第 6 轮。

**本项输出要求**：问题按"建议的 SPEC 修改"给出，每条注明条款编号、问题、为什么是问题、建议改法。只有 P-1 至 P-8 都没有实质问题时才判 PASS。

---

## 0. 通用约定（三项共用）

**G-1 计数与实体。**
- copy record：239 输出的一行，即一条查询在一个 contig 上的一条 HSP 链。
- 位点（site）：242 判定的一个插入拷贝所在的基因组位置，等级见 S-8。
- events（独立转座事件）：需要树上的极性，本 SPEC 三项都不产生。

任何计数字段名都须标明属于哪一种，不允许出现没有限定词的 "insertions" 或 "hits" 总数。

**G-2 坐标** [R1-16]。
- 每条 HSP 先取 `lo = min(sstart, send)`、`hi = max(sstart, send)`。
- 一条链的区间为 `[min(lo) − 1, max(hi))`，0-based 半开。这与 239 逐条 HSP 交换端点后再取包络的做法一致。
- 长度 `s1 − s0` 必须 > 0。

**G-3 方向** [用户 4]。
- `blast_strand`：copy record 的 BLAST 链，即 subject 相对该查询序列的方向。原样保留，任何阶段不得丢弃。
- `query_orient`：查询序列相对其 bag 参考序列的方向，取 `+`、`−` 或 `NA`。bag 参考序列即 bag 的原始 DNA 代表 `cluster_rep`。
  - 来源：A3 的 bag 内自比对 `selfaln_bag.tsv`。字段为 `qseqid sseqid pident length qstart qend sstart send qlen slen`（见 `A3_queries.sh` 的 outfmt），`sstart > send` 为反向。
  - 判定：在"查询 → 参考"与"参考 → 查询"两个方向的全部 HSP 中，取比对长度最长者的方向。
  - 以下情况取 `NA`：没有比对行；最长者并列且方向不同；两个方向的结论互相矛盾。
  - 不同 bag 之间没有共同参考；237 的等价表也不含链信息。因此跨 bag 的方向关系一律视为未知。
- `element_orient = blast_strand × query_orient`；`query_orient = NA` 时为 `NA`。`element_orient` 是相对**该 bag 参考序列**的方向，只在同一 bag 内可比 [R3-5]。
- 原始 `blast_strand` 相反**本身**不构成"不同拷贝"或"冲突"的证据，合并判定不得使用它。只有同一 bag、`element_orient` 均已知且相反时，才记为方向不一致（S-6）。

**G-4 可复现性。** 每个作业的每份输出（或其运行结果 JSON）必须记录：
- `git_commit_submit`：sbatch 时以 `--export=GIT_COMMIT=$(git rev-parse HEAD)` 传入。作业启动时该变量未设置即失败。
- `git_commit_runtime` 与 `git_dirty`：运行时实测。
- SLURM job id、array task id、主机名。
- 工具版本：运行时执行 `blastn -version`、`makeblastdb -version`、`python3 --version` 实测写入，不得硬编码。
- 完整参数串，以及全部输入文件的 sha256。

**G-5 继承的固定参数**（三项都不得改动；要改，另立 SPEC 项）：

| 参数 | 值 | 标注（G-8）与来源 |
|---|---|---|
| blast 任务 | `dc-megablast` | 实测：B2 基准（job 27721664）中 megablast 漏掉 4.35% 的位点，漏掉的是分化拷贝 |
| `-evalue` / `-perc_identity` / `-max_target_seqs` | `1e-30` / `80` / `1000000` | 继承：`C_sweep.sh`（job 27725308）；旧 181 查询扫描同值（`sweep.sh:46-47`） |
| 输出字段（14 个，顺序固定） | `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen` | 继承：`C_sweep.sh`；239 按此位置解析 |
| HSP 长度过滤 | 比对长度 `>= 200` bp | 继承：`isref_v1/jobs/sweep.sh:50`、`sweep2.sh:67`，由 `C_sweep.sh` 原样继承 |
| 链接与分级 | `full_qcov 0.90`、`part_qcov 0.50`、链接间隙 `<= qlen` | 继承：239 默认值，与 233 逐字节一致 |
| 编入目录判定 | 双向 reciprocal overlap `>= 0.50`，carrier 坐标 `insert_start0/insert_end0` | 继承：239 `--recip` 默认值，来自 233 |
| BLAST 版本 | 2.17.0 | 运行时实测（G-4）；`fnains_annot` 环境 conda-meta |

**G-6 数据不入库。** shard 输出、fna、blast db、日志都不进 git，由 `.gitignore` 排除。

**G-7 比较的类型规则** [R1-16]。
- 计数：两侧都必须是非空十进制整数。
- 摘要：64 位小写十六进制。
- 集合：先单独检查重复（除非条款另有规定，有重复即失败），再以排序去重后的列表比较。
- 浮点阈值：有限，且在声明的取值范围内。
- 字符串：逐字节相等。

任何一侧缺失、类型不符或越界，判失败（失败即关闭）。

**G-8 阈值标注** [R1-8, 用户 1]。每个数值阈值都标为以下之一，并写明依据：
- `继承`：沿用既有代码或作业的取值，写出处；
- `实测`：来自某次作业的测量，写 job id；
- `未校准操作值`：人为选定，未经校准。

SPEC 与任何输出都不得为任何阈值声称统计保证，包括错误率、置信度、p 值、概率上界。"与另一处取值一致"只说明一致，不说明经过校准。

**G-9 发布协议**（V2 与 242 共用；以提交记录取代 done marker）[R1-1, R1-5]。
- **运行目录**：每次运行把全部产物写入独立的 `<ns>/runs/<run_id>/`，其中 `run_id` 含单元名（如 `shard_17`）、`SLURM_JOB_ID` 与 array task id。下游永远不直接读取运行目录。
- **发布**：只有一种方式，即最后一步原子地（临时名 + `rename`）写出提交记录 `<ns>/committed/<unit>.json`。提交记录包含：
  - 单元名、`run_id`、UTC 时间；
  - 每个**数据产物**在运行目录中的相对路径、字节数、sha256；
  - **运行结果 JSON** 在运行目录中的相对路径、字节数、sha256（单独字段，不在数据产物列表中）；
  - 可信配置 `config/isref_all_trusted.tsv` 的 sha256；
  - `record_sha256`：对其余全部字段按规范化 JSON（键排序、UTF-8、无多余空白）计算。
- **两类文件** [R3-1]：
  - **数据产物**：下游消费的数据文件。它们的摘要与字节数同时出现在结果 JSON（按字段映射）和提交记录中。
  - **运行结果 JSON**：记录运行参数、记账与数据产物摘要。它**不**声明自身的摘要或字节数；它的摘要与字节数只出现在提交记录中。
- **发布前字节一致** [R2-1]：运行目录中的每个数据产物，必须与通过全部检查的那份字节完全相同。节点本地通过检查的文件复制进共享盘的运行目录后，重新计算 sha256 与字节数，必须等于检查时的记录；结果 JSON 中的数据产物摘要，取自通过检查的那份字节。结果 JSON 写入运行目录后，重算其 sha256 与字节数，写入提交记录。任何一处不等即失败，不写提交记录。
- **消费前核对**：下游消费任何产物之前，必须依次核对整条链，任一不符即拒绝消费：
  1. 重算 `record_sha256`；
  2. 结果 JSON **两方相等**：当前文件重算的 sha256 与字节数 = 提交记录中的值 [R3-1]；
  3. 对每个数据产物，**三方摘要相等**：结果 JSON 中按字段映射声明的摘要 = 提交记录中列出的摘要 = 当前文件重算的摘要；字节数同样三方相等 [R2-1]；
  4. 核对单元名与输入清单身份。

  字段映射见 V2-12（V2）与 S-16（242）。提交记录的数据产物列表与结果 JSON 按字段映射声明的数据产物集合必须完全相等：只在一侧出现的数据产物，都判核对失败。结果 JSON 不属于数据产物，不出现在数据产物列表中；提交记录缺少结果 JSON 的路径、字节数或 sha256，判核对失败。

  "提交记录 → 运行结果 JSON → 数据文件"这条核对链，是"已验证"或"已发布"身份的**唯一**来源。
- **未提交的运行目录**：
  - 对下游不可见；
  - 可识别：目录内的 `status` 文件为 `running` 或 `failed:<检查名>`；
  - 可回收：只在人工确认后删除，删除时不得触及任何被现行提交记录引用的目录。
- **失败**：不写提交记录。协议只承诺"不可见、不可消费"，不承诺"清理干净"。
- **替换**：同一单元的新提交记录以原子 `rename` 取代旧记录；被取代的运行目录保留，直到人工回收。

**G-10 审查前证据与退役组件** [用户 7]。
- V0（job 27769455）、V0b（27769564）、V1（27769473）的现有输出由未经审查、也未记录 commit 的代码产生，只作参考，不作验收证据。代码审查通过后，须用记录了 commit 的版本重跑。
- 以下组件退役，只保留溯源，不作任何验收依据：
  - `scripts/240_verify_shards.py`：HSP 跨 shard 中位数筛查。其历史输出（如 job 27725308 期间的 "OK 39 / SHORT 4"）不是验证结果。
  - `slurm/isref_all/C_retry.sh`：以"不同查询 id 数 ≥ 55,000"为门槛，已被否决。它在每个 shard 开始时就删除遗留输出（见 A-D5 与 `docs/legacy_shard_inventory.md`）。
- `slurm/isref_all/C_sweep.sh` 是遗留协议，只作溯源与比较基准。
- 可检查：V2、241、242 的代码及其驱动作业中，不得调用 `240_verify_shards.py`、`C_retry.sh`（grep）。

---

## 流程（路线 B；用户决定，不是验收条款）[用户 2]

1. 生产路径：全部 64 个 shard 最终都必须是 `VERIFIED_V2`（V2-14）。在此之前，任何位点计数都只能作为"部分结果"报告，并按 S-1 列出未纳入的 shard 及原因，不得作为最终结果。
2. 准备阶段（项 0）[R2-2]：
   1. V0 代码先按项 0 实现、做小样本测试、通过审查；
   2. 经用户明确批准，用 SLURM 运行 V0（这是项 1–3 审查完成之前唯一允许的 SLURM 操作）；
   3. 用户审阅 V0 输出后，提交 `config/isref_all_trusted.tsv` 入库。
3. 顺序：可信配置入库之后，项 1 → 项 2 → 项 3 依次实现、做小样本测试、通过审查。之后提出试点方案与修订后的资源估算，经用户批准才能 sbatch。
4. 重跑某个 shard 不会自动使它成为 `VERIFIED_V2`，只有满足 V2-14 的核对链才算。
5. 审查循环中禁止 sbatch、scancel，也禁止修改在运行作业所用的脚本（CLAUDE.md）。

---

## 项 0：V0 可信清单 [R2-2]

审查命令：`codex-review "V0 可信清单" review-hard`

**范围**
- `slurm/isref_all/V0_prep.sh`
- `slurm/isref_all/V0b_dbtest.sh`
- `scripts/248_trusted_manifest.py` 与 `tests/test_trusted_manifest.py`（待写：生成并自检可信清单）

**目的**：生成 V2 所依赖的可信期望值，让项 1 的审查与运行有一个先于它、单独审查过的锚点。

**验收标准**
- **V0-1 生成内容。** V0 生成 `config/isref_all_trusted.tsv` 的候选内容：
  - V2-2 表中的全部量；已有字面值的四项必须逐值相等，否则失败；
  - `all_queries.tsv` 的 sha256；
  - 基因组可信清单的 sha256。基因组可信清单含每个基因组的文件 sha256、contig 数、碱基数、contig id 列表的 sha256，以及**每条 contig 的完整 id（`ACC|CONTIG`）与长度**（V2-4c），存放在 scratch，由入库的 sha256 固定 [R4-2]。
  - **V0 清单须重建**：既有的 V0 输出（job 27769455）没有每条 contig 的长度，不能满足本条；按 G-10 本来也须用审查后的代码重跑。contig 长度表的 sha256 并入基因组可信清单的 sha256，不另设入库字段。每条 contig 的长度之和必须等于该基因组的碱基数，否则失败。
- **V0-2 清单划分。** 按 V2-4(a) 生成 64 份 shard 清单，检查两两不相交、并集等于 40,541、规模为 35×633 + 29×634；不满足即失败。
- **V0-3 不信任旧输出。** 所有值都从源文件重新计算，不读取任何先前 V0 运行的输出。
- **V0-4 工具行为自检。** outfmt 7 的正面记账测试（1 条命中、2 条不命中，必须全部宣告，尾行为 `processed 3 queries`），以及 V0b 的 `%t` / `%a` 测试，都在本次运行中重做，结果写入输出。
- **V0-5 可复现。** 输出记录 G-4 全部字段。可信配置在用户审阅之后才入库；入库前，项 1 的审查必然因缺少该文件而 ERROR（这是预期行为）。
- **V0-6 测试**（登录节点小样本）：248 对缺失文件、重复 id、清单相交、计数不符、contig 长度之和与碱基数不符都必须失败。V0-7 的测试用注入的错误模拟：一次瞬时错误后恢复（结果与无错误时逐字节相同，重读次数被记录）；持续错误直至该基因组失败；非重读类错误不重读；确认读取不一致即失败；节点累计达到上限即整体失败，失败记录写明节点名；`genome_reads.tsv` 被篡改时 `check` 失败。
- **V0-7 瞬时 I/O 错误的有限重读**（r6.1，用户第八次决定第 2 点）。
  - 只对读取基因组 FASTA 时的 `OSError` 中 errno 为 `EIO`、`ENOTCONN`、`ESHUTDOWN` 的三种重读；其他任何错误（文件不存在、权限、内容规则、读取字节数不足等）立即失败，不重读。
  - 每次读取都从第 0 字节开始完整重读：重新计算 sha256、重新解析，不保留失败那次读取的任何内容。每次读取都照常核对"读到的字节数 = 打开时的文件大小"以及全部内容规则。
  - 每个基因组最多重读 3 次（总计最多 4 次读取）。重读前的等待随错误次数递增：第 1 次错误后 10 秒，第 2 次错误后 30 秒。
  - **确认读取。** 出现过错误的基因组，必须再有两次完整读取在字节数、sha256 和每条 (contig id, 长度) 上完全一致才被接受，否则判 `genome_reread_mismatch` 失败。确认读取计入 3 次重读，所以一个基因组最多容忍 2 次错误，且都发生在第一次完整读取之前；某次错误之后，若剩余读取次数不足 2 次（一次完整读取加一次确认读取），该基因组立即失败（`genome_io_retries_exhausted`），不再等待和读取。
  - **节点上限。** 同一节点（V0 是单节点作业，即本次运行）累计 I/O 错误达到 5 次，整体失败，不再开始新的读取；失败记录 `FAILED` 写明节点名与 SLURM 作业号。任何原因的失败，`FAILED` 都写明节点名：248 写 `manifest/FAILED`；`V0_prep.sh` 在建立运行目录之后的任何非零退出（版本命令、自检、V0b、248）都写 `$RUN/FAILED`，含失败的行与命令、节点名、SLURM 作业号 [第 5 轮第 2 条]。
  - **记录。** 每个基因组的重读次数与错误类型写入 `genome_reads.tsv`（列 `shard assembly rereads io_errors`，按语料顺序，每个基因组一行）。`v0_result.json` 的 `io_rereads` 记录策略参数、节点名、节点错误总数，并单独列出重读次数 > 0 的基因组。`check` 核对该文件：行与语料一一对应；错误类型只能是上述三种；无错误时重读 0 次，有错误时重读次数 = 错误次数 + 1；节点错误总数 < 5。对已完成的运行，`check` 另要求 `v0_result.json` 存在且 `verified` 为真，各输出文件（含 `genome_reads.tsv`）的字节数与 sha256 等于其中的记录，`io_rereads` 汇总等于 `genome_reads.tsv` 推出的值 [第 5 轮第 1 条]。
  - **逐次记录**（r6.2，用户第九次决定，未送审）：每次上述错误发生时立即记录一条（时间、节点、进程、基因组、路径、errno、第几次读取、节点累计数）：写一行 `IO_ERROR {json}` 到 stderr（SLURM 日志），并加入父进程可读的共享列表（即使进程池被提前终止也不丢失）。失败记录写明各 errno 的次数及逐条的 `基因组:errno@时间#第几次`；成功时全部事件写入 `v0_result.json` 的 `io_rereads.events`。
  - 该文件**不属于**基因组可信清单，不进入 `genome_manifest.index`：可信清单的内容与 I/O 经历无关，有无重读都逐字节相同。
  - V0 作业以 `--exclude=n0052.cf1` 提交（SLURM 节点名带 `.cf1` 后缀；`n0052` 被 sbatch 拒绝为无效节点名）（job 27786542 在该节点上 6 个文件全部读取失败）。
- **第 5 轮的审查范围**（用户第八次决定第 3 点，特许）：第 4 轮（运行 r04）意见的修复（`load_expect` 的严格类型检查）与 V0-7 的实现、测试和 `V0_prep.sh` 的改动。其他条款照常适用。

---

## 项 1：V2 shard 验证

审查命令：`codex-review "V2 shard 验证" review-hard`

**范围**（必需材料见 `docs/review_deps.txt`）
- `slurm/isref_all/V2_verified.sh`
- `scripts/243_outfmt7_accounting.py` 与 `tests/test_outfmt7_accounting.py`（待写，V2-6）
- `scripts/244_publish.py` 与 `tests/test_publish.py`（待写，G-9 的写出与核对；V2 与 242 共用）
- 依赖：项 0 产出、经用户审阅后入库的 `config/isref_all_trusted.tsv`（V2-2）

job 脚本已从 scratch 原样复制入库，与原件 `cmp` 逐字节一致。

**目的**：重跑一个 shard，使"已验证"由逐条记账与守恒等式判定，而不是由计数推断。

### 验收标准

**V2-1 禁止统计性计数门槛；精确守恒必须执行** [R1-2]。
- 不得以下列量作为接受条件：命中中出现的不同查询数、HSP 总数、跨 shard 中位数及其倍数。这些量可以写入记录作参考。
- V2-6、V2-8 规定的逐条记账与守恒等式是精确的身份检查，不属于此禁令，而且必须执行。

**V2-2 可信期望值入库。** 期望值从仓库内已提交的 `config/isref_all_trusted.tsv` 读取，不得从与数据同处 scratch 的文件读取。缺失、为空或类型不符（G-7）即失败。

| 量 | 期望值 | 标注与来源 |
|---|---|---|
| `all_queries.fna` sha256 | `8d5bbd5a99988d34d94f8117dc712a5e0ba23aafa40701dccea661355d1e504a` | 实测：V0 job 27769455 |
| 排序后查询 id 列表 sha256 | `45dd8d78c6a990be93252741f4c76592126ce3aac9605ea439d84c87ff423caa` | 实测：同上 |
| 查询条数 N | 63,095，无重复 | 实测：同上 |
| `all_queries.fna` 字节数 | 96,988,793 | 实测：job 27768584 暂存时 |
| `all_queries.tsv` sha256 | V0 重跑时计算后入库 | V0 [R1-3] |
| 基因组可信清单（V2-4c）sha256 | V0 重跑时计算后入库 | V0 [R1-4] |

项 0 的 V0 运行后（流程第 2 步），表中已有字面值的项必须逐值相等，否则停止。

**V2-3 查询序列与元数据一致** [R1-3]。节点本地暂存副本须同时满足：
- (a) FASTA：sha256、排序 id 列表 sha256、条数 N 都等于期望值；id 无重复。
- (b) `all_queries.tsv`：
  - sha256 等于期望值；
  - 表头逐字为 `query_id bag_id seq_key len_nt`；
  - `query_id` 唯一，且 `query_id = "Q_" + seq_key`；
  - `query_id` 集合等于 FASTA 的 id 集合；
  - 每条查询的 `len_nt` 等于 FASTA 中该序列的长度；
  - `bag_id` 非空。
- (c) 交给 239 的 `--queries` 文件，必须是通过 (b) 的同一份字节，以 sha256 核对。

bag 映射的来源（Stage A 的 bag 表）以 `all_queries.tsv` 的 sha256 固定，本项不重新推导 bag。

**V2-4 基因组清单与读入完整性** [R1-4]。
- (a) shard i 的清单：`corpus_v2.tsv` 去掉表头后，第 k 行（从 0 起）满足 `k mod 64 = i`。这与 `C_sweep.sh` 的 `(NR−2)%64` 一致。
- (b) 64 份清单两两不相交，并集等于 40,541 个 assembly，规模为 35×633 + 29×634（实测：V0）。由 V0 检查，不满足即失败。
- (c) V0 为每个基因组记录文件 sha256、contig 数、碱基数、contig id 列表与每条 contig 的长度，构成基因组可信清单；清单本身的 sha256 入库（V2-2）[R4-2]。
- (d) V2 对每个基因组只读取一次字节流：同一次读取同时送进 sha256 计算和写入 shard FASTA 的转换（例如 `tee`）。须核对：
  - sha256 等于期望值；
  - 本次转换输出的 contig 数、碱基数、contig id 列表、每条 contig 的长度都等于期望值；
  - 基因组内 contig id 唯一；
  - 派生 id `ACC|CONTIG` 在 shard 内唯一。
- (e) 用语：(d) 证明写入的内容来自可信字节，且计数和 id 完整；V2-5 证明数据库的序列 id 集合与计数完整。二者都不单独称为"内容验证"。

**V2-5 数据库的 id 集合与计数** [R1-4]。
- `blastdbcmd -info` 的序列数等于 Σ 期望 contig 数；残基数等于 Σ 期望碱基数。
- `blastdbcmd -entry all -outfmt %t` 得到的完整 id（`ACC|CONTIG`）集合，等于由可信清单推出的完整 id 集合。比较完整 id，不只比较 accession。
- `makeblastdb` 参数与 `C_sweep.sh` 完全一致，不加 `-parse_seqids`。依据：V0b（job 27769564）实测不加时 `%a` 只返回 `BL_ORD_ID:n`，`%t` 返回完整 id。

**V2-6 逐查询、逐行记账（outfmt 7）** [R1-2]。由 243 解析，逐项检查：
- (a) 每个 `# Query: <id>` 开启一个块。宣告的 id 集合与清单集合双向相等，且无重复宣告。
- (b) 每个块恰有一行 `# N hits found`，N 等于该块中的数据行数。
- (c) 每条数据行恰为 14 个 tab 分隔字段，类型正确（整数、浮点），且 `qseqid` 等于所在块的 id。
- (d) 每条数据行的 `qlen` 等于 V2-3(b) 中该查询的 `len_nt`；`sseqid` 属于 V2-5 的完整 id 集合。
- (e) 尾行 `# BLAST processed N queries` 存在，N = 63,095，尾行之后没有数据行。
- (f) 零命中块数写入记录。
- (g) **坐标合法性** [R4-2]。下文 qs、qe、ss、se、L、mm、go 分别是 `qstart qend sstart send length mismatch gapopen`；`qlen`、`slen` 是行内自报值。
  - 长度来源：行内的 `qlen` 必须等于 V2-3(b) 中该查询的 `len_nt`（即 (d)）；行内的 `slen` 必须等于基因组可信清单中 `sseqid` 这条 contig 的长度（V2-4c）。以下边界都用这两个可信值。
  - 查询坐标：`1 ≤ qs ≤ qe ≤ qlen`（BLAST 的查询坐标总是 qs ≤ qe）。
  - subject 坐标：令 `lo = min(ss, se)`、`hi = max(ss, se)`（`ss > se` 即负链，是合法的），要求 `1 ≤ lo ≤ hi ≤ slen`。
  - 跨度：`q_span = qe − qs + 1`，`s_span = hi − lo + 1`。要求 `L ≥ q_span` 且 `L ≥ s_span`。
  - 缺口与失配：缺口列总数 `gap_cols = (L − q_span) + (L − s_span)`。要求 `mm + gap_cols ≤ L`，且 `gap_cols ≥ go`。
  - 这些都是 outfmt 6 字段定义的恒等关系（`length` 为比对列数，每列要么两侧都有碱基，要么一侧为缺口），不是阈值，不需要 G-8 标注。登录节点抽样核对：遗留 shard 0 与 31 的 `hits/shard_<s>.tsv.gz` 各取前 200,000 行，共 400,000 行（约一半为负链），没有一行违反。抽样只说明当前 BLAST 版本的输出与这些关系一致，不作他用。
  - 任何一行不满足 → 失败，不发布（V2-9）。

243 必须是独立、可单测的组件。负面测试（登录节点小样本）至少包括下列情形，每一种都必须失败：
- 缺一个块；
- 多一个未知 id；
- 重复宣告；
- 尾行缺失（模拟中途截断）；
- 尾行计数不符；
- 某块删去一行数据（与宣告数不符）；
- 13 字段的行；
- 数据行出现在错误的块下；
- `qlen` 与元数据不符；
- `sseqid` 不在数据库 id 集合中；
- 坐标为 0（`qstart = 0`；另一行 `send = 0`）[R4-2]；
- 查询越界：`qlen = 1000`、`qend = 2000` [R4-2]；
- subject 越界：`max(sstart, send) = slen + 1`；负链行 `sstart = slen + 1`、`send = 1` 同样失败 [R4-2]；
- `slen` 与可信清单的 contig 长度不符 [R4-2]；
- `length` 小于 subject 跨度；`mismatch + gap_cols > length` [R4-2]。

正面测试：完整输入通过；一条合法负链行（`sstart > send`，`1 ≤ send < sstart ≤ slen`，其余关系都成立）通过，作为上面负面测试的对照 [R4-2]。另一正面依据是 V0 中 3 条查询（1 条命中、2 条不命中）的测试。

**V2-7 每一阶段都要检查。** 凡失败会影响输出的命令都检查退出码；管道启用 `pipefail`；比较遵守 G-7。

**V2-8 过滤与链接守恒** [R1-2]。
- 设 R 为 V2-6 中全部块的数据行总数，F 为比对长度 ≥ 200 的行数，D = R − F。另外独立地再数一次长度 < 200 的行，结果必须等于 D，即每条被删除的行确属长度过滤。
- 交给 239 的文件恰好是这 F 行，行数与 sha256 写入记录。
- 链接守恒：239 在 stderr 汇报的 HSP 数等于 F；copies 各行 `n_hsp` 之和等于 F；F > 0 时 copies 行数 > 0。239 会静默跳过少于 12 列的行，V2-6(c) 与这条等式一起保证没有行被跳过。
- 每条 copies 行的 `query_id` 属于清单，`assembly` 属于本 shard 的清单。
- R、F、D、Σn_hsp、copies 行数全部写入运行结果 JSON。

**V2-9 发布。** 按 G-9 发布，提交记录最后写出。提交记录列出的数据产物：`copies.tsv`、`hits.tsv.gz`（F 行过滤后命中）、`hitq.txt`、`announced.txt`；运行结果 JSON 为 `result.json`，按 G-9 单独绑定 [R3-1]。V2-2 至 V2-8 任一失败，或 G-9 的发布前字节一致检查失败：不写提交记录，运行目录的 `status` 写为 `failed:<检查名>`。

**V2-10 命名空间隔离。** V2 只写 `$W/v2/` 下的 `runs/`、`committed/`、`locks/`。不得写入、移动或删除遗留 `copies/`、`hits/`、`shardlog/`、`verify/` 中的任何文件。

**V2-11 并发** [R1-1]。
- 写入任何东西之前，先用原子 `mkdir $W/v2/locks/shard_<i>` 取锁。锁已存在则以非零状态退出，且不写任何东西。锁由持有者在退出时释放；崩溃遗留的陈旧锁须人工处理（失败即关闭）。
- 启动时若 `squeue` 中仍有遗留扫描作业在运行（作业名 `isall_C`、`isall_Cr`），直接退出。
- 已有提交记录的 shard，先按 G-9 核对整条链：通过才跳过；核对失败则以非零状态退出并报告，不得静默跳过，也不得自动覆盖。

**V2-12 运行结果 JSON（`result.json`）的必需字段。**
- G-4 的全部字段，以及 `shard`；
- blast 完整参数串、`outfmt = 7`；
- 输入摘要：`query_sha256`、`query_ids_sha256`、`queries_tsv_sha256`、基因组可信清单 sha256；
- 查询记账：`queries_expected`、`queries_announced`、`blast_processed_trailer`、`queries_with_zero_hits`；
- 基因组与数据库：`genomes`、`contigs`、`bp`、`db_sequences`、`db_residues`；
- 守恒：`R`、`F`、`D`、`sum_n_hsp`、`copies_rows`；
- 产物与脚本摘要：`copies_sha256`、`hits_gz_sha256`、`hitq_sha256`、`announced_sha256`、`script_239_sha256`、`script_243_sha256`。G-9 的字段映射为：`copies.tsv` ↔ `copies_sha256`、`hits.tsv.gz` ↔ `hits_gz_sha256`、`hitq.txt` ↔ `hitq_sha256`、`announced.txt` ↔ `announced_sha256`，每个都附对应的 `*_bytes` 字段 [R2-1]。`result.json` 不含自身的摘要或字节数字段 [R3-1]；
- 起止时间，以及 `verified`。

**V2-13 大文件不上共享盘。** outfmt 7 原始输出与未过滤的 HSP 只放在节点本地盘。

**V2-14 `VERIFIED_V2` 的定义** [用户 2]。shard i 是 `VERIFIED_V2`，当且仅当以下三条同时成立：
- `$W/v2/committed/shard_<i>.json` 存在；
- 它引用的 `result.json` 含 `"verified": true`，且 V2-12 全部字段齐全、类型正确；
- 消费方在消费时按 G-9 重新核对整条链，并且通过。

重跑、退出码为 0、存在任何 marker，都不构成 `VERIFIED_V2`。只有 V2 写 `"verified": true`。

**V2-15 命中证据** [R1-11]。
- V2 在运行目录中，从已通过 V2-6、V2-8 的 F 行过滤后命中生成 `hitq.txt`：至少有一条 ≥200 bp 命中的查询 id，排序去重。
- 宣告列表另存为 `announced.txt`。
- 二者都列入提交记录。`announced.txt` 不得用作命中证据。

---

## 项 2：241 审计（遗留 shard 与 V2 的直接对照）

审查命令：`codex-review "241 审计" review-hard`（项名沿用历史名称）

**范围**
- `scripts/246_legacy_vs_v2.py` 与 `tests/test_legacy_vs_v2.py`（待写，A-D）
- `slurm/isref_all/V1_audit.sh`（遗留证据采集与身份记录，A-D1）
- `docs/legacy_shard_inventory.md`（64 个遗留 shard 的现状清单）

`scripts/241_audit_shards.py` 与统计筛查 A-S 已移出验收要求，只作为可选的探索性分析或历史记录保留，见附录 A [用户二 4]。

**角色。** 路线 B 下，生产只消费 `VERIFIED_V2`（S-1），所以本项**不是**生产验收，也不是门控。本项只做一件事：把每个遗留 shard 与它的 V2 重跑逐行比较，报告**观察到的**新旧差异。

### A-D：遗留输出与 V2 的直接对照

**A-D1 输入与资格。**
- 遗留输入是过滤后命中（outfmt 6，14 字段，≥200 bp），按 `docs/legacy_shard_inventory.md` 中的状态取文件：
  - `COMPLETED_LEGACY` 取 `hits/shard_<s>.tsv.gz`；
  - `CHAIN_FAILED`（FAILED marker 为 "chaining failed"，说明 blast 与过滤步骤已被执行）取 `hits/shard_<s>.tsv`。
- 遗留文件的 sha256 与字节数，由 V1 在一个获批的 SLURM 作业中记录；比较时复核，不符则不比较，判 `LEGACY_UNREADABLE`。
- V2 一侧须是按 V2-14 核对通过的提交记录中的 `hits.tsv.gz`。

**A-D1a 逐行格式与身份校验（先于任何比较）** [R3-4]。246 先完整读取两侧文件并逐行校验，校验全部通过后才进入 A-D2。不得跳过任何行后继续比较。
- 遗留一侧，以下任一情形都记为异常行：
  - 文件末尾缺少换行（末行残缺），或出现空行；
  - 字段数 ≠ 14（tab 分隔）；
  - 类型不符：`length mismatch gapopen qstart qend sstart send qlen slen` 须为非负十进制整数，`pident evalue bitscore` 须为有限浮点数；`qstart ≤ qend`；
  - `length < 200`（与遗留过滤条件矛盾）；
  - 坐标不合法 [R4-2]：按 V2-6(g) 的查询坐标、subject 坐标（min/max 处理负链）、跨度、缺口与失配各条检查，但边界取**遗留行自报的** `qlen` 与 `slen`，而不是可信值；
  - `qseqid` 不属于可信查询 id 清单（`config/isref_all_trusted.tsv` 所固定的 `all_queries.fna` 的 63,095 个 id）；
  - `sseqid` 不属于由基因组可信清单推出的 shard s 完整 id 集合（`ACC|CONTIG`，同 V2-5）。
- 有任何异常行 → 判 `LEGACY_MALFORMED`（A-D3），输出按上列类别分别计数的异常行数，以及第一条异常行的行号；该 shard 不输出任何比较结果。
- `qlen` 与可信元数据中该查询的 `len_nt` 不符，**不**算异常行：例如查询序列被截断读入时就会产生这种行，它是值得报告的差异，不应被格式检查吞掉。这类行照常参与比较，并单独计数 `legacy_rows_qlen_mismatch`。`slen` 与基因组可信清单中该 contig 的长度不符，同样不算异常行（例如遗留建库时基因组被截断读入），照常比较，单独计数 `legacy_rows_slen_mismatch` [R4-2]。
- 两类 shard 的坐标检查方式因此不同 [R4-2]：遗留一侧检查行内自洽（坐标在自报的 `qlen`、`slen` 之内），自报值与可信值的差别作为观察计数；V2 一侧要求自报值等于可信值（V2-6(d)(g)），坐标以可信值为界。两类输出都含 `qlen`、`slen` 列（`C_sweep.sh:67`、`V2_verified.sh:125`），所以不需要从别处补取长度；可信值分别来自 `all_queries.tsv` 的 `len_nt` 与基因组可信清单的 contig 长度。
- V2 一侧按 V2-6(g) 的严格规则校验（V2-6 已保证这些条件，再查一次是防御性的；`qlen`、`slen` 与可信值不符在 V2 一侧也算异常）。有任何异常 → 判 `NO_V2_TWIN`，原因 `V2_ROWS_INVALID`，并输出同样的分类计数。
- idx(q) 定义为 q 在 `all_queries.fna` 中的位置（从 0 起），由可信查询清单唯一确定。通过本条校验后，两侧每一行的查询都有 idx。

**A-D2 比较方法。** 对每条查询 q，把遗留与 V2 的数据行各自作为排序后的多重集合，逐字节比较。之所以按多重集合比较，是因为 BLAST 多线程下同一查询内的行序可能不同。

**A-D3 判词**（每个 shard 恰得一个；**全部是观察，不是原因结论**）[R2-3, R3-4]。按以下顺序，取第一个满足的：
1. 资格类判词，按 A-D1 先行判定：
   - `LEGACY_OUTPUT_DELETED`：遗留输出已被删除（清单状态 `LEGACY_OUTPUT_DELETED_BY_27768584`）；
   - `LEGACY_OUTPUT_ABSENT`：状态为 `CHAIN_FAILED`，但没有过滤后命中文件；
   - `LEGACY_RUN_INCOMPLETE`：既无 done 也无 FAILED marker（遗留任务没有跑完，例如 14 小时超时）；
   - `LEGACY_UNREADABLE`：遗留文件无法完整读取、gzip 校验失败，或 sha256 与记录不符；
   - `LEGACY_MALFORMED`：遗留文件可完整读取、摘要相符，但有 A-D1a 规定的异常行 [R3-4]；
   - `NO_V2_TWIN`：没有核对通过的 V2 提交记录（原因 `NO_VERIFIED_COMMIT`），或 V2 一侧未通过 A-D1a（原因 `V2_ROWS_INVALID`）。
2. 所有查询的多重集合都相同 → `NO_DIFFERENCE_OBSERVED`（比较范围：≥200 bp 的过滤后命中）。
3. 令 b 为按 idx 升序**第一个**遗留与 V2 多重集合不同的查询位置（因此唯一确定）。若同时满足：
   - (i) 在 b 处，遗留多重集合是 V2 多重集合的**真子集**（只缺行，没有多出或改变的行；可以为空）；
   - (ii) 所有 idx > b 的查询在遗留中都没有任何行，
   则 → `PREFIX_PATTERN_OBSERVED`，报告 b、b/N 以及 b 处缺少的行数。
4. 否则 → `OTHER_DIFFERENCE_OBSERVED`，报告遗留多出与缺少的行数、涉及的查询数与 assembly 数。

**A-D4 观察与结论必须分开** [用户二 4]。
- 246 只输出 A-D3 的观察判词与差异计数，不输出任何关于"截断"或其原因的判词。源码与输出中，不得把 `truncated`、`truncation`、`截断`、`cause`、`原因` 用作判词或判定字段。
- `PREFIX_PATTERN_OBSERVED` 与多种机制都一致：查询文件短读；原始输出或过滤输出被截断写入（遗留过滤步骤没有检查退出码）；作业中断。246 不区分这些机制。`OTHER_DIFFERENCE_OBSERVED` 也可能来自数据库内容不同、BLAST 版本或数字格式不同。
- 报告中任何关于截断或其原因的陈述，都必须标为"假设"，并写明支持它的另外证据。没有另外证据时，只能陈述观察。

**A-D5 限制**（必须写入输出与报告）：
- 遗留原始输出大多已被删除，只能比较 ≥200 bp 的过滤后命中；只影响 <200 bp 命中的差异看不见。
- 遗留作业没有记录 BLAST 版本。版本或格式不同时，差异会落入 `OTHER_DIFFERENCE_OBSERVED`；这个偏差的方向是"判为有差异"，不会误判为无差异。
- 曾被怀疑短的 shard 6、9、11、14，遗留输出已被作业 27768584 删除，只能判 `LEGACY_OUTPUT_DELETED`。仍能直接对照的短 shard 是 **58 与 63**：它们的 HSP 数分别为中位数的 2.4% 与 10.8%，遗留输出仍在（见清单）。
- 清单中记录的整数边界文件大小（shard 35 为 2³¹ 字节，shard 59 为 528 MiB）只是观察，A-D 不对其作解释。

**A-D6** `NO_DIFFERENCE_OBSERVED` 也不会使遗留输出成为生产输入。生产只用 V2 输出。

**A-D7 输出。** 原子写出（临时名 + `rename`），包含：G-4 字段；全部输入（遗留文件、V2 提交记录）的 sha256；每个 shard 的判词与差异计数；A-D1a 的分类异常计数、`legacy_rows_qlen_mismatch` 与 `legacy_rows_slen_mismatch`；`NO_V2_TWIN` 的原因；A-D5 的限制说明。

**A-D8 测试**（合成数据，登录节点小样本）：
- 完全相同 → `NO_DIFFERENCE_OBSERVED`；
- 同一查询内行序打乱 → 仍为 `NO_DIFFERENCE_OBSERVED`；
- 前缀模式（b 处遗留为 V2 的真子集，之后遗留为空）→ `PREFIX_PATTERN_OBSERVED`，且 b 等于第一个差异查询；
- 遗留缺少一个基因组的全部命中，构造成差异散布在多个查询、且其后的查询仍有遗留行 → `OTHER_DIFFERENCE_OBSERVED`；
- 某一行数字格式不同 → `OTHER_DIFFERENCE_OBSERVED`；只有最后一个有命中的查询发生格式变化 → 同样是 `OTHER_DIFFERENCE_OBSERVED`（b 处遗留不是 V2 的子集）；
- 遗留 sha256 与记录不符 → `LEGACY_UNREADABLE`；
- 无 marker → `LEGACY_RUN_INCOMPLETE`；`CHAIN_FAILED` 但没有文件 → `LEGACY_OUTPUT_ABSENT`；
- 遗留中出现一个不在可信清单中的查询 id（其余各行与 V2 完全相同）→ `LEGACY_MALFORMED`，未知查询计数为 1，**不得**输出 `NO_DIFFERENCE_OBSERVED` [R3-4]；
- 遗留末行只有 9 个字段且没有结尾换行（sha256 与记录相符）→ `LEGACY_MALFORMED`，末行残缺计数为 1 [R3-4]；
- 遗留某行 `qlen` 与元数据不符 → 不判畸形，照常比较（结果为差异类判词），`legacy_rows_qlen_mismatch` = 1；
- 遗留某行 `qend` 大于该行自报的 `qlen`；另一行负链 `sstart` 大于自报的 `slen` → `LEGACY_MALFORMED`，坐标类异常计数为 2；一条合法负链行作为对照，不计异常 [R4-2]；
- 遗留某行 `slen` 与可信 contig 长度不符、但坐标在自报 `slen` 之内 → 不判畸形，`legacy_rows_slen_mismatch` = 1 [R4-2]；
- V2 一侧出现异常行 → `NO_V2_TWIN`，原因 `V2_ROWS_INVALID`；
- 输出中不出现 A-D4 禁用的词。

---

## 项 3：242 位点合并

审查命令：`codex-review "242 位点合并" review-hard`

**范围**
- `scripts/242_collapse_sites.py`
- `scripts/245_query_orientation.py` 与 `tests/test_query_orientation.py`（待写，G-3 的 `query_orient` 表）
- `slurm/isref_all/S1_sites.sh`（待写，驱动作业）
- `tests/test_collapse_sites.py`（待写）

**目的**：把 copy record 归并成位点，即跨查询去重。输出是**候选**插入拷贝位点：既未经验证（验证属 Stage E），也不是 events。

### 验收标准

**S-1 输入门控与部分结果** [用户 2]。
- 生产模式只消费 `VERIFIED_V2` shard，每次消费前按 G-9 核对。
- 每条输出记录都带 `shard_id`。
- 每个 shard 的位点结果按 G-9 发布为单元 `sites_shard_<i>`。
- 汇总报告只有在全部 64 个 shard 都是 `VERIFIED_V2`，且 64 个 `sites_shard_<i>` 都已提交时，才可称为最终结果。否则标为 `PARTIAL`，列出纳入的 shard，并逐个写明未纳入的 shard 及原因 [R2-7]。原因只能取：`NOT_RUN`（尚无 V2 运行）、`NO_COMMIT`（有运行目录但没有提交记录，包括运行失败）、`CHAIN_CHECK_FAILED`（V2 提交记录未通过 G-9 核对）、`SITES_NOT_COMMITTED`（位点单元尚未提交）、`SITES_CHAIN_CHECK_FAILED`（位点单元未通过 G-9 核对）。
- 遗留 shard 不进入任何位点计数。

**S-2 "shard 间基因组不相交"前提必须被检查**，否则致命：
- (a) 64 份清单两两不相交，并集等于 40,541（一次性检查）；
- (b) 每条输入行的 assembly 属于该 shard 的清单；
- (c) 已消费的 shard 之间没有共享的 assembly。

只有这个前提成立，"按 shard 独立归并即全局精确"才成立。

**S-3 分组中间文件与输入格式** [R1-12]。
- 239 按 (query, subject) 首次出现的顺序输出，同一 (assembly, contig) 的行并不连续。驱动作业须从核对过的 copies 文件生成按 (assembly, contig) 分组的中间文件，保留表头。
- 须核对数据行数守恒，并把中间文件与来源 copies 文件的 sha256 都写入元数据。
- 242 读取中间文件时，以下情形致命：
  - 首行不逐字等于 239 的 16 列表头：`query_id bag_id assembly species contig s0 s1 strand n_hsp qcov pident copy_kind catalogue_status catalogue_occurrence catalogue_bag same_bag`；
  - 任何一行字段数 ≠ 16、s0/s1 不是整数，或 s1 − s0 ≤ 0（机器写出的文件出现畸形行，意味着截断或损坏）；
  - 任何一行 `s0 < 0`，或 `s1` 大于基因组可信清单中 (assembly, contig) 这条 contig 的长度 [R4-2]。按 G-2，合法的 V2 输入总满足 `0 ≤ s0 < s1 ≤ contig 长度`（V2-6(g) 保证每条 HSP 的 subject 坐标在 `[1, slen]` 内）；这里再查一次，防止中间文件出错。
  - 某个 (assembly, contig) 键在其分组处理完之后再次出现。

**S-4 记录分类与"同一元件"关系。**
- **锚记录**：`copy_kind = full_copy`（qcov ≥ 0.90，继承自 239），且长度满足 0.8 × qlen ≤ s1 − s0 ≤ 1.2 × qlen。0.8 和 1.2 为未校准操作值；qlen 取自核对过的 `all_queries.tsv`。
- **超长记录**：`full_copy`，但长度超出上述范围。例如 239 把相邻串联拷贝链成了一条，或拷贝内部被另一个插入打断。超长记录不作锚，见 S-7。
- **非锚记录**：`partial_copy` 与 `fragment`。
- **同一元件 same(q1, q2)**：两条查询属于同一 bag；或者二者的代表序列在 237 的 `query_equivalence.tsv` 中被直接列为一对。代表序列取自 236 的 `query_decisions.tsv`：被吸收者取其 `absorbed_by`，其余取自身。不取传递闭包。

**S-5 锚之间的关系（核心规则）** [R1-13, 用户 3]。设两条锚 a、b，len(a) ≥ len(b)，ov 为重叠长度，c_b = ov / len(b)，c_a = ov / len(a)（因此 c_a ≤ c_b）：

| 条件 | 判定 |
|---|---|
| ov ≤ 0 | 不相交 |
| same(a, b) 且 c_b ≥ 0.80 | **同一拷贝**：b 作为 a 所在拷贝的支持，不另成位点。包括同一元件的较短变体完全落在较长变体之内的情形 |
| same(a, b) 且 0 < c_b < 0.80 | `OVERLAP_UNRESOLVED` |
| 非 same，c_a ≥ 0.80 且 c_b ≥ 0.80 | **同一物理拷贝，元件身份不定**（`ELEMENT_AMBIGUOUS`）：计一个位点，列出候选 bag |
| 非 same，c_b ≥ 0.90 且 c_a < 0.80 | **嵌套候选**（`NESTED_CANDIDATE`）：只说明 b 落在 a 内，本身不足以证明 b 是另一个插入；是否单独计数见 S-5a |
| 非 same，其余重叠情形 | `OVERLAP_UNRESOLVED` |

- 阈值 0.80、0.90 为未校准操作值。
- 依据：物理上不同的拷贝占据互不重叠的碱基，所以两条全长锚高度重叠，说明它们描述的是同一段碱基。成对重叠只作为必要条件，**不作为生物学证据**。
- 不再使用上一版的 0.50：一个 1 kb 锚落在 2 kb 锚内时，c_a 恰为 0.50，因此 0.50 无法区分"同一拷贝"与"嵌套"。

**S-5a 嵌套候选的升级条件** [用户二 5, R2-5, R2-6, R3-2, R3-3]。"不同 bag（非 same）+ 包含度 ≥ 0.90"只构成嵌套**候选**。候选必须**同时**满足 E1 与 E2，才升级为单独计数的嵌套位点 `NESTED`。
- **E1 元件不同**：把内层锚的查询 q_b 与外层锚的查询 q_a 比对（G-5 参数），q_b 被 q_a 覆盖的合并长度 < 0.50 × len(q_b)（未校准操作值）。覆盖 ≥ 0.50 **只表示未满足 E1**：b 保持候选，既不升级，也不因此被吸收。
- **组成部分（COMPONENT）须有正面证据**：b 只有同时满足下面两条，才判为外层元件的组成部分，作为 a 所在单元的支持、不另成位点：
  - (i) q_b 被 q_a 覆盖 ≥ 0.90 × len(q_b)（未校准操作值）；
  - (ii) b 的碱基本身比对在外层元件上，且中间没有插入 [R3-3]：在外层记录的 HSP 链（按下面 E2 第 1–2 步重建）中，存在**单条** HSP h，满足：
    - h 的 subject 区间覆盖 b 区间 ≥ 0.90 × len(b)（未校准操作值）；
    - h 内两个方向的缺口列总数**分别** ≤ 50 bp（未校准操作值）[R4-1]：`L − q_span ≤ 50` 且 `L − s_span ≤ 50`，记号同 V2-6(g)。`L − q_span` 是 query 一侧为缺口的列数，即 subject 相对 query 多出的碱基总数；`L − s_span` 是另一方向。两者分别约束，插入与删除不能互相抵消。

    依据：outfmt 6 的 `length` 是比对列数，所以缺口列总数可以由端点与 `length` 精确算出，不需要缺口位置，可以机械检查（V2-6(g) 已保证两者非负）。同一条 HSP 连续覆盖 b 的 90%，且 subject 一侧多出的碱基总数 ≤ 50 bp，说明 h 覆盖的这段碱基与外层查询连续共线，其中没有超过 50 bp 的额外插入。它不声称 b 精确对应 q_a 的哪一段。
    仍然存在的误吸收情形（偏少方向，写入 S-8）：b 中未被 h 覆盖的不足 10% 部分不受约束；外层查询 q_a 自身含有与 q_b 同源的片段、且真实插入恰好与之等长对齐时，本条无法区分。

  缺任一条，或链重建不可评估（E2 第 1 步），b 保持候选。
- **E2 插入边界**（绑定到具体的外层记录；50 bp 均为未校准操作值）：
  1. **链重建**：从 V2 已核对的 `hits.tsv.gz` 中取出 (q_a, 该 contig) 的全部 HSP，用与 239 `chain()` 逐字节一致的实现重建链。取方向等于外层记录 `strand`、包络等于外层记录 `[s0, s1)`（G-2）、且 HSP 数等于外层记录 `n_hsp` 的链。找不到恰好一条这样的链 → E2 不可评估。
  2. **坐标**：每条 HSP 的 subject 区间为 `[min(sstart, send) − 1, max(sstart, send))`，query 区间为 `[qstart − 1, qend)`（BLAST 输出中 qstart ≤ qend）。
  3. **分块**：链中 subject 终点 ≤ b.s0 + 50 的 HSP 为左块，subject 起点 ≥ b.s1 − 50 的为右块。两块都必须非空，且每条 HSP 恰好属于一块：链中任何 HSP 都不得覆盖 (b.s0 + 50, b.s1 − 50) 内的碱基。
  4. **边界 HSP**（必须唯一）[R3-2]：左边界 HSP h_L 是左块中 subject 终点最大的那条，右边界 HSP h_R 是右块中 subject 起点最小的那条。最大值或最小值出现并列 → E2 不可评估。记 L1 = h_L 的 subject 终点，R0 = h_R 的 subject 起点。要求 |L1 − b.s0| ≤ 50，且 |R0 − b.s1| ≤ 50。
  5. **整条链的查询单调性** [R3-2]：把链中全部 HSP 按（subject 起点，subject 终点）升序排列；两条 HSP 的 subject 区间完全相同 → E2 不可评估。对每一对相邻 HSP（前者 h、后者 h′）：
     - `+` 链：h′ 的 query 起点 − h 的 query 终点 ≥ −50；
     - `−` 链（subject 增大时 query 减小）：h 的 query 起点 − h′ 的 query 终点 ≥ −50。

     即沿 subject 方向，query 坐标单调推进，相邻 HSP 的 query 重叠不超过 50 bp。任何一对不满足 → E2 不成立。
  6. **边界处共线** [R3-2]：只用两条边界 HSP 核对 b 两侧的查询间隙。`+` 链要求 −50 ≤ (h_R 的 query 起点) − (h_L 的 query 终点) ≤ 50；`−` 链要求 −50 ≤ (h_L 的 query 起点) − (h_R 的 query 终点) ≤ 50。不使用块内 query 区间的并集。
- E1 或 E2 所需的数据缺失、或不可评估时，候选保持候选，不升级。
- 升级时：b 为 `NESTED`，记录 `nested_in`；a 所在单元标 `interrupted_by`。
- 未升级的候选：不另计位点；外层单元标 `has_nested_candidate`；候选数单独报告。
- 数据来源：E1 与组成部分判定的 (i) 来自新组件 `scripts/247_query_pair_containment.py`（只对候选对做查询对查询的比对，记录覆盖度与 q_a 坐标区间）；E2 与 (ii) 来自 V2 提交记录中的 `hits.tsv.gz`。

**S-6 拷贝单元的构建与串联桥接** [R1-13, R2-4]。
- **处理顺序固定**：
  1. 按 S-4 给记录分类；
  2. 按 S-7 解析超长记录：满足条件者成为"被打断锚"，长度检查改用 len(q_a) + len(b)；用于解析的 b 直接记为 `NESTED`；
  3. 用锚与被打断锚构建拷贝单元（下列规则）。被打断锚之间、被打断锚与锚之间都按 S-5 判定，所以多条查询各自解析出的同一外层拷贝，只要 same 且 c_b ≥ 0.80，就归入同一单元，只计一次；
  4. 按 S-5a 评估其余嵌套候选；
  5. 按 S-7 分配非锚记录。
- 锚按（长度降序、s0 升序、s1 升序、query_id 升序）排序后依次处理。一条锚加入**第一个**被 S-5 判为"同一拷贝"或 `ELEMENT_AMBIGUOUS` 的种子单元；没有可加入的单元时自成种子。
- 成员只与种子比较：不做传递合并，种子区间也不因吸收成员而扩张。嵌套候选与嵌套位点的父单元，取包含它的最短单元。
- 单元内若有两条**互不重叠**的全长锚（例如种子 [0,2000) 内的 [0,1000) 与 [1000,2000)），该单元标为 `AMBIGUOUS_STRUCTURE`（可能是串联或复合结构）：既不静默计为一个确定位点，也不拆成两个位点。
- 单元内同一 bag、`element_orient` 均已知且相反的成员 → 该单元标 `ORIENTATION_INCONSISTENT`，计入不确定结构。`element_orient` 未知时不判冲突（G-3）。
- 输出分别给出种子坐标（`seed_s0, seed_s1`）和成员包络（`env_s0, env_s1`）。

**S-7 非锚记录与超长记录** [用户 3]。
- 非锚记录与某个拷贝单元种子区间的包含度 ov / len(记录) ≥ 0.80（未校准操作值）时，作为该单元的支持，**不另成位点**。落在多个单元内时，优先取 same 元件的最内层单元，没有就取最内层单元；记录二者是否 same。
- 0 < 包含度 < 0.80 时，记作该单元的边界记录，同样不另成位点。
- 与任何单元都不重叠的非锚记录，按 S-6 的种子规则（以最长者为种子）归成"仅部分支持"位点组 `PARTIAL_ONLY`，单独报告。
- 每条超长记录形成一个 `AMBIGUOUS_STRUCTURE` 实体，落在其内的拷贝单元保留，并标 `inside_ambiguous_structure`。唯一的例外：超长记录内含一条非 same 的锚 b，b 与该记录满足 S-5a 的 E1 与 E2（E2 的链即该超长记录自身的链），且记录长度与 len(q_a) + len(b) 相差 ≤ 100 bp（未校准操作值）。此时该记录解析为**被打断锚**，在 S-6 第 3 步与其他锚一起构建单元（同一外层拷贝的多条被打断锚因此只计一次）；所在单元计为 `COPY`（标 `interrupted_by`），b 为 `NESTED` [R2-4]。

**S-8 位点等级与计数口径** [用户 3, 用户二 5]。
- `COPY`：拷贝单元，元件唯一，没有不确定标记（包括被 `NESTED` 打断、标 `interrupted_by` 的外层拷贝）。
- `COPY_ELEMENT_AMBIGUOUS`：同一物理拷贝，元件身份不定。
- `NESTED`：满足 S-5a（E1 与 E2）的嵌套位点。单独计数，同时计入标题。
- `NESTED_CANDIDATE`：未满足 S-5a 的嵌套候选。单独报告，**不计入标题**。
- `AMBIGUOUS_STRUCTURE`：S-6 的串联或方向不一致单元、S-7 中未解析的超长记录实体，以及涉及 `OVERLAP_UNRESOLVED` 的单元。
- `PARTIAL_ONLY`：只有部分或片段支持的位点组。

**标题计数 = COPY + COPY_ELEMENT_AMBIGUOUS + NESTED**。`NESTED_CANDIDATE`、`AMBIGUOUS_STRUCTURE`、`PARTIAL_ONLY` 永远单独报告，不并入标题计数。

须同时声明两个方向的偏差：
- **偏少**：
  - 真实拷贝被归入不确定结构（串联拷贝被链成超长记录时，只得到一个 `AMBIGUOUS_STRUCTURE`）；
  - 真实的嵌套插入缺少 E1 或 E2 证据时只作候选，不计入（例如 E2 需要的 HSP 边界因序列分化而不清楚）；
  - 同一拷贝的两条同元件锚因长度变异 c_b < 0.80，被判 `OVERLAP_UNRESOLVED`；
  - 部分拷贝紧贴全长拷贝时，被并为边界记录；
  - 真实的内层插入在 S-5a COMPONENT (ii) 的剩余盲区内（见该条）被误判为组成部分 [R4-1]。
- **偏多**：
  - 239 的链接间隙规则把一个拷贝拆成两条锚时，会算成两个位点；
  - 两个相同元件的独立插入恰好首尾相接、又被判为不同元件（无 237 等价对）时，S-5a 可能误把其一升级为嵌套。E1、E2 同时成立使这种情形少见，但无法排除。

**S-9 方向信息保留** [用户 4, R3-5]。每个位点输出：
- 成员 `blast_strand` 的 + / − 计数；
- **按 bag 分别**给出 `element_orient` 共识 [R3-5]：对位点内每个 bag，取该 bag 成员中已知的 `element_orient`：全为 `+` → `+`；全为 `−` → `−`；两者都有 → `MIXED`；没有已知值 → `UNKNOWN`。`MIXED` 只在同一 bag 内部使用；
- 位点方向 `site_orient`（"已知共识"指 `+`、`−` 或 `MIXED`）：
  - 恰有一个 bag 给出已知共识、其余 bag 都为 `UNKNOWN` → 取该 bag 的共识（可以是 `MIXED`，此时它来自同一 bag 内部），并记录所依据的 bag；
  - 两个或更多 bag 给出已知共识 → 因为不同 bag 的参考序列之间没有共同方向基准（G-3），取 `UNKNOWN`，不得汇总为 `MIXED`；
  - 没有任何 bag 给出已知共识 → `UNKNOWN`；
- 每条锚的 (query_id, bag_id, blast_strand, query_orient, element_orient)，另存为明细文件 `orient_detail.tsv`。

合并判定不使用 `blast_strand`。

**S-10 确定性。**
- 种子顺序如 S-6；代表锚按（qcov 降序、pident 降序、长度降序、query_id 升序）选取。
- 组内输入行序打乱后，输出逐字节相同。

**S-11 复杂度。**
- 候选只在与种子重叠的记录中搜索（区间索引，或排序后二分查找），不得对每个 contig 做全两两比较。
- 验收：行数最多的已消费 shard 在单个 cf1 节点上 4 小时内完成（未校准的预算）。
- 背景：遗留 shard 平均约 4,900 万条 copy record，约每个基因组 7.9 万条。

**S-12 计数恒等。** 每条输入数据行恰好归入一个去向：锚成员（含被打断锚、嵌套候选与嵌套位点单元）、支持（含有正面证据判为组成部分的内层锚）、边界记录、`PARTIAL_ONLY` 成员、未解析的超长实体。各去向之和等于输入行数，否则致命。该等式写入元数据。

**S-13 与旧 loci 数字的历史对照** [R1-15]。在同一输入上，另报按旧规则 `isref_v1/jobs/loci.sh` 得到的计数。旧规则：只取 full_copy；按 s0 排序；`s0 < ce − 50` 时合并，且 `ce = max(ce, s1)`；传递。同时须：
- 另报只取旧 181 条查询序列的子集计数。按序列 sha256 与 `is_like_queries.fna` 匹配，找不到的逐条列出。
- 写明消费的 shard 范围（全部 64 个才与旧扫描的基因组范围一致）、查询集合，以及其余差异：旧扫描的 blast 参数与本扫描相同（`sweep.sh:46-47` 对 `C_sweep.sh`）；旧扫描本身也受过 scratch 读故障影响。
- 该计数定位为历史算法对照，不用于校准位点数，也不作验收门槛。

**S-14 判词与目录。**
- `catalogue_status` 只取 `catalogued` 或 `not_in_catalogue`。位点只要有至少一个成员记录被 239 判为 catalogued，就是 `catalogued`。
- 242 的输出不得出现 "validated" 或 "new site"。

**S-15 已知限制必须声明** [R1-14]。
- 环状 contig 上跨原点的拷贝：239 不携带 contig 长度，本项不处理。两端的片段各自可能因 <200 bp 而被过滤，结果可能是零个、一个或两个记录，偏差方向不固定。
- 元数据的 `limitations` 字段须列出这一条，以及 S-8 的偏差说明。

**S-16 发布与元数据。** 按 G-9 发布：数据产物为位点文件 `sites.tsv` 与方向明细 `orient_detail.tsv`；运行结果 JSON 为元数据 `meta.json`，按 G-9 单独绑定 [R3-1]。元数据包含：
- G-4 字段；
- 输入（提交记录、copies、中间文件）的路径、字节数、sha256；
- 全部阈值及其标注；
- 各等级的计数、S-12 的恒等式、S-13 的对照；
- `limitations`；
- 未消费的 shard 及其原因（取值见 S-1）[R2-7]。

`meta.json` 在 G-9 中充当运行结果 JSON：它必须声明两个数据产物的摘要与字节数（`sites.tsv` ↔ `sites_sha256`、`orient_detail.tsv` ↔ `orient_detail_sha256`，各附 `*_bytes`），数据产物按 G-9 三方核对 [R2-1]；`meta.json` 不含自身的摘要，按 G-9 与提交记录两方核对 [R3-1]。

**S-17 测试**（小样本 fixture，登录节点）：
- (a) 1.3 kb 的非同元件锚完全落在 5 kb 锚内，二者长度都正常，没有 E1/E2 数据 → 一个 `COPY` 加一个 `NESTED_CANDIDATE`，标题计数为 1；
- (a2) 同上，但 E1 成立（q_b 不被 q_a 覆盖），且 5 kb 元件的 HSP 恰在 1.3 kb 锚两侧被隔开（E2 成立；外层为长度 ≈ len(q_a) + 1.3 kb 的超长记录）→ 外层 `COPY`（`interrupted_by`）加 `NESTED`，标题计数为 2；
- (a3) q_b 被 q_a 覆盖 ≥ 0.90，且外层链中单条 HSP h 覆盖 b ≥ 0.90 × len(b)、`L − q_span ≤ 50` 且 `L − s_span ≤ 50`（组成部分的正面证据，与 S-5a COMPONENT (ii) 一致）→ 内层作为支持，标题计数为 1 [R5-1]；
- (a3c) 插入与删除互相抵消 [R4-1]：q_b 被 q_a 覆盖 ≥ 0.90，单条 HSP 覆盖 b，subject 跨度与 query 跨度相等，但 `length` 比两者都大 1,300（h 内一处 1,300 bp 的 subject 插入与一处 1,300 bp 的删除）→ `L − q_span = 1300 > 50`，(ii) 不成立，b 保持 `NESTED_CANDIDATE`；
- (a3b) q_b 被 q_a 覆盖 ≥ 0.90，但外层链在 b 处断开（没有单条 HSP 覆盖 b 的 90%），或覆盖 b 的那条 HSP subject 长度比 query 长度多 1.3 kb → b 保持 `NESTED_CANDIDATE`，标题计数为 1 [R3-3]；
- (a4) q_b 与 q_a 共享超过一半序列（覆盖 0.60），并且有 E2 边界证据 → E1 未满足，b 保持 `NESTED_CANDIDATE`，既不升级也不吸收，标题计数为 1；
- (a5) 两条同 bag 的 5 kb 查询各自对同一个被 1.3 kb 插入打断的外层拷贝产生约 6.3 kb 的超长记录，二者都可解析 → 一个外层 `COPY` 加一个 `NESTED`，标题计数为 2，不是 3；
- (a6) E2：合法的负链情形 → 升级；左右块来自相反方向、或来自不是该外层记录自身的链（伪边界）→ 不升级；
- (a7) E2 伪边界反例（第 3 轮审查给出）[R3-2]：全部为 `+` 链，三条 HSP 的 (subject, query) 半开区间为 `([0,1500),[1000,2500))`、`([1500,2500),[0,1000))`、`([3800,6300),[2500,5000))`，b = `[2500,3800)` → 第 5 步失败（第二条相对第一条 query 倒退 2,500 bp），第 6 步也失败（h_L 的 query 终点 1000，h_R 的 query 起点 2500，间隙 1,500）→ 不升级；作为超长记录也不得经 S-7 解析为被打断锚；
- (a8) E2 边界 HSP 并列：左块两条 HSP 的 subject 终点相同 → 不可评估，不升级 [R3-2]；
- (b) 同一 bag 的两条近同锚 → 一个位点；
- (c) 同一元件的 1.3 kb 变体落在 2.5 kb 变体内 → 一个位点（作为支持），不判嵌套；
- (d) 种子 [0,2000) 内有两条互不重叠的同元件锚 [0,1000)、[1000,2000) → 一个 `AMBIGUOUS_STRUCTURE` 单元，标题计数为 0；
- (e) 超长记录（长度 2 × qlen）→ `AMBIGUOUS_STRUCTURE`，不作种子；
- (f) 片段落在锚内 → 作为支持，位点数不变；
- (g) 方向：两条同 bag 锚 `blast_strand` 相反、`query_orient` 也相反 → 一个 `COPY`，不标冲突；`query_orient` 为 `NA` → 不标冲突；`element_orient` 均已知且相反 → `ORIENTATION_INCONSISTENT`；
- (h) 非同元件锚互相覆盖 ≥ 0.80 → 一个 `COPY_ELEMENT_AMBIGUOUS`；
- (h2) 跨 bag 方向 [R3-5]：(h) 的两条锚分属 bag X、Y，X 的共识为 `+`，Y 的共识为 `−` → 各 bag 共识照实输出，`site_orient = UNKNOWN`，不标 `MIXED`，不标 `ORIENTATION_INCONSISTENT`；只有 X 的共识已知 → `site_orient` 取 X 的共识并记录依据 bag X；
- (i) 打乱行序 → 输出相同；
- (j) 缺表头、键重现、畸形行、assembly 不在清单、`s0 = −1`、`s1` = contig 长度 + 1 → 致命；`s1` 恰等于 contig 长度 → 通过 [R4-2]；
- (k) 计数恒等成立；
- (l) 跨原点边界 fixture（一端 <200 bp）→ 按 S-15 记为限制，不报错。

---

## 项 M1-B：Arm B 结构等位发现（审查与修复）[用户 2026-11-11 第十次决定]

审查命令：`CODEX_REVIEW_MODEL=gpt-6-astra CODEX_REVIEW_TIMEOUT=2400 codex-review "M1 Arm B" review-astra`。模型与推理强度由用户指定：gpt-6-astra，high。本项属于 M1 模块，与项 0–3（M5 全语料扫描）相互独立；审查者只按本节审，第 0 节"通用约定"中的 G-4（记录提交）与"登录节点只跑小样本"适用，其余 M5 条款不适用。

**范围**
- 主审：`scripts/30_armB_graph.py`、`scripts/lib/tsd.py`、`scripts/lib/util.py`、`scripts/lib/fasta.py`。
- Arm B → Layer 2 的交接：`tools/discovery_pipeline.py`、`tools/discovery_evidence.py`。
- 测试：`tests/test_armB_calls.py`、`tests/test_discovery_pipeline.py`，以及本项新增的测试。
- 只读参照：`scripts/70_allele_reconstructor.py`（Layer 2，2026-08-28 锁定），只用于核对交接契约。对 70 本身的意见记入 `docs/DECISIONS.md`，本项不修改它；要改须用户另行批准。
- 文档：`docs/DISCOVERY_PIPELINE.md`。
- 只读参照（被上述文件引用，提供全文以免缺材料；不审，意见可列出但标 `[范围外]`）：`scripts/lib/paf.py`、`scripts/lib_insert.py`、`tools/audit_armB.py`、`tools/audit_armB_layer2.py`、`tools/compare_pangraph_mummer.py`、`tools/verify_refactored_local.py`、`tools/local_debug.py`、`tools/local_armB_layer2.py`（本地 R6/TIGR4 个案工具）。

**不在范围**
- 过去的运行结果（catalogue / database v2、v3 等）不复核、不重跑。本项的改动只影响今后的输出；v3 所用代码的版本记录在 `docs/CATALOGUE_V3_PROVENANCE.md`。
- Arm A（20）、stage 32、50、60–64、86、87、96–99。

**目的**：Arm B 是主路线的第一层，它的事件表、序列和接合坐标直接成为 Layer 2 的种子。这里的错误会以格式正确的行流到下游，因此每一条都要求"出错即可见"。

**验收标准**
- **B-1 失败不得伪装成结果。**
  - 任一外部工具失败、输入不合格或解析失败时，退出码非零，并给出 FATAL 原因。
  - 最终文件名下的输出（`armB_events.tsv`、`armB_large_events.tsv`、各 FASTA 等）只在全部写完后才出现，可用临时文件加重命名，或完成标记。中途崩溃不得留下可被下游当作已完成的非空最终文件。
  - `--reuse` 只能复用与当前图一致的产物：图重新构建后，不得复用旧的 call 文件。
- **B-2 解析不猜测。**
  - `gfatools bubble` 与 `minigraph --call` 的每一行，要么按文档格式完整解析，要么计为无法解析，并在日志和汇总中报告数目。不得用"取第一个整数"之类的回退去猜等位长度。
  - 缺失调用（`*`、`.`）与长度为 0 的调用必须区分。
- **B-3 等位分类不静默丢弃。**
  - 每个泡上的每个基因组只能归入以下之一：空、满、中间（长度落在空与满的阈值之间）、缺失。中间与缺失的数目写入事件表。
  - 分类阈值（`--len-tol`、`--min-insert`）只在一处定义。
- **B-4 尺寸门与停放。**
  - `insert_size` 在 [min_insert, max_insert] 内的进入主表；大于 max_insert 的停放到 `armB_large_events.tsv`，并保存序列。
  - 没有泡因为尺寸被静默丢弃；小于 min_insert 的不是事件，但要计数报告。
  - `parent_large_event_id` 列要么正确填写，要么从表中移除并在文档中说明，不得长期写入无意义的占位值。
- **B-5 序列与调用一致。**
  - `armB_inserts.fna` 中的插入序列必须来自被判为空与满的那两个等位，或经核对与之一致。
  - 导出的插入长度 `insert_len_direct` 与按调用得到的 `insert_size` 相差超过 `--len-tol` 时，必须在 flags 中标出。
  - 有两种以上不同长度的泡（多等位）必须能被识别。
- **B-6 接合坐标的坐标系。**
  - `junction_L` / `junction_R` 以空等位参考基因组（`empty_ref_genome`）的 contig 坐标给出。坐标约定（0 起始、半开）在代码和文档中写明，并有测试；正链与负链都必须正确。
  - 计算所用的锚长度必须来自与该参考基因组自身序列一致的序列。若该基因组在此泡处的序列与用来求锚的序列（或其反向互补）不一致，不得给出坐标：标为无坐标系，或另行标记，而不是套用不对应的偏移。
- **B-7 TSD 只作注释。** 任何代码路径都不得按 `tsd_*` 过滤事件，或改变事件是否输出，包括交接给 Layer 2 的部分。
- **B-8 交接给 Layer 2。** 每个主表事件都出现在候选表中：可播种的给出种子（区间为接合 ± 200 bp，且落在 contig 内）；不可播种的保留并写明原因。缺少基因组或 contig 时给出原因，而不是抛出异常。
- **B-9 零泡默认失败。** 退出码 3；用 `--allow-zero-bubbles` 可改为记录。
- **B-10 测试。** 在登录节点上用小样本运行，可用桩程序代替 minigraph 和 gfatools。至少覆盖：
  - B-1：工具失败、中途崩溃、reuse 与图不一致；
  - B-2：无法解析的行；
  - B-3：中间等位；
  - B-4：>5 kb 停放、<min_insert 计数；
  - B-5：导出长度与调用不一致被标出；
  - B-6：正链、负链，以及参考基因组序列与锚序列不一致；
  - B-8：缺 contig。

**轮数**：同一项最多 3 轮（CLAUDE.md）。对 Layer 2（70）的意见只记录，不在本项修改。

---

## 附录 A：A-S 截断迹象筛查（可选的探索性分析；不属于任何验收要求）[用户二 4]

本附录保留第 1 轮之后修订过的统计筛查设计，供需要时做探索性分析，或作为历史记录。它**不是**任何项的验收标准：审查项 2 时不审查本附录，`scripts/241_audit_shards.py` 也不是项 2 的必需材料。若将来要使用，须另立 SPEC 项审查。A-S 的任何判词都不得作为 shard 的验收、重做或排除依据。


**定义**
- N = 63,095。idx(q) 为 q 在 `all_queries.fna` 中的位置（从 0 起）。
- H(s)：遗留 shard s 的过滤后命中中，至少出现一行的查询集合。来自 V1 证据 `hitq_<s>.txt`。
- 参考集：只取 `VERIFIED_V2` shard 的 `hitq.txt`，经 G-9 核对链获得 [R1-7, R1-11]。被审计的遗留数据不得用于构建 R。
- R：在 ≥ 90%（未校准操作值）的参考 shard 中出现的查询集合；p_r 为查询 r 的出现率。
- M(s) = max{ idx(q) : q ∈ H(s) }。
- T_bag(s)：R 中 idx > M(s) 的成员所属的**不同 bag** 数 [R1-6]。
- B(s)：R 中 idx ≤ M(s) 的成员里，不在 H(s) 中的比例。

**A-S1 判词表封闭。**
- 判词只能是 `SUSPECT_TRUNCATION`、`SUSPECT_OTHER`、`NO_SIGNAL`、`CORRUPT`、`NOT_AUDITED`、`UNDECIDABLE` 之一。
- `NO_SIGNAL` 的意思是"没有发现迹象"，不是"完整"。
- 241 的源码与输出中，不得把 `complete`、`COMPLETE`、`OK`、`verified`、`IDENTICAL` 用作 A-S 的结论。

**A-S2 证据采集与身份绑定** [R1-1, R1-10]。
- V1 对每个遗留 shard 记录：命中 gz 与 copies 文件的字节数、mtime、**sha256**，以及采集状态 `READ_OK` 或 `READ_FAILED:<原因>`。
- 证据先写临时名，再 `rename`。
- 采集失败时：原子写出状态记录 `READ_FAILED`，并删除该 shard 旧的 `hitq_<s>.txt`，不得沿用旧证据。

**A-S3 判定优先级**（每个 shard 恰得一个判词；按以下顺序取第一个满足的）[R1-9, R1-10]：
1. 采集状态为 `READ_FAILED`，或 gzip 校验失败 → `CORRUPT`；
2. 没有遗留 done marker、没有证据，或证据记录的 sha256 与当前文件不符 → `NOT_AUDITED`；
3. 全局前提（A-S5）不满足 → `UNDECIDABLE`；
4. H(s) 中出现清单之外的查询 id → `CORRUPT`（证据与查询集不一致）；
5. H(s) 为空 → `SUSPECT_TRUNCATION`，原因 `EMPTY`；
6. T_bag(s) ≥ 5（未校准操作值）→ `SUSPECT_TRUNCATION`；
7. B(s) 的分母为 0（R 中没有 idx ≤ M(s) 的成员）→ `UNDECIDABLE`；
8. B(s) > median_ref(B) + k × 1.4826 × max(MAD_ref(B), 0.001)，k = 5 → `SUSPECT_OTHER`。k、常数 1.4826 和下限 0.001 都是未校准操作值；
9. 否则 → `NO_SIGNAL`。

median_ref 与 MAD_ref 在参考 shard 上用留一法计算：每个参考 shard 的 B，用排除它自身后构建的 R 计算 [R1-8]。

**A-S4 不作统计保证** [R1-6, R1-8, 用户 1]。
- T_bag 与 B 的阈值都是未校准操作值。SPEC 与输出不得给出误报率、族错误率或任何概率上界。
- 独立性不成立：同一 bag 的近同查询、同一质粒上的元件会一起缺席。若输出 ∏(1 − p_r) 一类的量，只能标为"假设独立时的描述量"。
- 须同时输出 k 取 3、5、7 时各 shard 的判词（敏感性），判词以 k = 5 为准。

**A-S5 全局前提**（任一不满足，所有可审计 shard 都判 `UNDECIDABLE`，并写出带原因的判词文件）[R1-7, R1-9]：
- (a) 参考 shard 数 ≥ 8（未校准操作值）；
- (b) R 所含的不同 bag 数 ≥ 50（未校准操作值）；
- (c) 把 idx 按十分位划分，R 在每个十分位的密度 ≥ 各十分位平均密度的 50%（未校准操作值）；
- (d) 清单 id 无重复。

参考集来自 V2，所以不需要、也不允许"剔除截断 shard 后迭代重算 R"。

**A-S6 检测限制**（必须写入输出与报告）[用户 1]：
- L1 截断点之后必须还有足够的 R 成员才能探测到：截断点晚于 R 中倒数第 5 个不同 bag 的最早成员时，T_bag < 5，检测不到。
- L2 只能探测 R 中的常见查询。稀有查询（多数单例 bag）不在 R 中，它们的缺失看不见。
- L3 证据只含 ≥200 bp 的过滤后命中。
- L4 数据库一侧的损失（例如某个基因组没有进入遗留数据库）很少使常见查询在整个 shard 中消失，A-S 基本看不到，只有 A-D 能发现。
- L5 不假设独立性，所以 `SUSPECT_*` 可能出现在完整的 shard 上；结论以 A-D 为准。
- L6 没有足够的 `VERIFIED_V2` 参考 shard 时，A-S 不作任何判定。
- L7 证据已被删除的 8 个 shard（A-D5）判 `NOT_AUDITED`。

**A-S7 输出。** 原子写出（临时名 + `rename`），包含：
- G-4 字段；
- 全部输入（含参考提交记录）的 sha256；
- 全部阈值及其标注；
- |R|、R 的不同 bag 数；
- 每个 shard 的全部判定量与判词；
- A-S6 的限制说明。

**A-S8 测试**（合成数据，登录节点小样本）。A-S：
- (a) 在 34% 处前缀截断，且已读部分有正常缺席 → `SUSPECT_TRUNCATION`；
- (b) 共同尾部截断：64 个遗留 shard 中有 7 个在最后 1% 之前截断，参考集来自独立的 V2 → 这 7 个都判 `SUSPECT_TRUNCATION` [R1-7]；
- (c) 同一 bag 的 5 个 R 成员都在尾部缺席，其余正常 → T_bag = 1，不判截断 [R1-6]；
- (d) 全程散布缺席 30% → `SUSPECT_OTHER`；
- (e) gzip 损坏 → `CORRUPT`；
- (f) 证据 sha256 与当前文件不符 → `NOT_AUDITED`；
- (g) H 为空、清单外 id、参考 shard 不足、R 的 bag 数不足、B 分母为 0 → 各自按 A-S3、A-S5 判定；
- (h) 每个 shard 恰得一个判词；
- (i) 输出中不出现 A-S1 禁用的词。

---

## 不在本 SPEC 范围内

- Stage E：侧翼提取、占位 / 空位验证、ISCompare。"已验证的新位点"只能由 Stage E 产生，另立 SPEC 项。
- events（独立转座事件）。
- 修改 239：它在基线中冻结，遗留 shard 与 V2 都依赖它。
- 试点与生产运行的资源估算：流程第 3 步，单独提交用户批准。

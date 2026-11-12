# fna_based_mgefinder_project

## 审查循环（GPT-6.1-Sol 审查）
每完成 docs/SPEC.md 中的一项：
1. 运行 `CODEX_REVIEW_TIMEOUT=2400 codex-review "<项名>"`，用 run_in_background 在后台运行（不受工具 900 秒上限约束；2026-11-11 用户批准）；核心判定逻辑用 `CODEX_REVIEW_TIMEOUT=2400 codex-review "<项名>" review-hard`
2. 读 .review/review.md 首行：
   - PASS → git commit，message 注明 "reviewed: N rounds"
   - FAIL → 修复你认同的问题；不认同的写理由进 docs/DECISIONS.md；进入下一轮
   - ERROR → 立即停止，向我汇报，不计入轮数
3. 同一项最多 3 轮；仍 FAIL 则停下，汇报分歧，等我
4. 禁止在循环中 sbatch、scancel 或修改正在运行作业所用的脚本；需要算力或要改在用脚本时，停下问我
5. 每个作业的输出须记录提交时的 git commit hash
6. 登录节点只跑小样本测试

### 审查基础设施（2026-11-11 用户批准）
- wrapper v4：`~/.local/bin/codex-review`，仓库副本 `tools/codex-review`，两者须逐字节一致。测试：`TMPDIR=<scratch> bash tests/test_codex_review_wrapper.sh`（封闭测试，用桩程序，不调用真实 codex）。
- 每项的必需材料在 `docs/review_deps.txt` 中声明。wrapper 逐个检查，任一文件缺失、为空、为二进制、超过 200,000 字节或被 .gitignore 忽略，直接输出 `VERDICT: ERROR`，不调用审查者。
- `docs/review_deps.txt` 中的 `DIFF <旧文件> <新文件>` 行：两侧文件都按必需材料检查，但只发送 `diff -u`，旧文件不全文发送。用于上一轮被审的 SPEC 快照，避免把历次存档整份发给审查者。
- 例外：M1 模块的审查用 gpt-6-astra（profile `review-astra`，推理强度 high；用户 2026-11-11 指定），命令为 `CODEX_REVIEW_MODEL=gpt-6-astra CODEX_REVIEW_TIMEOUT=2400 codex-review "M1 Arm B" review-astra`，wrapper 同样核对 profile 与实际模型。
- 审查模型固定为 gpt-6.1-sol（`review` 与 `review-hard` 两个 profile 都是），两个 profile 都设 `model_reasoning_summary = "concise"`（2026-11-11 用户批准，永久生效），使长时间推理期间连接上有数据流出。wrapper 运行前核对 profile，运行后核对 codex 日志中的实际模型，不符即 ERROR、本轮结果作废。
- 每次运行（含 ERROR）自动存档到 `docs/reviews/<slug>/rNN.md`（审查文本，加输入 sha256、git 快照、逐文件 sha256、模型）和 `rNN.SPEC.md`。原始输入和 codex 日志只留在本地 `.review/history/`。
- `rNN` 是**运行序号**（含 ERROR 和 INCOMPLETE），不是审查轮数。轮数看存档中的 `review_round` 字段，只计 `completed_round: yes`（PASS / FAIL）的运行。
- wrapper 内部超时默认 2400 秒（`CODEX_REVIEW_TIMEOUT`），wrapper 在后台运行。启动 codex 之前，先把本次尝试登记为 `INCOMPLETE`；超时或收到 TERM/INT/HUP 时，只终止本次尝试的进程组，保留部分输出，存档为 ERROR。每次尝试使用独立目录 `.review/attempts/<attempt_id>/`，绝不沿用上一次的结论。

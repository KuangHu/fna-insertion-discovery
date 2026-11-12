# DECISIONS

审查循环中的处理记录：不认同的意见及理由、被推迟的意见、用户特许的处理。

## 项 0（V0 可信清单）第 5 轮（运行 r05，FAIL，2 条），2026-11-11

依据用户第八次决定第 3 点："PASS，或 FAIL 但不影响清单内容正确性 → 修复、记入 DECISIONS.md、提交代码，并提交 V0 作业。FAIL 且影响正确性 → 停下等我。"

两条意见都接受、都已修复。修复**未再送审**（第 5 轮是特许的最后一轮）。

| # | 意见 | 是否影响清单内容正确性 | 修复 |
|---|---|---|---|
| 1 | `check_dir()` 只查 `genome_reads.tsv` 的内部一致性；把某行 `0 / -` 改成 `2 / EIO` 仍能通过 `check`，与 `v0_result.json` 矛盾 | 否。`genome_reads.tsv` 是 I/O 经历的记录，不属于基因组可信清单，不进入 `genome_manifest.index`，也不进入入库的可信配置；清单本身的篡改由 index 与逐项检查发现（既有测试） | `check` 子命令（已完成的运行）要求 `v0_result.json` 存在且 `verified` 为真，逐个核对 5 个输出文件的字节数与 sha256，并核对 `io_rereads` 汇总与 `genome_reads.tsv` 一致。build 内部的检查在写 `v0_result.json` 之前运行，不要求它。新增 3 个测试（一致性篡改被摘要发现、汇总被篡改、缺少结果文件），共 53 个，全部通过 |
| 2 | `V0_prep.sh` 在启动 248 之前失败（版本命令、outfmt 7 自检、V0b）时不写失败记录，取不到节点名与作业号 | 否。失败的运行不产生 `v0_result.json`，不会被当作结果使用；这条只关系失败时的诊断记录 | `set -E` 加 ERR trap 记录失败的行与命令；EXIT trap 在非零退出时写 `$RUN/FAILED`（行、命令、节点名、SLURM 作业号、提交时的 commit；248 失败时指向 `manifest/FAILED`）。运行目录已存在时在建立 trap 之前退出，不写入别人的目录。用替换了路径的副本在 scratchpad 中实测：BLAST 路径不存在时写出 `FAILED`，含行号与节点名；运行目录已存在时原有 `FAILED` 不被覆盖 |

SPEC V0-7 相应补写这两点（标 `[第 5 轮第 n 条]`）。

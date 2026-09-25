---
name: fna-ins-discovery-bench60-result
description: "Stage 60 result 2026-08-27 — Arm A recall vs ISEScan tracks copy number, is 100% on IS110, and 0% on low-copy families"
metadata: 
  node_type: memory
  type: project
  originSessionId: f0e46982-bc7a-4c4f-8e44-dd8e6ec7df19
  modified: 2026-08-28T02:29:36.560Z
---

First real stage 60 run (SLURM job 25306752, ISEScan on 20 of the 200 E. coli, 0 failures, 10–25 min/genome). Outputs at `/global/scratch/users/kh36969/fna_ins_discovery/bench60_{contingency,per_is,novel_pool}.tsv`.

Contingency: ArmA+/ISEScan+ = 95 families; **ArmA+/ISEScan− = 55 families** (the discovery pool); ArmA−/ISEScan+ = 502 multi-copy + 314 singleton IS calls.

Recall on multi-copy **complete** IS: **64.3%** (541/842). Per-copy 53.8% all multi-copy; partial IS only 12% (expected — Arm A needs intact multi-copy structure). Family-level (per genome × ISEScan family, did Arm A hit ≥1 instance): **65.0%** (76/117).

**The finding: recall is a monotone function of copy number, not length.** All families sit in Arm A's 300–5000 bp window, so length explains nothing.

| ISEScan family | family-level recall | mean ISEScan ncopy |
|---|---:|---:|
| IS5, IS4, IS3, IS200/IS605, **IS110**, IS91 | **100%** | 7–18 |
| ISAS1 | 94% | — |
| IS1 | 64% | — |
| IS66 / IS30 / IS21 | 33–60% | — |
| IS6 / IS256 | **0%** | 2.7–4.5 |
| ~~ISNCY~~ | *excluded — not a real family* | — |

**IS110 is 100% (6/6 family units, 43/43 copies)** — the target family, and the one ISEScan is structurally weakest on.

**ISNCY must be excluded from the accounting — it is not a family.** ISEScan's "not classified yet" catch-all lumps unrelated elements under one label. Verified directly: extracting the ISNCY instances from each of 3 genomes and aligning them all-vs-all (`minimap2 -x asm5 -c -N 200 -p 0 --secondary=yes`) gives **zero mutually aligned pairs** in every genome, while IS3 in the same genomes gives 109–482 pairs at 98.7–100% identity. So Arm A's 0/20 on ISNCY is correct behaviour, not a blind spot. (A first attempt with `-x asm20 --secondary=no` showed ~0 pairs for IS3 too and was simply the wrong preset for ~1 kb fragments — do not reuse it.)

Corrected headline, multi-copy + complete, ISNCY excluded: **family-level 78.4% (76/97)**, per-copy **69.3% (541/781)**. With ISNCY wrongly included these read 65.0% and 64.3%.

Misses are real, not threshold artifacts: **758 of 816** missed calls have *zero* overlap with any Arm A copy; only 11 fall in the 0.25–0.50 band below `--min-overlap`. Every one of the 20 genomes produced Arm A families (4–13 each), so this is not a pipeline failure.

Novel pool (55 families): median length ~3.0 kb, notably longer than the recovered known IS (1,314 bp). Verdicts 21 `REPEAT_NO_BOUNDARY` (likely rrn), 17 `STRONG_MOBILE_CANDIDATE`, 16 `MOBILE_CANDIDATE`, 1 `SEGMENTAL_DUP_LIKE`.

P(Arm A detects | ISEScan ncopy), family-level, complete IS, ISNCY excluded — a clean monotone curve that defines Arm A's operating regime:

| ncopy | 2 | 3 | 4 | 5 | 6–9 | >=10 |
|---|---:|---:|---:|---:|---:|---:|
| recall | 10.0% | 36.8% | 73.3% | 93.3% | 95.2% | 88.9% |

With ISNCY included the n=5 bin falsely dips to 46.7% (ISNCY alone is 0/15 there), which is what first exposed the labelling problem.

**Why:** this defines Arm A's operating regime quantitatively — it is a high-copy engine, ~90%+ from 5 copies up and near-useless at 2 — and validates the IS110 design claim.
**How to apply:** the 17 STRONG novel families still need stage 50 (barrnap + ORF fraction + HMM) before any of them can be called a discovery. If low-copy recall matters, the lever is Arm A's `min_copies: 3`, not the length window.

See [[fna-ins-discovery-stage60]], [[fna-ins-discovery-state]].

---
name: fna-ins-discovery-union-recall
description: Recall(A u B) = 85.7% family-level on the 20-genome benchmark — meets the target, but Arm B rescues less than hoped
metadata:
  type: project
---

The user's 第二步, run 2026-08-27 via new `scripts/62_union_recall.py`. Arm B was run on the same 20 benchmark genomes: stage 10 gave 3 clusters (6/6/3, 15 of 20 genomes clustered), stage 30 gave **251 events + 103 oversized bubbles parked**.

Recall on multi-copy **complete** IS, clustered genomes only, ISNCY excluded:

| level | Arm A | Arm B | **A ∪ B** | n |
|---|---:|---:|---:|---:|
| **family (the honest one)** | 80.5% | 81.8% | **85.7%** | 77 units |
| per-copy (Arm B optimistic) | 69.7% | 89.3% | 92.7% | 590 calls |

**Meets the user's 85–90% target at family level**, and the ArmA+/ISEScan− pool stayed clean (see [[fna-ins-discovery-stage50-result]]).

**Critical caveat, now printed by the script itself.** Arm A is matched to ISEScan *positionally*; Arm B *by sequence*, because bubble coordinates are in backbone space and the backbone is often not the genome carrying the IS. One Arm B insert therefore matches every instance of that element in its cluster: **527 calls matched from only 63 distinct (genome, family) units**, and cl0002 matched 90 calls from just 15 available inserts. So per-copy Arm B is an upper bound, not a detection rate — compare arms at family level only.

**Arm B rescued only 4 family units Arm A missed** (IS1 ×2, IS30 ×2) — much less than hoped. The two arms overlap heavily at family level (80.5% and 81.8%, union only 85.7%).

**IS6 and IS256 were NOT rescued.** Still missed by both arms: 11 units — IS6 ×3, IS66 ×2, IS21 ×2, IS256 ×2, IS30 ×1, IS1 ×1. The hypothesis that Arm B would cover Arm A's low-copy blind spot did not hold for these families on this data.

**Why:** this is the measurement that decides the two-arm architecture's value, and it half-confirms it — the union clears the target, but the arms are more correlated than the design assumed, so the residual misses are a shared blind spot rather than complementary ones.
**How to apply:** to find out why IS6/IS256 survive both arms, check whether they sit in the unclustered 5 genomes, at contig breaks, or in bubbles that exceeded `max_insert`. Only 15 of 20 genomes were clusterable, which caps Arm B a priori.

See [[fna-ins-discovery-bench60-result]], [[fna-ins-discovery-state]].

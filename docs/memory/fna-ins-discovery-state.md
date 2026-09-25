---
name: fna-ins-discovery-state
description: "Status of the FNA-only mobile-element pipeline as of 2026-08-27 — current numbers, what is stale, what never ran"
metadata: 
  node_type: memory
  type: project
  originSessionId: f0e46982-bc7a-4c4f-8e44-dd8e6ec7df19
  modified: 2026-08-28T00:07:01.220Z
---

Reference-free insertion-element discovery from assemblies only (no reads/BAM). Arm A = within-genome multi-copy, Arm B = between-genome empty/filled, unioned at stage 40, annotation LAST (stage 50). Primary catalogue 300 bp–5 kb.

**The architecture is settled and already built.** Do not re-derive or re-propose it; see [[fna-ins-discovery-dont-rebuild]]. `README.md` holds the rationale and is reliable EXCEPT for the Arm B and stage 40 numbers below, which it predates.

Current numbers (2026-08-27, supersede the README):
- Arm A, 200 genomes: 1,648 families, ~1.2 s/genome. IS110 recovered at 93% STRONG, 0/13 events with a TSD. (README still accurate here.)
- Arm B, mini3 cluster: **141 events, 33 with TSD** — README says 183/43.
- Stage 40 → `stage40_v2/`: **272 elements — 19 BOTH_ARMS / 198 ARM_A_MULTICOPY / 55 ARM_B_EMPTY_FILLED**. README says 310 / 21 / 202 / 87, and `stage40_test/` is the stale run behind it. E00000 remains the IS1-like anchor: 768 bp, 21 copies, 124 genomes, 9 bp TSD.

Changed in the code on 2026-08-27:
- Stage 30 no longer deletes oversized bubbles. They go to `armB_large_events.tsv` + `armB_large_inserts.fna` (ids `<cid>.L%06d`), and `armB_events.tsv` gained `parent_large_event_id`. On mini3 this recovered **42 bubbles** (median 15 kb, range 5.5–66 kb) — the bin holding the known IS110-inside-33kb/36kb-cargo cases.
- `architecture` is now categorical (`asymmetric_candidate` / `clean_or_graph_slack`), magnitude only in `junction_span`. It used to embed the span (`asymmetric_candidate_+1066`), making every asymmetric event its own class and breaking grouping in stage 61. mini3 now 64 asymmetric / 77 clean.
- Stage 60 rewritten as a contingency table; see [[fna-ins-discovery-stage60]]. **It has now been run** — results in [[fna-ins-discovery-bench60-result]].

Open items, roughly in priority order:
1. **The 70 read-validated junctions have never been located.** Stage 61 has only ever run on synthetic data in `bench61_synth/`. This blocks the FNA-first decision and only the user can supply it.
2. **Stage 32 (nested sub-5 kb search inside the retained large bubbles) is designed but NOT written.** This is what makes keeping the 42 bubbles pay off.
3. **Arm B cl0000 (9 genomes) returns 0 bubbles** while mini3 (3 genomes) returns 248; its `graph.gfa` is absent from the output dir. Undiagnosed.
4. ~~Stage 50 has never been run~~ — **run 2026-08-27**, rewritten to report separate evidence axes. Result in [[fna-ins-discovery-stage50-result]]: no novel MGE in the 55-family pool. Arm A parameters are LOCKED by user decision: `300 <= L < 5000`, `min_copies = 3`, TSD optional — do not tune them.
5. ~~第二步 Arm B union recall~~ — **done**, result in [[fna-ins-discovery-union-recall]]: A u B = 85.7% family-level, target met, but IS6/IS256 rescued by neither arm. Next open question is why those survive both.
6. **`target_site_recorder` still not written** — design in [[target-site-recorder-design]]. Real empty-allele data now exists in `armB_bench/stage30/*/` to build against.
6. ISEScan has been run on **20** of the 200 genomes (`isescan_run/tsv/`). Extending the 2x2 to all 200 is a resubmit of `slurm/isescan_array.sh` with a bigger list.

**Why:** the 2026-08-27 chat was lost with no transcript and no memory, and none of this is derivable from the README.
**How to apply:** before quoting any Arm B or stage-40 number, check file mtimes against `scripts/30_armB_graph.py`.

See [[fna-ins-discovery-paths]], [[long-jobs-need-slurm]].

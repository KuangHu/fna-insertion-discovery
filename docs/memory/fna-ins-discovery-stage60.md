---
name: fna-ins-discovery-stage60
description: "Stage 60 is a 2x2 contingency table against ISEScan, not a recall number — and it must not grade its own novel cell"
metadata: 
  node_type: memory
  type: project
  originSessionId: f0e46982-bc7a-4c4f-8e44-dd8e6ec7df19
  modified: 2026-08-28T00:07:13.905Z
---

The user asked for stage 60 to be a contingency table, explicitly not a single recall figure, because ISEScan is not truth: it finds IS by transposase pHMM + terminal repeats and so carries the family bias the pipeline exists to avoid.

```
                 ISEScan +              ISEScan -
  Arm A +   known IS recovered    ANNOTATION-FREE CANDIDATES  <- the point
  Arm A -   sensitivity failure   background (not enumerable, reported "-")
```

Design decisions in `scripts/60_benchmark_armA_vs_isescan.py` (rewritten 2026-08-27), each load-bearing:
- Genomes lacking either side are dropped. Counting an Arm A family from a genome ISEScan never saw would inflate the ArmA+/ISEScan− cell — the one cell that matters.
- ArmA−/ISEScan+ splits into multi-copy (a real failure) vs singleton (a declared blind spot: Arm A needs ≥3 copies at distinct loci and never claimed singletons).
- It deliberately does **not** grade the novel pool. Deciding whether an ArmA+/ISEScan− family is real needs the stage 50 screens (barrnap rRNA, ORF fraction, transposase HMM); doing it inside the benchmark would put annotation back in front of discovery. Stage 60 only isolates and describes the pool into `_novel_pool.tsv`.

**Why:** ArmA+/ISEScan− is expected to be enriched, not wrong — IS110/IS1111 has neither TIR nor TSD, the two features ISEScan keys on, so an annotation-free method should win there. Reading that cell as false positives would invert the result.
**How to apply:** outputs are `<out>_{contingency,per_is,novel_pool}.tsv`. Feed `--isescan-glob 'isescan_run/tsv/*.tsv'`. Grade the novel pool with stage 50, separately.

See [[fna-ins-discovery-state]], [[long-jobs-need-slurm]].

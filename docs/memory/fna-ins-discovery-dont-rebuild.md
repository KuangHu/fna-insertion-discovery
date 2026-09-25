---
name: fna-ins-discovery-dont-rebuild
description: On this project the user re-states the whole pipeline design as context — implement only the delta, not the recap
metadata:
  type: feedback
---

The user writes long messages that restate the entire FNA-only pipeline design (skani → PanGraph → NucDiff/SVIM-asm → ISEScan, annotation last, TSD not a filter, tiered benchmark, contig-break QC). Most of it is already implemented. When they say it again, find the one or two genuinely new asks buried in it and do only those.

Corrected mid-turn on 2026-08-27 with: "and you have already do pipeline setup and test, dont follow all" — after I began re-implementing the tiered stage-61 benchmark that was already fully built.

**Why:** re-implementing settled design burns the turn and risks overwriting working code.
**How to apply:** diff the message against `README.md` and the current scripts first; state briefly what is already done, then implement only the remainder. In that message the real deltas were just two — keep >5 kb bubbles instead of deleting them, and make the ISEScan comparison a 2×2 rather than a single recall number.

See [[fna-ins-discovery-state]].

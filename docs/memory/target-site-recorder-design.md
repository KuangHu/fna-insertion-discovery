---
name: target-site-recorder-design
description: Design decision — target_site_recorder is a separate module anchored on homologous anchors, never on a POS
metadata:
  type: project
---

User decided (2026-08-27) to build a **separate `target_site_recorder` module**, kept out of the Arm A / Arm B discovery logic. Not yet written.

The motivating case: the empty allele is often NOT `left_flank + right_flank`. It can carry a short interval that the insertion replaced —

```
filled:  [flank_L][      INSERTION      ][flank_R]
empty:   [flank_L][GTC][flank_R]
```

Reconstructing the empty site as a plain concatenation destroys that `GTC`, which may be exactly the target-site grammar a future RNA-guide model needs.

**Core object is not a POS.** It is `left homologous anchor + variable interval + right homologous anchor`, with both alleles expressed in that frame. So the transformation recorded is `GTC -> IS_ELEMENT`, never `empty -> IS_ELEMENT`.

Requirements:
- Record `empty_interval` and `filled_interval` separately, plus anchor identities and lengths.
- `event_class` from sequence transformation only, no biological assumption: `pure_insertion`, `replacement`, `local_deletion_plus_insertion`, `target_duplication`, `target_retention`, `asymmetric_replacement`, `complex`, `unresolved`.
- Support **multiple empty alleles** — store all observed empty sequences with counts, a consensus and an entropy, rather than picking one genome as truth. Independent genomes agreeing on the same 3 bp is strong evidence against assembly error.
- Store pre-insertion target **context windows** (50/100/200 bp each side of the empty interval), because the filled element's own flanks may have been altered by recombination.
- Never hard-code "empty" as ancestral: emit `allele_A` / `allele_B` + `ancestral_state_confidence`, polarity only from an outgroup.

**Why:** naive `left+right` reconstruction silently loses the pre-insertion target, and a single POS cannot express an asymmetric junction where the left and right breakpoints differ (already observed in the Nanopore audit: dL=+59/dR=0, dL=-56/dR=0).
**How to apply:** build it downstream of Arm A/B as a refinement over candidate events; it consumes anchors, not coordinates. Best built once Arm B has run on the 20 benchmark genomes so there is real empty-allele data to test against.

See [[fna-ins-discovery-state]].

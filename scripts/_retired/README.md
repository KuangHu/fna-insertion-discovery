# Retired 2026-09-02 by docs/SIMPLIFICATION.md

Kept for reproduction of historical numbers, NOT on the pipeline path.

| script | why retired |
|---|---|
| `80/81/82/83/85_*` | Layer 3 polarity chain — answers ancestry, not target site |
| `84_polarity_ingroup_refiner.py` | focal clade exists only to serve Layer 3 |
| `88_anchor_alignability_sweep.py` | outgroup-distance study, Layer 3 machinery |
| `89_layer3_regression.py` | gates a module no longer on the path |
| `92_structural_polarity.py` | `junction_ambiguity_bp` reproduces it exactly (76/76, lengths identical, conflicts 49→0) |
| `30b_armB_pangraph.py` | second graph engine, minigraph path validated |
| `31_armB_refine.py` | `svim-asm`/`nucdiff` absent; the 200 bp pad means Layer 2 relocates the junction from sequence (verified 100/100) |
| `40_merge_catalog.py` | unions the two arms; `99_mobility_join.py` joins them |
| `10_skani_cluster.py` | single-linkage clusters are not clades; `91_seed_panel_sampler.py` replaces it. Needed to reproduce the 234-locus series |

**Before restoring any of these, read `docs/SIMPLIFICATION.md` §3 and §7** —
their removal was a scope decision with the cost stated, not an oversight.

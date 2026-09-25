---
name: fna-ins-discovery-stage50-result
description: Stage 50 first real run 2026-08-27 — the 55-family ArmA+/ISEScan- pool contains no novel MGE; it is Rhs toxins and prophage
metadata: 
  node_type: memory
  type: project
  originSessionId: f0e46982-bc7a-4c4f-8e44-dd8e6ec7df19
  modified: 2026-08-28T02:57:46.423Z
---

First real stage 50 run, on the 55 `ArmA+/ISEScan-` families from [[fna-ins-discovery-bench60-result]]. Evidence axes: barrnap, prodigal ORF fraction, ISEScan's own ISfinder pHMM library (`$fnains_annot/bin/pHMMs/clusters.faa.hmm`, 266 profiles), and full Pfam-A with `--cut_ga` (43 s for 152 proteins). Output `stage50_novel/element_annotation.tsv`.

| class | n |
|---|---:|
| known_non_MGE_repeat | 29 |
| known_other_MGE (prophage) | 19 |
| ribosomal_repeat | 3 |
| unclassified | 2 |
| segmental_duplication | 1 |
| known_IS_missed_by_ISEScan | 1 |

**Result: no novel mobile element in this pool.** The 17 `STRONG_MOBILE_CANDIDATE` families resolve to 12 known-non-MGE, 3 rrn, 2 prophage — zero novel. The pool is dominated by **Rhs / YD-repeat polymorphic toxin systems** (`RHS_repeat`, `TEN_YD-shell`, `DUF6531`), which are genuinely repeated with divergent flanks and therefore look exactly like Arm A's target, but are not transposons. Second largest component is prophage structural genes (`Phage_GPD`, `Phage_base_V`, `Gp5_trimer_C`, tail/terminase).

Both classes are invisible to ISEScan *by construction*, so ArmA+/ISEScan− being large was never evidence of novelty. The one `known_IS_missed_by_ISEScan` (`GCA_002057355.1.A0006`) is a composite — IS3 pHMM hits *and* phage tail domains in a 4.8 kb `REPEAT_NO_BOUNDARY` — i.e. an IS3 inside a prophage, not a clean missed IS.

3 rrn families were promoted to `STRONG_MOBILE_CANDIDATE` by discovery, not caught by `REPEAT_NO_BOUNDARY`. That residual false-positive mode is precisely what stage 50 exists to catch.

Two bugs fixed while doing this, both pre-existing:
- **`longest_orf_bp` was 180 for every element.** Prodigal wraps protein FASTA at 60 aa and the parser took the max single *line* (60 × 3). Now read from prodigal's own `# start # end #` header coordinates. Median `orf_frac` went 0.05 → 0.59.
- The single collapsed `class_call` if/elif chain let one axis mask another. Every axis is now its own column; the class is derived and never a substitute for reading them.
- Added `known_non_MGE_repeat`: a confident non-DUF Pfam identity that is not an MGE domain still *identifies* the element. Without it, 10 Rhs families were being reported as `unknown_mobile_element_candidate` and would have inflated the discovery pool with a well-characterised family.

**Why:** this is the first evidence on whether annotation-free discovery yields anything ISEScan cannot, and on this subset the answer is no — the honest read is that the pool is real biology that ISEScan ignores, not new biology.
**How to apply:** the pool came from only 20 genomes. The full Arm A run has 1,648 families over 200; extending ISEScan to all 200 gives a far larger pool and is the way to actually test for novelty. Do not quote "55 novel families" as a discovery result.

See [[fna-ins-discovery-stage60]], [[fna-ins-discovery-state]].

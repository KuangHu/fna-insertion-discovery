# FNA-only Mobile Element Discovery — Full Project Summary

**Status as of 2026-08-31.** Layers 1–4 are built and closed; the blocker is a
resource-allocation question, not an algorithmic one.

---

## 1. What this project is

Reference-free, annotation-free discovery of insertion elements by comparing
assembled bacterial genomes. No SRA download, no read alignment, no BAM.

The founding constraint, which every design decision follows from:

> **Never start from a transposase.** Discover structural insertion alleles
> first; ask what they are only afterwards.

Starting from an IS library can only find what the library already contains.
This project inverts that ordering, so it can in principle find noncanonical or
RNA-guided elements that no HMM would have retrieved.

**Verified, not assumed.** An audit of every executable line in the discovery
chain (stages 10, 20, 30, 70, 80, 81, 83, 84, 86) found the only external
binaries invoked are `skani`, `minimap2`, `minigraph`, `gfatools`, `mafft`,
`iqtree`, `samtools`. Zero calls to `hmmsearch`, `phmmer`, `isescan`,
`prodigal`, `barrnap`, or Pfam. The only mentions of "transposase"/"IS" in those
files are docstring prose stating the constraint. Candidate loci are minigraph
bubble IDs; the sole `is_` column in the event table is `is_inversion`. Genome
selection was the first 20 alphabetically, not IS-content driven.

**One honest qualification:** the <5 kb scope cap is IS-*motivated*. It is
applied as a pure length filter with no sequence model, so it adds no family
bias among sub-5 kb elements, but it does bound what can be found. Oversized
bubbles are parked in `armB_large_events.tsv`, never discarded.

---

## 2. Architecture

```
Layer 1   candidate discovery          Arm A (within-genome multi-copy)
                                       Arm B (between-genome graph bubbles)
   ↓
Layer 2   allele reconstruction        L_anchor | allele | R_anchor      LOCKED
   ↓
Layer 3   polarity                     local-flank genealogy → P1 tier   CLOSED
   ↓
Layer 4   pre-insertion target         combine L2 structure + L3 direction  LOCKED
   ↓
Layer 5   annotation                   ISEScan / Pfam / HMM — strictly last
```

The layers are deliberately separable. Layer 2 answers *what changed*; Layer 3
answers *which direction*; Layer 4 only combines them and introduces no new
boundary estimate.

---

## 3. Development history, in order

### 3.1 Arm A — within-genome multi-copy discovery

Self-alignment finds repeats present at ≥3 distinct loci with high internal
identity and divergent flanks. On 200 *E. coli*: 1,648 families, 1.2 s/genome,
64% `STRONG_MOBILE_CANDIDATE`.

**Bug found — minimap2 masking (major).** `self_align()` ran minimap2 without
`-f`, so the default masking of the most frequent 0.02% of minimizers applied.
That default is correct for read mapping and **exactly backwards here**: the
highest copy-number repeat in a genome is the one whose minimizers are most
frequent, so the default was deleting precisely the elements Arm A exists to
find.

On `GCA_011404755.1` (60 IS1 copies):

| | alignments genome-wide | covering an IS1 locus |
|---|---:|---:|
| default `-f` | 138 | **0** |
| `-f 0` | 3,442 | **57** |

Fixed with `-f 0`. Effect across 200 genomes:

| | pre-fix | post-fix |
|---|---:|---:|
| families | 1,648 | 1,610 |
| copies | 8,563 | 9,847 |
| families ≥20 copies | 18 | 57 |
| families ≥50 copies | **0** | **8** |
| max copies | 37 | **109** |

Before the fix, *no family in 200 genomes had ≥50 copies*.

### 3.2 Arm B — between-genome graph bubbles

minigraph + `gfatools bubble` over an ANI cluster; empty/filled alleles fall out
of the bubble directly. Oversized bubbles (>5 kb) are parked rather than
dropped — 42 recovered on a 3-genome test, median 15 kb, including the bin that
holds the known IS110-in-33 kb-cargo cases.

### 3.3 Layer 2 — allele reconstruction

Cuts a 500 bp anchor each side of a candidate locus, maps both anchors into
every genome, and extracts the interval **between** them verbatim — including
when it is 0 bp. Alleles are labelled A, B, C… by carrier count. Nothing is
called "empty", because polarity is not decidable from alleles alone.

Beyond the original design it adds `decompose()`: explain the long allele as the
short one plus an insert, via longest common prefix/suffix. A real record:

```
locus GCA_002057355.1.A0000.c5    6 genomes, 4 alleles
  shortest A =   40 bp   (3 genomes)
  longest  D = 1319 bp
  lcp=37  lcs=3  inserted=1279 bp  ambiguity=0  lost=0
```

The element landed **37 bp into a 40 bp target interval**. A naive
`left_flank + right_flank` reconstruction destroys all 40 bp. The dominant
inserted length across the run was **1279 bp — exactly the IS110 insert size Arm
B measured independently by a different route**.

**Bugs found:**
- When all alleles have equal length, `min`/`max` by length return the *same*
  sequence, and the decomposition silently reported a perfect clean insertion
  for 23 loci that were substitution-only.
- Negative "ambiguity" was being reported as ambiguity when it means the
  opposite (target bases destroyed).

**Exact matching was too brittle.** At strain-level divergence a single SNP
truncates the prefix/suffix match. Median `(lcp+lcs)/len(short)` was 1.00 for
clean decompositions but **0.26** for failures. Replaced with an
alignment-based `tolerant_decompose()` behind a strict gate (coverage ≥0.85,
identity ≥0.95, dominant insert ≥3× the next indel). Scorable events **61 → 135**
with the originally-clean 61 unchanged.

**Placement bug.** `place_locus()` maximised raw `nmatch` under a 50 kb cap, so
on repeat-rich loci the two anchors could land on different copies of a repeat —
fabricating a **42,367 bp** allele. Rewritten as anchor-pair scoring with a
length-normalised score, a 7 kb interval cap, a seed-genome localisation prior,
and a runner-up margin; ambiguous placements are reported, not forced.

### 3.4 Layer 2 validation

Benchmarked against Arm B on 251 events across 3 clusters (loci seeded with a
200 bp pad so the junction had to be relocated from sequence, not handed over).

| | value |
|---|---:|
| loci resolved | 234 |
| inserted-length exact, exact method | **91.7%** |
| inserted-length exact, tolerant method | 74.0% |
| all decomposed | 82.0% |
| reverse-complement invariance | **96.2%** |
| catastrophic placements | **0** |

**A bug in the audit nearly produced a false 19.2% error rate.** Arm B stores
its insert in the backbone's orientation, Layer 2 in the pre-event carrier's;
forward-only comparison returns ~0.50 identity — the score two unrelated
sequences get. **6 of 9 apparent errors matched `revcomp` at 1.000 identity.**
All sequence comparison in this project is now orientation-agnostic.

Of 26 residual disagreements, **91.7% are Arm B extent imprecision** (Arm B
simply drew a wider box). The 2 unexplained are recorded as
`validation_flag = assembly_crossmethod_discordant` — 0.9%, no repeated pattern.

**LOCKED** with 8 frozen gates in `73_layer2_regression.py`, all passing.

### 3.5 Layer 3 — polarity

The hard part, and where most of the project's effort went.

**Majority voting was tested and destroyed.** Re-running outgroup selection at
ANI ceilings 99.0 / 98.5 / 98.0 / 97.5:

| max-ANI | callable ≥2 | unanimous | majority | split |
|---|---:|---:|---:|---:|
| 99.0 | 175 | 75 | 65 | 35 |
| 98.5 | 155 | 65 | 54 | 36 |
| 98.0 | 127 | 51 | 48 | 28 |
| 97.5 | 108 | 46 | 42 | 20 |

**55.7% of loci change their supported allele across the grid.** Including loci
that were *unanimous* at 99.0 and flipped to a different real allele at 98.5.
Close relatives share the derived event; a unanimous vote among them is not
ancestry.

**Replaced with local-flank genealogy.** Per locus: 1 kb left and right flanks
(variable interval excluded), MAFFT alignment, IQ-TREE ML trees for L, R and
concatenated LR, Fitch parsimony for the discrete allele character on the LR
tree, UFBoot-backed support. Hard gates: focal-clade monophyly, L/R concordance,
≥2 informative outgroups, unique ancestral state.

**Discovery cluster ≠ polarity clade.** Stage 10 chains members in
(single-linkage), so the "cluster" is not a clade — cl0000's internal ANI spans
98.92–99.94 while its outgroups reach 98.97, i.e. **the bands overlap**. A
separate `84_polarity_ingroup_refiner.py` carves a complete-linkage focal clade
with an explicit separation gap. P1 went **4% → 28%**, monophyly failures 18 → 8.

**Two measurement bugs made every gap negative:** NCBI carries every assembly
twice (`GCA_x` and `GCF_x` at 100.000% ANI) and re-deposits strains at ~99.99%;
and the gap must be measured against the **selected outgroups**, not the whole
pool — with 13,027 genomes the ANI continuum is full, so a whole-pool gap is
unachievable by construction.

**Topology vs state separation.** A taxon has two independent qualifications:
`topology_informative` (its flanks place) and `state_informative` (it also
matches a known allele confidently). A distant outgroup is routinely the first
without the second. Previously such taxa were dropped from the tree entirely —
so far outgroups contributed nothing, not even rooting. Now they enter with
`state=None`. **NOVEL is topology-useful, not useless.**

**The result that matters:**

$$P(\text{P1} \rightarrow \text{different real allele}) = \mathbf{0.0\%}$$

Across every panel tested (NEAR / FAR / MIXED), shared P1 loci agree on the
ancestral allele — 19/19 and 25/25. 47% become *unresolved* under resampling,
which is the gates withdrawing a claim as evidence thins. Compare 55.7% for
majority voting on the same data.

> **Changing taxon sampling may make a P1 call disappear.
> It must never make a P1 call change its answer.**

**CLOSED** with 5 frozen gates in `89_layer3_regression.py`, all passing.

### 3.6 Layer 3 waterfall (234 loci)

| step | n | % | drop |
|---|---:|---:|---:|
| Layer-2 resolved loci | 234 | 100.0 | |
| ≥2 informative outgroups | 162 | 69.2 | −72 |
| flanks placed in enough taxa | 162 | 69.2 | −0 |
| focal clade locally monophyletic | 103 | 44.0 | −59 |
| no L/R genealogy conflict | 98 | 41.9 | −5 |
| unique ancestral state | 83 | 35.5 | −15 |
| bootstrap support pass → **P1** | **72** | **30.8** | −11 |

Anchor placement costs nothing (−0). The two bottlenecks are outgroup
information and local monophyly.

### 3.7 Layer 4 — pre-insertion target catalogue

Combines Layer 2 structure with Layer 3 direction; introduces no new boundary
heuristic. Two gates that **refuse rather than repair**: if Layer 3 says the
*long* allele is ancestral, the event is emitted as `deletion_candidate` with no
target (P1 qualifies ancestry, not insertion-ness); and the allele ID is
resolved back to the Layer 2 table and md5-checked before any sequence is
written.

From 72 P1 loci: 20 insertion targets, 14 deletion candidates, 12 complex
polarity, 17 non-insertion decomposition, 9 length-inconsistent.

**Downstream-ready = 20/72 (27.8%).** Quality: `target_bases_lost == 0` on
100%, junction ambiguity 0 on 55%, median inserted length 1072 bp with modes at
**777 bp (IS1)** and **1279 bp (IS110)** — both recovered structurally, with the
HMM applied only afterwards.

**A framing error caught mid-build.** The first version reported a 409 bp
"pre-event interval" as the target site. It is not — loci were seeded with a
200 bp pad, so the anchored interval is ~400 bp of ordinary flanking sequence
with the junction at offset `lcp` inside it. Context is now re-centred on the
junction and **re-extracted from the carrier genome**, recovering all 7
truncated loci (13 → 20). Layer 2's `target_context.tsv` was renamed
`anchored_interval_context.tsv` with a warning, and the Layer 2 lock re-verified
after the rename.

**End-to-end regression:** `pre[:offset] + insert + pre[offset:]` must
regenerate the observed derived allele. Median identity **1.0000**, min 0.9915
(residual is strain-level SNPs in retained flanks, since the two alleles come
from different genomes).

### 3.8 Novelty baseline

`87_novelty_positional.py` maps each insert back into a filled carrier and
intersects with that genome's ISEScan calls. **126/234 known (53.8%), 108
unmatched (46.2% upper bound on novelty).** Families: IS1 44, IS4 23, **IS110
20**, IS3 18.

**Three bugs, all inflating apparent novelty:** an earlier size-only heuristic
gave 58.6%; it `break`s at the first size-compatible call so family attribution
followed iteration order (which is why IS110 appeared absent); and **ISEScan
writes comma-separated content into `.tsv` files**, so a tab-assuming reader saw
506 intervals instead of 1,411 and reported 85.5%. Only **46.2%** is valid.

---

## 4. Scaling work

### 4.1 Data inventory

| | |
|---|---:|
| E. coli files on disk | 13,027 |
| **unique accession cores** | **7,723** (40.7% were GCA/GCF duplicates) |
| median contigs | 3 (84% ≤5 contigs — complete-level, not draft) |
| GenBank E. coli complete+chromosome | 8,958 |
| **already covered** | **7,720** |
| **Step 1 download delta** | **1,238 genomes (~6 GB)** |

Metadata downloaded (2.9 GB, all to scratch): GenBank + RefSeq assembly
summaries, GTDB bac120, PDG pathogen-detection clusters, ATB metadata, Pfam clan
tables.

### 4.2 Density — the pool is not sparse

All-vs-all ANI over 7,723 genomes, **59,034,423 pairs**:

| neighbours (ANI ≥99, AF ≥85) | genomes | % |
|---|---:|---:|
| 0 | 165 | 2.1 |
| 1–2 | 167 | 2.2 |
| 3–5 | 172 | 2.2 |
| 6–10 | 237 | 3.1 |
| 11–50 | 638 | 8.3 |
| **>50** | **6,344** | **82.1** |

**94.1% have ≥5 partners.** This overturned two earlier conclusions: the
PDS-based estimate (15%) understated by ~6× because SNP-surveillance clustering
is far stricter than ANI; and the benchmark's small 6/6/3 clusters were a
**sampling artefact** of drawing 20 genomes alphabetically from 7,723, not a
property of E. coli.

**Consequence: do not download more E. coli for density. The lever is
selection.**

### 4.3 Degree does not predict panel size

Panels grown as **cliques** (a genome joins only if it clears the threshold
against *every* member, not just the seed):

| | median seed degree |
|---|---:|
| panels reaching k=48 | 440 |
| panels falling short | 430 |

Essentially identical. `GCF_029916805.1` has degree 360 and stalls at k=10;
`GCA_029717345.1` has degree 663 and stalls at 23; a seed with degree 118
reached 48. **Median usable fraction of a neighbourhood: 0.104.** Any figure
derived from raw degree overstates capacity ~10×.

Panel quality where they form (73/80 reached target k): at k=24, **23 distinct
PDS clusters**, largest PDS fraction 0.04 — diversity scales with size, so these
are phylogenetically distributed panels, not clonal expansions.

### 4.4 Infrastructure lessons

- **`skani triangle` buffers all output.** A 12 h wall kill destroys 100% of the
  work. Measured: 4,400/7,723 queries after 10h47m, needing ~8 h more. Cancelled
  and resharded as 40 × `skani dist` blocks — ~4 h, restartable.
- **Compute nodes have outbound network** (NCBI 200 in 0.49 s). `lrc-xfer`
  refuses non-interactive auth, so batch jobs on `cf1` are the transfer path.
- **es1/es0/es2 are GPU partitions** — a CPU-only job pends forever on
  `QOSMinGRES`. Cost 8 h of wall-clock once. `cf1` is the reliable CPU choice.
- Compute-node `curl` is **7.61.1**, predating `--retry-all-errors`.
- Scratch holds 9.66 M inodes already — never unpack to per-genome files.

---

## 5. Current status

| layer | component | status | gate |
|---|---|---|---|
| 1 | Arm A / Arm B discovery | in use | — |
| 2 | `place_locus()` | **LOCKED** | `73_layer2_regression.py` (8 gates) |
| 2 | exact + tolerant decomposition | **LOCKED** | ″ |
| 2 | placement-ambiguity reporting | **LOCKED** | ″ |
| 3 | focal-clade definition | **LOCKED** | — |
| 3 | topology/state separation | **LOCKED** | `89_layer3_regression.py` (5 gates) |
| 3 | P1 polarity | **QUALIFIED** | ″ |
| 3 | outgroup architecture | **CLOSED** | ″ |
| 4 | pre-insertion target | **LOCKED** | end-to-end reconstruction ≥0.99 |
| — | nucleotide-resolution junction | **NOT CLAIMED** | needs read gold |

**Capability boundary.** Layer 2 claims: given a candidate locus and a cluster
of homologous assemblies, reconstruct the alternative alleles between homologous
anchors, preserve zero-length and retained target intervals, and decompose
high-homology allele pairs into insertion/replacement representations with
explicit placement and junction ambiguity. It does **not** claim a
nucleotide-resolved biochemical junction — Arm B's own `junction_span` has
median 35 bp, so no assembly-only comparison can certify ±2 bp. The 66.4%/88.2%
figures are *coarse concordance*, never nucleotide precision.

---

## 6. The current problem

The k = 6/12/24 saturation sweep (5 seeds × 3 k, nested panels, pre-registered)
returned **14/15 cells**. It did not produce a saturation curve.

**Yield (P1 / downstream-ready):**

| role | k=6 | k=12 | k=24 |
|---|---|---|---|
| high_coherence | 1/1 | 3/2 | 6/3 |
| high_pds_diversity | 11/7 | 15/12 | 25/18 |
| low_coherence_pass | 37/25 | 51/38 | **14/5** |
| many_unlabelled | 0/0 | 8/0 | *(running)* |
| typical | 7/6 | 9/8 | **38/33** |

Pre-registered ratios $Y_{12\to24}/Y_{6\to12}$: **0.50, 0.60, −1.27, 6.25**.
Median 0.55 nominally clears the "run k=48" bar, but **the median is meaningless
here** — one seed collapses, another jumps 4×. Applying the rule to it would be
exactly the post-hoc reasoning the pre-registration existed to prevent.

### The structural cause

**The focal clade does not grow with the panel.**

| role | k=6 | k=12 | k=24 |
|---|---|---|---|
| high_coherence | 3/6 | 5/12 | 5/24 |
| **low_coherence_pass** | 3/6 | 3/12 | **3/24** |
| high_pds_diversity | 3/6 | 3/12 | 4/24 |
| typical | 4/6 | 4/12 | 3/24 |

`84_polarity_ingroup_refiner.py` selects the best-gap complete-linkage clade,
and that clade stays at **3–5 genomes regardless of panel size**. So increasing
k adds near-external taxa and loci but **not** focal-clade members, while making
monophyly of a fixed 3-genome clade harder inside a larger, more divergent tree.

Measured on the collapsing seed: `local_focal_clade_nonmonophyly` goes
**17 → 30 → 162** (of 78/131/197 loci). Loci nearly triple; P0 swallows almost
all of them. That is a **changed inference regime**, not diminishing returns.

`typical` moved the opposite way (nonmonophyly 13 → 4, P1 9 → 38), so panel
growth is not uniformly harmful — it depends on whether added genomes tighten or
loosen the local genealogy.

### What predicts productivity better than k

ANI **gap** and **floor**. `many_unlabelled` (gap 0.17–0.45, floor 99.14) gave
0–8 P1 and **zero** ready targets at any k. `high_pds_diversity`
(gap 1.12–1.19, floor 99.98) yielded steadily across all k.

### Recommendation

**Do not run k=48.** Not because yield is saturated, but because k is the wrong
variable — more k on this design would amplify a known pathology.

The open question: **should the focal clade grow with the panel?** That is a
change to a LOCKED module and needs its own regression, not an ad-hoc tweak.

---

## 7. Open items

1. **Focal-clade/panel decoupling** — the current blocker (§6).
2. **The 70 read-validated junctions have never been located.** Stage 61 has
   only ever run on synthetic data. This is the sole route to a
   nucleotide-precision claim and the one input only the user can supply.
3. **Stage 32** (nested sub-5 kb search inside the parked large bubbles) is
   designed, not written.
4. **Arm B `cl0000` returned 0 bubbles** from 9 genomes while a 3-genome set
   returned 248 — undiagnosed.
5. **The 16 persistent `local_focal_clade_nonmonophyly` loci** are retained
   unresolved with full diagnostics, deliberately not rescued (re-drawing the
   clade to fit the tree and then reading ancestry off that tree is circular).
6. **`--allele-margin` frozen at 0.005.** Distant outgroups give smaller
   margins, but that means the sequence no longer distinguishes the alleles.
   A divergence-scaled margin needs truth calibration first.
7. **Rejected on evidence:** Pool O (sister species at 90–93% ANI sit below the
   95% anchor-placement cliff — 30.8% placement); more E. coli for density
   (94.1% already have ≥5 partners).

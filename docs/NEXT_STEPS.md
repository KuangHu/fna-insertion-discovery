# NEXT STEPS — dual-channel architecture, tiers, and the work queue

**Written 2026-09-01.** Supersedes the panel-size line of enquiry. Read
`docs/PROJECT_SUMMARY.md` first for how the pipeline got here.

---

> ### SCOPE OF THE 234-LOCUS BASELINE — applies to every number derived from it
>
> **The 234 loci (cl0000/cl0001/cl0002) come from `10_skani_cluster.py`,
> single-linkage.** §3.5 showed single-linkage clusters are NOT clades — cl0000
> spans 98.92–99.94 internally while its outgroups reach 98.97, two overlapping
> bands. Every downstream figure quoted from that set — **234 loci, 133
> decomposable, 117 clean, 85 TSD band, 83 M1, 77 usable, 72 P1, 34 catalogued,
> 20 targets** — inherits that construction.
>
> **K. pneumoniae and every later species run through
> `91_seed_panel_sampler.py`** (mutually comparable cliques at a calibrated
> floor), and §9b's gap/floor criterion is defined only for `91_`. **The 234
> series is therefore NOT comparable across species.** For cross-species
> comparison use the 5 k=24 sweep panels: **519 loci, 373 decomposable, 90 P1,
> 68 catalogued** — same species, `91_` construction, 5 seed backgrounds.
> Report the sweep baseline per panel, never as a total: its P1 ranges 38/25/7/6/14
> across the five, so a sum is dominated by panel quality.


## 1. The dual-channel architecture

```
Layer 2 (LOCKED)  →  decomposition, 91.7% exact inserted length
   ├─ STRUCTURAL channel:  TSD + TIR
   │     asserts "arrived by insertion at some point"
   │     two genomes, no panel, microseconds
   └─ GENEALOGY channel:   Layer 3 (CLOSED)
         asserts "which state was ancestral at the ingroup MRCA"
         needs a panel, expensive
tier = which channels fire + whether they agree
Layer 5 annotation still runs last
```

### Why this is not a re-reading of a failed criterion

A pre-registered grid asked whether to **replace** Layer 3 with a structure-first
route (`C >= 0.40`, `K >= 0.95`, overlap `n >= 20`). It returned C = 32.5%,
K = 88.9%, n = 18 — **all three fail, and that proposition is rejected.** Layer 3
is not demoted.

The dual-channel proposal is a different claim and rests on a direct observation,
not on `K`: of **76 clean TSD calls, 55 were P0 and 3 were P2** — Layer 3 reached
only 18. The two channels are largely a **union, not a competition**. A
substitution rule could not express that outcome at any threshold, which is
itself why the pre-registration was misspecified (§2).

---

## 2. Two failed pre-registrations, archived verbatim

Kept so neither is later quoted as if valid. Both failed by **misspecification**,
not by a badly chosen threshold.

**k=48 rule (2026-08-31).** "Run k=48 if `Y(12→24)/Y(6→12) >= 0.50`."
Per-seed ratios: **0.50, 0.60, −1.27, 6.25**; median 0.55.
*Failure:* assumed a monotone dose-response, so the ratio could change sign. One
seed collapsed 38→5 targets, another jumped 4×. A median over sign-changing
quantities is meaningless. The decision was made on the structural cause instead
(the focal clade does not grow with the panel).

**TSD pivot grid (2026-09-01).** "Pivot to structure-first if `C >= 0.40` and
`K >= 0.95`, K untestable at n < 20." Results: 32.5%, 88.9%, n=18 — all fail.
*Failure 1:* assumed the channels were **alternatives**; they are mostly a union.
*Failure 2:* **K is not a validity measure.** TSD asserts "arrived by insertion at
some point"; Layer 3 asserts "which state was ancestral on this branch". Ancient
insertion followed by lineage-specific loss makes both true, so disagreement ≠
error — confirmed: both disagreements were already emitted as
`deletion_candidate` by Layer 4.

**Rule for future pre-registrations:** state the assumed *relationship* between
the quantities, not only the thresholds. Ask: can this metric change sign? Are
the options exclusive or additive? Does the comparator measure the same thing as
the candidate?

---

## 3. Tier definition

| tier | condition |
|---|---|
| **T_A** | both channels agree · junction ambiguity 0 · ≥2 carriers each side · junction ≥50 kb from nearest contig end |
| **T_B** | one channel only. If structural-only, **marked `per_locus_unverified`** |
| **T_C** | no structural signal (IS110-class) → genealogy channel only |
| **T_D** | channels disagree — the ancient-insertion-then-loss enrichment set, archived separately |

**Family copy number is demoted from gate to recorded field** (§5).
**Junction-to-contig-end distance enters the gate.**

### T_A is not scalable, and T_B is the shipping product

T_A requires both channels, and the genealogy channel requires a panel. There are
currently **18 overlapping loci, of which only 5 are clade-polymorphic**. T_A
cannot be produced at volume.

> **The bulk of output is T_B, and T_B is per-locus unverified.**

> **The structural channel has TWO SYSTEMATIC BLIND SPOTS, and the gaps are not
> random.**
>
> 1. **No-TSD families are invisible by construction.** IS110/IS1111 leave no TSD
>    and no TIR — measured 0/13 in this project's own Arm B output — and IS110 is
>    one of the two dominant length modes here (1279 bp). Those loci fall through
>    to Layer 3 unchanged.
> 2. **4–5 bp TSD families are invisible by POWER**, established 2026-09-02 from
>    the detector-matched null (§9c): chance placement already produces a 4 bp
>    repeat 47% of the time, so a 4 bp peak is undetectable at any realistic n
>    (>140 for 80% power) and 5 bp needs n≈100. **IS3, IS5, IS21 and IS1380 are
>    therefore not recoverable by the spike test.**
>
> Together these remove a *specific, named* slice of the element spectrum — the
> no-TSD families and the short-TSD families — not a random sample of it.
> **Therefore no statement about element spectrum, family composition, or
> "what is present in a genome" may be claimed as complete from the structural
> channel alone.** This sits alongside the T_B caveat above: both must appear
> wherever output is described, or a reader will read a census where none exists.

Any accuracy claim therefore rests on **distribution-level evidence** (§4) plus a
small T_A core as an anchor. This must be stated wherever output is described;
without it a reader will assume T_A's stringency covers everything shipped.

---

## 4. What actually supports the accuracy claim

**The spike test, not the concordance figure.** TSD length distribution: a
**9 bp peak** — IS1's 9 bp, recovered by a detector that has never seen an IS
library. This needs no external reference.

> **CORRECTED 2026-09-02 — the 4–5 bp shoulder is WITHDRAWN.** It was read as
> IS3's 3–4 bp. Against the detector-matched null (§9c), which is U-shaped and
> puts **47% of chance calls at 4 bp**, the E. coli histogram scores:
>
> | bin | 4 | 5 | 6 | 7 | 8 | **9** | 10 | 11 | 12 | 13 | 14 |
> |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
> | z | **−5.94** | +0.35 | −0.82 | +2.01 | +3.71 | **+10.18** | +3.01 | −1.07 | +0.11 | −1.78 | −1.40 |
>
> **9 bp is the only bin that clears the critical value (~3.9).** The 5 bp bin is
> **indistinguishable from chance** (z = +0.35) and the 4 bp bin is significantly
> **DEPLETED** (z = −5.94) — the detector calls *fewer* 4 bp repeats than random
> placement produces. The earlier reading compared the histogram against a flat
> expectation, which does not describe this detector.
>
> **The claim is one peak, IS1's 9 bp. Nothing about IS3.** 8 bp (z = +3.71,
> IS6/IS66/ISL3/IS256 territory) is suggestive but does not clear the threshold
> and is not claimed either. A detector reporting noise does NOT yield a flat
> distribution here — it yields the U-shape — so "not flat" was never the right
> test.

It is **distributional evidence that the detector is not noise. It does not
certify any individual call.** Of 76 clean TSD calls, only 18 have an independent
reference.

**Conservative concordance: 85.0% on 85 calls.** The `--max-tsd 15` guard has a
genuine prior basis (real TSDs are 2–15 bp) but its effect on concordance
(85.0 → 88.9%) was found after the fact, so 88.9% is not clean out-of-sample.

**Chance matching is real:** without the ≥4 bp floor, C would be 38.5%, inflated
by 10 loci at 2 bp and 4 at 3 bp. Report both; decide on 76.

---

## 5. Collapse risk — RETRACTED, and the covariates that replace it

The earlier claim "18 of 20 ready targets are HIGH repeat-collapse risk" is
**withdrawn**. The metric (same-size family copy number in a carrier) measures
whether the family is high-copy, which for IS1/IS110 is true *by definition* —
every genuine IS polymorphism scores HIGH. Specificity ≈ 0.

**The mechanism was impossible.** Collapsing `TSD–insert–TSD` into a clean empty
allele would require the ~9 bp direct repeat to act as a collapsible unit, but
9 bp is far below any assembler's k (31–127): k-mers spanning the TSD are unique
and the de Bruijn graph never forms that loop. What high-copy IS actually causes
is a **contig break** at the element — and a break would put Layer 2's anchors on
different contigs, where `place_locus()` fails or reports `ambiguous`.

**Locus-specific covariates, now standard fields, measured on all 20 targets:**

| covariate | source | result |
|---|---|---|
| `junction_to_contig_end` | carrier `.fai` + `junction_L_abs` | **min 45,508 bp; median ~1.2 Mb; 0 within 5 kb** |
| `carrier_n_contig` | PDG `asm_stats_n_contig` | **1 or 3** |
| `carrier_contig_n50` | PDG `asm_stats_contig_n50` | complete-level |
| `assembly_method` | PDG `assembly_method` | **PacBio SMRT, HGAP2/3, CLC — all long-read** |
| `family_copy_number` | Arm A `-f 0` | recorded only, **not a gate** |

**Output:** `catalog/collapse_covariates.tsv`.
**Criterion:** flag `collapse_risk = elevated` only if
`junction_to_contig_end < 5000` **or** (`assembly_method` is short-read-only
**and** `carrier_n_contig > 50`). On current data that flags **0 of 20**.

---

## 6. Genealogy channel — new job: manufacture verification overlap

Because the structural channel needs no panel, **§6 of PROJECT_SUMMARY (should the
focal clade grow with the panel?) is no longer blocking.** It affects only how
fast the verification set grows.

> The genealogy channel's purpose is now **verification overlap, not yield.**

Overlap today is 18 loci (5 clade-polymorphic). The cheapest way to strengthen the
accuracy claim is **not** a larger k sweep producing more P1 in general, but
panels built *deliberately around the 76 TSD-called loci*.

**Target: overlap 18 → 60–70.** Same compute, buying verification instead of
duplicate discovery.

**Output:** `verification/overlap_panels.tsv`, `verification/concordance.tsv`.
**Criterion:** report concordance **stratified by clade-polymorphic vs
clade-constant** — an invariance test run on a set that is 80% clade-constant
passes almost automatically and certifies nothing.

---

## 7. The 93 undecomposable loci — PENDING, and probably the most interesting set

> ### RETRACTION 4 (2026-09-01) — the independent-acquisition result is withdrawn
>
> **Everything below about "10 of 11", "5 genome pairs" and "two independent
> genomic contexts" is WRONG at the inclusion step, and the corrections that
> shrank it 23 → 11 → 7 → 5 → 2 were all downstream of a broken entry
> condition.** The finding has been re-derived correctly; see
> "Corrected scan" below.
>
> **What was wrong.** The test compared *the two longest alleles* at a locus and
> called their differing middles "two elements at one target site". It never
> checked that **both** alleles were actually long. At **8 of the 11** loci there
> is only ONE long allele — the second-longest is the same length as the site's
> shortest:
>
> | locus | allele lengths | excess over shortest |
> |---|---|---|
> | cl0000.B000014 | 597, 597, 597, **1152** | 0 and 555 |
> | cl0000.B000016 | 504, 504, **2564** | 0 and 2060 |
> | cl0001.B000023 | 644, **2439** | 0 and 1795 |
> | cl0001.B000109 | 578, **1295** | 0 and 717 |
> | cl0001.B000121 | 525, 527, 527, **5524** | 2 and 4999 |
> | cl0001.B000126 | 596, 596, **1152** | 0 and 556 |
>
> So the "middle-vs-middle" comparison was **the inserted element against the
> un-inserted target-site residual** — an ordinary single insertion, not two
> elements.
>
> **The test could not have come out negative.** Scored on the 82 plain
> single-insertion loci in the same data (exactly one long allele), the old
> comparison gives a median identity of **0.908**, and **23 of 82 (28%) land
> inside the 0.73–0.90 "diverged homologue" band** — including
> `cl0001.B000073`, `cl0000.B000014` and `cl0000.B000016` themselves. The large
> flank-minus-middle gap follows automatically: flanks match at ~99% because they
> are the same locus, and the element does not match the empty target because it
> is an insertion. **Every ordinary insertion in the dataset produces this
> signature.** This is the session's failure shape a fifth time — an inclusion
> rule that admits the confirming case — and it is the one that got furthest
> before being caught.
>
> **Survivors of the old set: 3 of 11** — `cl0000.B000048`, `cl0000.B000052`,
> `cl0001.B000061` — which do have two genuinely long alleles.
>
> **Fix, now enforced in code** (`scripts/96_shared_target_scan.py`): both of the
> two longest alleles must exceed the site's SHORTEST allele by at least
> `--min-middle`. A 2-allele locus is excluded automatically, since its
> second-longest IS the shortest.

### Corrected scan + both nulls — 2026-09-01

`scripts/96_shared_target_scan.py` (entry condition, measured floor, NA-not-zero,
genome-derived flanks, uniform coverage cutoff, component-level effective n) and
`scripts/97_shared_target_nulls.py` (coincidence permutation, divergence shape).

| | original 234 | sweep 524 (5 k=24 panels) |
|---|---:|---:|
| two genuinely long alleles | 67 | 127 |
| rejected as short-vs-long | 143 | — |
| measured noise floor | 0.473 | — |
| diverged homologue | 11 | 13 |
| **at DIFFERENT positions** (null 1) | **29 of 67** | **61 of 127** |
| nested (null 2) | 51 of 67 | 111 of 127 |
| diverged + well-covered + same-site + not nested | **3** | **6** |

**Null 1 fires hard.** 29 and 61 loci are blocky, well-covered and at genuinely
different positions — 45, 66, 114 bp apart. Every one would have counted under
the old framing.

**TEST 1 HAD NO DISCRIMINATING POWER — do not cite its p-value as support.**
The permutation null was "random positions", so rejecting it selects against
chance but not *among the non-chance explanations*. Common descent produces
identical offsets **deterministically**: one ancestral insertion inherited by two
genomes sits at the same base with probability 1. So p<0.0001 is exactly what
vertical inheritance predicts, and Test 1's result is fully explained by the
mechanism Test 2 found. It never bore on the question. The framing error was
symmetric with the original one — landing outside a null distribution was treated
as evidence for the hypothesis, when the null only excluded one alternative.
Recorded numbers, for completeness only: survivors 3 of 3 at delta==0 against
NULL B mean 0.05, p<0.0001 (sweep p=0.0002), 20,000 draws, NULL B resampling the
pooled observed relative-offset distribution.

**TEST 2 — divergence shape is what settles it, and the global identity was
hiding a mosaic.** Windowed identity along the alignment, 100 bp windows. The
MEDIAN window decides; blockiness alone conflates two opposite readings:

| locus | global | median window | reading |
|---|---:|---:|---|
| cl0000.B000052 | 0.866 | **0.995** | same element + divergent block — ONE arrival |
| cl0000.B000048 | 0.822 | 0.750 | intermediate |
| cl0001.B000043 | 0.836 | 0.885 | intermediate |
| GCF_009618015.1k24.B000016 | 0.872 | **0.960** | same element + block — ONE arrival |
| GCF_009618015.1k24.B000147 | 0.828 | **0.955** | same element + block — ONE arrival |
| GCF_041319995.1k24.B000051 | 0.855 | **1.000** | same element + block — ONE arrival |
| GCF_009618015.1k24.B000146 | 0.636 | **0.520** | two different elements |
| GCF_041319995.1k24.B000076 | 0.541 | **0.480** | two different elements |
| GCF_041319995.1k24.B000151 | 0.699 | **0.620** | two different elements |

`cl0000.B000052` reads 0.866 globally but its window profile is
`1.000 x16, then 0.65 0.60 0.57 0.57 0.58 0.51 0.46 0.52 0.64, then 0.97...` —
one element, essentially identical, with a discrete recombined block. **Identical
offset there needs no coincidence at all: one ancestral insertion inherited by
both genomes.** Four of the nine survivors are this.

**Result: 0 of 234 original loci, and 3 of 524 sweep loci, are consistent with
two different elements at the identical base.** And the flank control thins those
three further:

| locus | mid id | flank id | gap | verdict |
|---|---:|---:|---:|---|
| GCF_009618015.1k24.B000146 | 0.636 | 0.965 | 32.9 | **clean** |
| GCF_041319995.1k24.B000076 | 0.541 | 0.862 | 32.1 | marginal — flanks only 86% |
| GCF_041319995.1k24.B000151 | 0.699 | **0.728** | 2.9 | **fails** — divergent region, no contrast |

**HONEST STATE: 1 clean locus and 1 marginal one, out of 758 E. coli loci across
7 genome backgrounds.** That is indistinguishable from nothing.

**SCOPE OF THE NEGATIVE — precise in both directions.** What was tested is
**within-cluster divergence in E. coli**: genomes at >=99% ANI, where very little
time has elapsed for a second element to arrive at an already-occupied site. The
statement is therefore **"not supported at within-cluster divergence in E. coli"**
— *not* "does not occur". Equally, this is **not** a reason to reopen: the two
places the phenomenon would be more visible are wider ANI, already rejected on
the ~95% anchor-placement cliff ([[outgroup-distance-tradeoff]]), and the parked
oversized bubbles (stage 32, still unwritten). §7 is CLOSED, not paused.

### A retracted negative

An earlier check asked whether a third, shorter allele decomposed cleanly against
both long alleles, and reported **0/93** — treated at the time as closing the
question. **That test was circular.** These loci entered the bucket *because*
decomposition against the shortest allele failed; re-running "decompose against
the shortest" restates the inclusion condition, so the headline "80/93, not even
one decomposes" was necessarily true. The design is **confirm-only**: a positive
would have been sufficient, a negative is uninformative, because the shortest
allele need not be an empty allele at all.

### The test that can refute — and what it found

Compare the two **long** alleles against each other and ask about the *shape* of
divergence. Coverage 1.000 at identity 0.854 is compatible with two opposite
structures: divergence spread uniformly (two drifted homologues — hypothesis
wrong) or concentrated in a block with near-identical flanks (two different
elements at one site — hypothesis right).

**This design could have refuted: 26 loci returned uniform divergence**, the
outcome that kills the hypothesis. It therefore meets the §11 bar, unlike the
three-allele check it replaced.

| divergence shape | n |
|---|---:|
| flanks match, middle also similar | 38 |
| **uniform divergence — hypothesis wrong** | **26** |
| BLOCKY: flanks match, middle divergent | **29** |
| other / mixed | — |

*(A first pass reported 23 BLOCKY. The windowed identity returned 0.000 when a
window had **no aligned columns**, which is not 0% identity. Fixed; 23 → 29. The
apparent "middle identity 0.000", read at the time as the strongest evidence, was
partly this artefact.)*

### Three follow-up tests, and what survives

**1. Orientation — clean.** IS elements insert in either direction, and this
project has been burned once by a forward-only comparison. Every divergent-middle
pair was scored against `revcomp` as well: **0 cases of "same element, reverse
complemented"**. The revcomp scores also give an empirical noise floor
(median identity **0.488**), which is what an unrelated pair actually scores —
not 0.000.

**2. Reclassification against that floor.** Of 26 loci with two substantial
middles:

| forward comparison | n |
|---|---:|
| diverged homologue (0.73–0.87 identity, high coverage) | 11 |
| same element, lengths differ (≥0.90 identity) | 5 |
| **genuinely unrelated (≈ noise floor)** | **6 strict, 8 loose** |
| ambiguous (high identity, low coverage) | 2 |

**The "diverged homologue" class may be the STRONGEST evidence, not a weaker
one — this reclassification probably demoted it wrongly.** These genomes sit
inside an ANI cluster at ≥99% identity, with focal-clade floors up to 99.98. If
an element had been inserted once in their common ancestor and inherited
vertically, the two copies at that locus would be ~99% identical like the rest of
the genome. They are **73–87% identical**. Vertical inheritance in place cannot
produce that gap; two related but distinct elements arriving independently at the
same site can. Homology between the middles says they are the same family; the
identity gap says they arrived separately — which is exactly the shared-target
claim, with more information than an unrelated pair carries.

The decisive comparison is per-locus middle identity against a control that
describes the same genome pair. **RUN 2026-09-01, three checks — the finding
survives, at a smaller effective n and with one control replaced.**

**Effective n is 7 genome pairs, not 11 loci.** The 11 loci are not 11
independent observations: one pair, `GCA_001559675.1 / GCA_002741575.1`, supplies
**5** of them, and `GCA_000258145.1` participates in 3 more. A single divergent
genome pair produces elevated middle divergence at *every* site it shares, so
those 5 loci are largely one observation repeated. The headline is **"confirmed
in 7 genome pairs"**, not 10 of 11.

**The within-locus control had to be rebuilt.** Flank identity taken from the
Layer-2 allele is **circular** — the "flanks" *are* the longest common prefix and
suffix, so decomposition defines them as exactly matching and the measurement
returns 100.00 on every locus, always. The valid control extracts **2 kb from
each carrier genome outside the allele boundary entirely** and aligns those. That
is independent of the decomposition: same locus, same pair, same estimator.

| locus | mid id% | mid cov% | ext-flank id% | flank cov% | gap |
|---|---:|---:|---:|---:|---:|
| cl0000.B000014 | 85.60 | 97.7 | 98.92 | 100.0 | 13.33 |
| cl0000.B000016 | 74.41 | 96.4 | 99.27 | 100.0 | 24.86 |
| cl0000.B000048 | 82.23 | 98.4 | 98.67 | 99.7 | 16.44 |
| cl0000.B000052 | 86.58 | 98.3 | 98.88 | 100.0 | 12.29 |
| cl0001.B000005 | 73.37 | 67.2 | 98.90 | 100.0 | 25.53 |
| cl0001.B000023 | 89.63 | **55.5** | 98.20 | 100.0 | 8.57 |
| cl0001.B000061 | 78.20 | **66.2** | 97.98 | 98.9 | 19.77 |
| cl0001.B000073 | 86.52 | 98.2 | 99.57 | 100.0 | 13.06 |
| cl0001.B000109 | 76.88 | 78.7 | 96.77 | 99.1 | 19.89 |
| cl0001.B000121 | 73.52 | 79.8 | **86.25** | 96.6 | 12.73 |
| cl0001.B000126 | 83.13 | 87.1 | 99.55 | 100.0 | 16.42 |

External-flank-minus-middle gap: median **16.4** points, range 8.6–25.5, n=11.
The gap is not an artefact of the wrong control — it holds against local flanking
sequence, not only against core-genome ANI, and the two estimators agree
(genomic ANI 99.02–99.94; external flanks 96.8–99.6 on ten of eleven).

**What the external flanks bought that was not claimed: the recombination reading
is now weakened.** The standing caveat was that identity alone cannot exclude
recombinational replacement of one element by another. But the flanks are at
~99% only **2 kb out** from middles at 74%. Any imported tract would therefore
have to be confined to roughly element length; homologous recombination tracts in
*E. coli* are typically longer than that, so a generic recombination import does
not fit the observed geometry. And the version that does survive — a short,
precisely bounded swap at the element itself — is element-mediated at a shared
site. **Both surviving readings still have two elements using one target site**,
which narrows the four-way interpretation set the previous round left open.


**Coverage, with a stated cutoff.** Middle alignment coverage median **87.1%**,
minimum 55.5%. **Cutoff: a locus is well-covered at coverage >= 75%.** Three fall
below it and are set aside — `cl0001.B000023` (55.5%), `cl0001.B000061` (66.2%),
`cl0001.B000005` (67.2%): high identity over a partial alignment, the previous
round's "ambiguous" category. That leaves **8 well-covered loci**.

One well-covered locus is weakened differently: `cl0001.B000121` has coverage
79.8% but its own external flanks are only **86.25%** identical, so it sits in a
genuinely divergent region and the within-locus contrast is small.

**Validity limit on the flank control:** flank matching was an *inclusion
criterion* for BLOCKY, so the control is sound for the loci in hand but cannot
estimate how often such loci occur.

The **5 "same element, ≥0.90"** are genuinely weaker — a truncated copy is
consistent with one arrival plus a deletion.

**Headline number — pairs within the well-covered subset.** Restricting to the
8 well-covered loci collapses the pair count from 7 to **5**:

| genome pair | well-covered loci |
|---|---:|
| GCA_001559675.1 / GCA_002741575.1 | **4** |
| GCA_000258145.1 / GCA_001886935.1 | 1 |
| GCA_000258145.1 / GCA_002057355.1 | 1 |
| GCA_000258145.1 / GCA_009664205.1 | 1 |
| GCA_001886935.1 / GCA_003288415.1 | 1 |

**The precise statement is two independent genomic contexts, not five pairs.**
The pair graph has exactly **2 connected components**, and the 8 loci split 4+4
across them:

| component | genomes | pairs | loci |
|---|---:|---:|---:|
| GCA_000258145.1, GCA_001886935.1, GCA_002057355.1, GCA_003288415.1, GCA_009664205.1 | 5 | 4 | 4 |
| GCA_001559675.1, GCA_002741575.1 | 2 | 1 | 4 |

This understates nothing and overstates nothing: the raw count of 5 pairs hides
that four of them share genomes (`GCA_000258145.1` in 3, `GCA_001886935.1` in 2),
while "5 non-independent pairs" hides that the evidence is **not** concentrated
in one strain background. **Two contexts is the honest ceiling.**

Defensible count: **6–8 unrelated pairs, plus same-family independent acquisition
in 8 well-covered loci across 5 genome pairs drawn from 7 genomes.** This still
reverses the earlier demotion of the diverged-homologue class — the mechanism
argument is threshold-free and does not depend on the count — but state it as a
5-pair result, never as 10 of 11.

**3. Per-middle TSD — negative, but the window was probably wrong.** If each
middle arrived by insertion, each should carry its own 4–15 bp direct repeat.
Result: **0 of 8 carry a TSD on either middle** — but the test searched strictly
*inside* the middle, and the middle boundary is where the two long alleles stop
agreeing, which is **not** the element boundary. If both elements sit at the same
site, the shared target sequence is in the shared *prefix*, so a TSD would fall
on or outside the boundary, not inside the searched window. The corrected test
searches ±20 bp around each boundary (queue item 12). Even corrected, a negative
would be weak: IS110 leaves no TSD and is one of the two dominant families.

### Honest status

**What this finding is: a lead, not a result.** 11 loci, 8 well-covered, 5
non-independent genome pairs from 7 genomes, 2 clusters, one species, out of 234
Layer-2 loci. It is the closest the pipeline has come to the founding question
and it is **annotation-free**, which is the part that matters. Whether it
generalises is what the §9 per-species census tests — not anything further that
can be done to this set.

The **observation** is real and survives every check: at 6–8 loci, two alleles
share near-identical flanks and carry mutually unrelated middles of 300–4,500 bp,
in the same orientation.

The **interpretation is not established.** There is no independent structural
evidence that either middle arrived by insertion. At least four readings fit:

- two different elements at one target site (the original hypothesis)
- a second element inserted into the first (nested — would show one TSD, not two)
- ~~one element replaced by another via recombination~~ — **weakened 2026-09-01**:
  flanks at ~99% just 2 kb out bound any imported tract to about element length,
  shorter than typical *E. coli* recombination tracts; the surviving short-swap
  version is element-mediated at a shared site anyway
- one element plus a large independent indel in the same window

Two of the four readings still place two elements at one target site, so the
interpretation set is narrower than the previous round recorded — but not closed:
the nested and independent-indel readings remain live.

**§7 records the observation, not the interpretation.** The multi-allelic module
stays **pending**: 6–8 loci is a thin basis, and the TSD negative means the
cheapest confirmation route is already exhausted. Separating the four readings
needs annotation (Layer 5, legitimately — these are candidates, not discovery) or
read data.

**Also open:** the 3 loci where exactly one long allele decomposes against the
shortest — the shape of "one element recognised, the other diverged". They were
wrongly counted inside the 80 as failures.

**Output:** `catalog/shared_target_candidates.tsv` — locus, both allele ids and
lengths, shared prefix/suffix, both middle lengths, forward and revcomp
identity/coverage, per-middle TSD, per-pair genomic ANI, and the four-way
classification.

### A note on the noise floor, and its one circularity

Using revcomp scores as an empirical null is worth carrying forward as standard
practice: revcomp preserves base composition exactly, so it is a
composition-matched null rather than an arbitrary cutoff. It is what exposed that
the 0.90 identity cutoff was wrong.

**The circularity:** "0 revcomp cases" and "revcomp defines the floor" are read
off the same comparisons. A genuine reverse-complement pair would both inflate
the floor and be a missed orientation call. Max 0.683 against a median of 0.488
is a wide enough gap to check — queue item 13.

## 8. Event-level deduplication

Production will run many overlapping panels, so the same event will be
rediscovered repeatedly.

**Key:** `md5(inserted_sequence)` + `md5(pre_insertion_context_200)`.
Both are already emitted by Layer 4. Insert identity alone is insufficient — the
same element at two different sites is two events; the same site with two
elements is two events.

**Output:** `catalog/events_dedup.tsv` with `event_key`, `n_panels_observed`,
`panels`, `tier_best`.

---

## 9. Per-species census and parameter calibration

**Do not hard-code thresholds tuned on E. coli.** Per species, derive:

| parameter | rule |
|---|---|
| ingroup ANI floor | quantile of the observed within-species distribution, **not** a constant 99.0 |
| outgroup band | anchored to the measured anchor-placement cliff (E. coli: usable to ~95–96%, collapses below 95%) |
| `--min-af` | **not** a proxy for anchor placement — measured 38 genomes below AF 80 still placed 80% of anchors |
| 5 kb cap | keep, with the `armB_large_events.tsv` parking mechanism intact |

**CENSUS RUN 2026-09-01** (`scripts/95_species_census.py`, 13 s over 3,282,482
GenBank rows, 107,887 species_taxids). Ranked by **usable = complete +
chromosome**, because Arm B compares whole assemblies and a contig-level pool
contributes fragmentation, not alleles:

| species | usable | total | usable% |
|---|---:|---:|---:|
| Escherichia coli | 8,961 | 423,326 | 2.12% |
| **Klebsiella pneumoniae** | **5,916** | 132,633 | 4.46% |
| Staphylococcus aureus | 4,209 | 130,846 | 3.22% |
| Salmonella enterica | 3,545 | 626,212 | 0.57% |
| Pseudomonas aeruginosa | 2,085 | 57,269 | 3.64% |

**194 species have >= 50 usable genomes; 99,274 usable genomes exist in total.**
So the substrate for a multi-species study is real, and E. coli is 9% of it.

**Ranking by total would have picked the wrong pilot.** Salmonella enterica is
the largest species in GenBank at 626,212 assemblies but only **0.57%** usable —
3,545, fewer than half of E. coli's. Genome count and discovery substrate are
different quantities, and only the second one matters here.

**Species selection — resolved by the pre-registered rule:** top of the census
by usable count, `--species` to override. E. coli is #1 and already done, so
species 2 is **Klebsiella pneumoniae** (5,916 usable). Download submitted
2026-09-01, job 25457001, 16-way array on cf1 ->
`genomes/kpneumoniae/`. `Staphylococcus aureus` is the stronger generalisation
test (Gram-positive, outside Enterobacteriaceae) and is the intended species 3;
K. pneumoniae is the near test that also checks whether the E. coli-tuned
parameters survive a congeneric move at all.

**Output:** `census/genbank_species_census.tsv` (all 107,887 rows),
`census/kpneumoniae_manifest.tsv`.

### 9a. Dedup — NOT needed for an NCBI-sourced pool (recorded so it isn't inherited)

K. pneumoniae: **0 accession-core duplicates, 5,916 GCA / 0 GCF, 3 duplicate
BioSamples.** E. coli's 13,027 -> 7,723 (40.7%) loss was a **local-pool artifact**
— that pool mixed sources and carried GCA/GCF pairs. `assembly_summary_genbank`
lists only GCA, so a pool built from it needs no accession dedup. Do not inherit
a phantom cleanup step. The key is the full versioned accession matching the
filename; panel builder, Arm B and Layer 2 must all use it (item 11's ANI lookup
failed precisely because dedup preferred GCF while carriers were named GCA).

### 9b. ANI threshold — quantile is a DEFAULT, the gap/floor test is the calibration

Measured, full pools, same estimator:

| | E. coli | K. pneumoniae |
|---|---:|---:|
| median reported ANI | 97.35 | **99.02** |
| quantile of ANI 99.0 | 0.9266 | **~0.24** |
| neighbours at 99.0 / AF85 | tens | median **2,383** |

**ANI 99.0 sits on K. pneumoniae's mode** — it would retain ~76% of reported
pairs instead of the top 7%, lump clonal complexes, and balance the threshold on
a peak where 0.05 ANI moves millions of pairs. Transferring E. coli's *quantile*
(0.9266) gives **ANI 99.59**.

**But the quantile is not a calibration, and must not be recorded as one.**
skani reports only pairs above ~80% ANI, so these are quantiles of a **truncated,
pool-composition-dependent set**: q=0.9266 is not the same object in two species.
K. pneumoniae's reported set is dominated by within-species pairs at high ANI;
E. coli's pool had a different composition. Quantile beats constant as a starting
value; it is not evidence the threshold is right.

**THE PROPERTY THAT MUST TRANSFER IS THE ONE LAYER 3 DEPENDS ON:** a focal clade
with an internal ANI floor and a separation gap from its outgroups. §6 established
that gap and floor predict productivity better than any panel-size parameter, so
they are the target. E. coli reference, 5 k=24 panels:

| | median | range |
|---|---:|---|
| focal-clade ingroup ANI floor | **99.80** | 99.62–100.00 |
| ANI gap to closest outgroup | **0.98** | 0.65–1.19 |
| focal clade size (of 24) | 3 | 3–5 |

**Decision rule, pre-registered:** run panels at ANI 99.59; if the resulting
median floor and gap are comparable to the above, keep it. **If not, the
threshold moves on the gap/floor criterion, not on quantile match.** Record the
value, its quantile, AND the achieved floor/gap in `census/<species>_params.tsv`.

**Cost warning:** at median 2,383 neighbours, the §4.3 greedy clique growth is far
more expensive than in E. coli (median usable neighbourhood fraction 0.104).
**Time ONE panel before launching the full set.**

### 9c. TSD spike test — PRE-REGISTERED, written before any K. pneumoniae data

The informative outcome is **where the peak sits**, not whether one exists.
E. coli peaked at **9 bp** — IS1, recovered by a detector that never saw an IS
library, and matched to a family only *after the fact*. (The 4–5 bp "IS3
shoulder" is withdrawn: see §4. It is chance under this null.) Post-hoc matching to a family whose TSD length is already known is weak
evidence, so the admissible outcomes are named below **before looking**.

#### The null is NOT uniform — it is U-shaped, and measured

A uniform null over 4–15 bp is too generous: chance repeat length is k-dependent,
so a detector with no family signal still produces a sloped histogram. The null
used here is **detector-matched and composition-matched**: the same `find_tsd`,
on the same real alleles, at **random interval placements** instead of the true
junction pair (720 draws in-window, true junction excluded by +/-25 bp).

| bin | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| null | **.471** | .131 | .062 | .036 | .051 | .042 | .046 | .036 | .037 | .040 | .047 |

**66% of chance calls land at 4–6 bp**, and the raw search also piles up at the
15 bp ceiling (67% of all random hits before censoring — that bin is excluded as
censored, and is what `--max-tsd` exists to catch). Consequences:

- a **9 bp peak is MORE significant** than the uniform null implied (null 4.2%,
  not 9.1%) — E. coli's 27.6% modal bin is a strong result;
- a **4–6 bp peak is far LESS significant** than it looks.

#### Power is position-dependent, and 4 bp is a BLIND SPOT

Statistic: **max standardized excess over the null proportions** (not max raw
bin), critical value ~3.9 = 99th percentile of the null. n for 80% power against
an E. coli-shaped peak (27.6% modal), 4,000 trials:

| peak at | 4 bp | 5 bp | 6 bp | 8 bp | 9 bp | 12 bp |
|---|---:|---:|---:|---:|---:|---:|
| n needed | **>140** | 100 | 40 | 30 | 30 | 25 |

> **A 4 bp peak is not detectable by this instrument at any realistic n** — the
> null already puts 47% of chance calls there. **IS3, IS5, IS21 and IS1380
> (TSD 3–5 bp) are therefore invisible to the spike test**, a blind spot
> alongside the known IS110/no-TSD one. This was not known before the null was
> measured; the uniform-null answer of "n >= 50" concealed it in both directions.

#### PRE-REGISTERED OUTCOMES

> **Floor: n >= 30 clean TSD calls.** Below 30, report **UNDERPOWERED** — never
> "the detector does not transfer". Report achieved power for the observed peak
> position alongside the histogram, since the floor covers 6–14 bp only.
>
> **O1 — peak at 9 bp** (same as E. coli): consistent with transfer, but the
> **weakest** positive. Enterobacteriaceae share IS1-like content, so this may
> reflect shared element content rather than independent per-species recovery.
>
> **O2 — peak at a different canonical length: {6, 8, 11, 12, 13}** (IS1182 6;
> IS6/IS66/ISL3/IS256 8; IS4 11–13). **The strongest outcome**: the detector
> transfers AND resolves a real species difference with no annotation.
>
> **O3 — peak at a non-canonical detectable length: {7, 10, 14}.** Interesting —
> either a novel target-site preference or an artifact. A follow-up, never a
> claim on its own.
>
> **O4 — peak at 4–5 bp:** NOT CALLABLE. Declared a blind spot in advance; do
> not report it as a finding in either direction.
>
> **O5 — flat with n >= 30:** a real negative. Report it as one, and apply no
> tier to the species until it is explained.

> **THE NULL DOES NOT TRANSFER EITHER — recompute it per species.** The U-shape
> above, its critical value (~3.9) and the n >= 30 floor were all derived from
> **E. coli's allele length distribution**, because chance repeat frequency
> depends on interval length and composition. K. pneumoniae's alleles differ, so
> its null shape, critical value and n floor must be **recomputed from its own
> Layer-2 output** before the histogram is read. This is the same class of error
> as carrying ANI 99.0 across species (§9b): a value measured on one species is
> not a constant.

Canonical list fixed in advance (modal TSD): IS1 9; IS3 3–4; IS4 11–13; IS5 4;
IS6 8; IS21 4; IS30 2–3; IS66 8; IS110 none; IS200/IS605 none; IS256 8–9;
IS630 2; IS982 3–4; IS1182 6; ISL3 8; Tn3 5; ISAs1 9; IS1380 4–5; IS91 none.
Of the 9 detectable bins, 6 are canonical and 3 are not, so O2/O3 is a real
distinction rather than a formality.

---

## 10. Read gold — two pre-written branches

The 70 read-validated junctions are **not on scratch** (searched: nanopore dirs
are IS110 length/consensus work; `empty_junctions.fa` holds 167,255
*assembly-derived* junctions; no VCFs anywhere; T1/T2/T3 tiering exists only in
the synthetic bench). Shell history shows rsync from `igi.biotite.berkeley.edu`.

Priority is **one notch lower than previously stated**, because the retraction in
§5 removed the claim that read data was the only way to resolve 18/20 targets. It
blocks precision claims; it does not block deployment.

**Branch A — obtainable.** Run tier calibration against it. Language becomes
*"tier calibrated"*; per-locus nucleotide-precision claims permitted within the
measured tolerance.

**Branch B — not obtainable.** All precision language stops at
distribution-level evidence. Permitted: *"TSD length distribution matches the
canonical bacterial modes"*, *"85.0% concordant with the genealogy channel on 85
overlapping calls"*. **Not permitted:** any per-locus nucleotide-precision claim,
and any restatement of the coarse concordance figures as junction precision.

---

## 11. Work queue

| # | task | cost | blocking |
|---|---|---|---|
| 1 | collapse covariates → standard fields | done, needs wiring | no |
| 2 | tier assignment implementing §3 | small | no |
| 3 | verification panels, overlap 18 → 60+ | moderate | no |
| 4 | event dedup key | small | production |
| 5 | per-species census + calibration | moderate | multi-species |
| 6 | read gold branch decision | **user** | precision claims |
| 7 | Arm B `cl0000` 0-bubble bug | small | upstream of everything |
| 8 | stage 32 nested sub-5 kb search | moderate | no |
| 9 | shared-target candidates: long-vs-long profile for all 93 (§7) | small | no |
| 10 | the 3 "one allele decomposes" loci, examined separately (§7) | small | no |
| 11 | ~~middle-vs-middle identity vs per-pair genomic ANI~~ **DONE 2026-09-01: 10/11 independent acquisition** | done | — |
| 12 | TSD re-search in a ±20 bp window around middle boundaries, on the reclassification survivors (§7) | small | no |
| 13 | coverage on the max-0.683 revcomp pair — if high-coverage it is not noise, and both the floor and the orientation call are affected (§7) | trivial | no |

**Not on the queue, rejected on evidence:** k=48 sweep; Pool O (sister species
sit below the 95% anchor cliff); more E. coli for density (94.1% already have ≥5
partners); loosening `--allele-margin` or the 0.95 identity gate.

Every item on that list was ruled out by a test that *could* have gone the other
way. The multi-allelic module was previously listed here and has been **removed**
— it was tested by a confirm-only design, and §7 now shows 23/93 loci supporting
it. Mixing an underpowered result into this list would force the whole list to be
re-audited later.

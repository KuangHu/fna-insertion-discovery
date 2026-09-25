# Simplification — the pipeline is discovery plus reconstruction

**Drafted 2026-09-02.** Supersedes `PIPELINE_RECTIFICATION.md` §2 (the S×D tier
grid), §4 (the CT core test) and §5 (the genealogy channel's verification role),
and the whole of `NEXT_STEPS.md` §9c (the spike test pre-registration).

**The goal, restated plainly:** for each locus where two closely related genomes
differ by an inserted segment under 5 kb, emit the empty target sequence and the
insertion offset, reasonably accurately, at volume.

Nothing else is in scope. Ancestry is not in scope. Family assignment is not in
scope. Nucleotide-certified junctions are not in scope.

---

## 1. What the pipeline is

Two algorithms and one join.

```
#1  DISCOVERY      flank[target+insert]flank  vs  flank[target]flank
                   minigraph bubbles over a set of comparable genomes
                   insert < 5 kb; oversized parked, never dropped

#2  RECONSTRUCTION cut a 500 bp anchor each side, map both into every genome,
                   extract the interval BETWEEN them verbatim -- including at 0 bp
                   decompose long against short -> offset + signed overlap
                   -> flank[target]flank  +  insertion_point_offset

M   MOBILITY JOIN  is the insert present at >=3 dispersed loci? (Arm A)
                   -> a field on the record
```

`pre[:offset] + insert + pre[offset:]` regenerates the observed derived allele.
That reconstruction, and the verbatim target between the anchors, is the product.

---

## 2. Removed

| module | disposition | reason |
|---|---|---|
| `80/81/83/85_*` Layer 3 chain | removed from the pipeline | answers "which allele is ancestral at the ingroup MRCA" — a phylogenetic claim, not a target site |
| `84_polarity_ingroup_refiner.py` | deleted | the focal clade exists only to serve Layer 3 |
| `89_layer3_regression.py` | retired with Layer 3 | gates a module no longer on the path |
| `92_structural_polarity.py` | deleted | `junction_ambiguity_bp` reproduces its calls exactly — 76/76, lengths 100% identical, conflicts 49 → 0 |
| `31_armB_refine.py` | already off the path | `svim-asm` and `nucdiff` absent; the 200 bp pad means Layer 2 relocates the junction from sequence anyway. Verified 100/100 on E. coli |
| `30b_armB_pangraph.py` | deleted | the minigraph path is validated; a second engine buys nothing |
| `40_merge_catalog.py` | superseded by `99_mobility_join.py` | it unions the two arms; what is needed is a join |
| `10_skani_cluster.py` | demoted to historical reproduction | single-linkage clusters are not clades (cl0000 spans 98.92–99.94 while its outgroups reach 98.97) |
| `86_pre_insertion_target_catalog.py` | rewritten, not deleted | stops being a polarity gate, becomes the catalogue writer |

### What Layer 3 was actually protecting — measured before deleting it

The Layer 4 catalogue holds 34 rows: 20 insertion + 14 `deletion_candidate`.
Without Layer 3 those 14 would be emitted as targets. Read against S1:

| | n | |
|---|---:|---|
| **S1 already stops** | **8** | `target_bases_lost` 84–592 bp — imprecise excision, exactly what S1 catches |
| S1 admits | 6 | `lost == 0`, overlap 0–7 — Layer 3 was the only thing blocking these |

So Layer 3 was uniquely protecting **6 of 234 loci**, and for the other 8 it was
redundant with a structural field. Recorded as what is being given up, not as a
gate.

---

## 3. Why removing polarity does not cost target accuracy

Two things must be true for an empty target site to be right:

1. **The middle really is an element**, not host sequence deleted in one genome.
   → the **M line** covers this. Host deletion is a single-copy event; an insert
   with ≥3 dispersed copies is a mobile element. Layer 3 says nothing about it.
2. **The short allele faithfully represents the sequence before arrival.**
   → the **structural fields** cover this: `target_bases_lost`, signed overlap,
   placement unambiguity.

Layer 3 answers a third question that is neither of these.

And getting direction wrong does not automatically corrupt the target: for a
scarless element or a precise excision, the short allele is still the true target
sequence. Only *imprecise* excision corrupts it — and that is exactly what
`target_bases_lost > 0` catches, as the 8 rows above demonstrate.

---

## 4. Fields, never gates

| field | source | meaning |
|---|---|---|
| `insertion_point_offset` | #2 | where the insert sits inside the target |
| `overlap` (signed) | #2 | `> 0` a direct repeat of that length (a TSD, when it is one); `= 0` point insertion; `< 0` target bases destroyed |
| `target_bases_lost` | #2 | `max(0, -overlap)` |
| `decomposition_method` | #2 | `exact` or `tolerant_alignment` |
| `placement_margin` | #2 | runner-up margin from `place_locus()` |
| `insert_copies_in_genome` | M | Arm A copy count |
| `insert_distinct_loci` | M | dispersion |
| `insert_md5` | #2 | dedup key component |

**TSD is a field, not a stage.** It is `overlap > 0`, derived from the
decomposition, free the moment #2 finishes. Nice when present, absent for whole
families (IS110 excises scarlessly). No pre-registration, no null model, no power
floor, no length histogram, no family matching. The column is there if anyone
wants to plot it later.

> **A caveat that survives, because it scopes a field rather than a claim:**
> `shortest_allele_len` is the **200 bp pad**, not a biological target interval —
> it equals the seed interval width in **95.5%** of E. coli loci (median
> difference 0, range 300–500 bp for 126 of 133, **zero loci below 100 bp**).
> `target_retained_bp` is therefore the pad too and measures nothing biological.
> `target_bases_lost` is unaffected: 95 non-zero values up to 5,517, computed by
> the tolerant alignment.

---

## 5. The only quality gate: S1

```
S1  =  target_bases_lost == 0
   AND overlap <= 15                  (above 15 is a repeat array or segmental
                                       duplication, not a junction feature)
   AND placement unambiguous
```

E. coli: **117 of 133** decomposable loci pass. With the M line, **77 usable**.

`overlap == 0` is **not** required. The old `junction_ambiguity == 0` was
backwards — it excluded every element that produces a direct repeat, which is
the majority (115 of 133 have `overlap > 0` after the cap fix).

---

## 6. Tests retained — Layer 1 and Layer 2 only

**Layer 1** (`30_armB_graph.py`), all verified:

| assertion | behaviour |
|---|---|
| every manifest entry exists | exit 2 |
| every entry above `--min-genome-bytes` (default 1 MB) | exit 2 |
| `graph.gfa` non-empty with ≥1 S line | exit 2 |
| bubble count 0 | exit 3 unless `--allow-zero-bubbles` |

**Layer 2** (`73_layer2_regression.py`), nine frozen gates: the eight frozen
2026-08-28 plus the **overlap-distribution gate** added 2026-09-02, because 99
loci changed `junction_ambiguity_bp` and all eight originals returned
bit-identical numbers — proof no gate covered the field.

**Plus, in any join:** assert a nonzero match rate, or print matched/total and
fail loudly below a stated floor. Four silent-negative failures this session —
the item-11 GCA/GCF key, `cat | awk 'NR>1'` counting 199 headers as data,
`split(".")[-1]` returning M1 = 0, and Arm B's empty input writing a header and
exiting 0. **All four were caught by human domain intuition; none by code.** At
census scale that intuition is not available.

That is the complete test surface. Nothing else.

---

## 7. Tests removed

| removed | was |
|---|---|
| TSD spike test, U-shaped null, critical values, per-species *n* floor, O1–O5, canonical-length list | a validation programme built on a field that is now just a column |
| CT core / target-motif test | an external precision claim; precision is out of scope |
| per-species mobility false-positive measurement **as a requirement** | run it if curious; not a gate |
| clade-rule invariance gate, stratified concordance | Layer 3 machinery |
| verification-overlap panel construction (18 → 60–70) | Layer 3 machinery |
| the read gold on `igi.biotite` | the only route to per-locus precision, which is out of scope. **Closed, not pending.** |

**What this costs, stated once and not revisited:** accuracy claims stop at
*structurally consistent and mobility-supported*. The catalogue is not per-locus
verified, precision is not quantified, and E. coli's 1.2% mobility false-positive
rate does not transfer to other species. Every field needed to calibrate later is
written to the record, so this is recoverable if the goal changes. **It is a
scope decision, not an oversight.**

---

## 8. Panel construction collapses

The largest simplification, and easy to miss.

Every panel requirement existed to serve Layer 3: the focal clade needed local
monophyly plus a separation gap, so panels had to be phylogenetically shaped
objects. Without Layer 3 a panel needs one property:

> the genomes are close enough for anchors to place — above the ~95% ANI cliff.

- §9b gap/floor calibration — **not needed**. Only an ANI threshold above the
  placement cliff, from each species' own distribution.
- §4.3 clique growth, the 0.104 usable fraction, degree-does-not-predict — **no
  longer constraints**.
- §6 focal-clade/panel decoupling — **irrelevant, permanently**.
- The 3 K. pneumoniae panels that missed k=24 — **not a problem**.

Panels stop being designed and become "a set of sufficiently close genomes".
82.1% of E. coli genomes have >50 neighbours at ANI ≥99, so selection is not a
bottleneck once the phylogenetic requirement is gone.

---

## 9. Per-species stages

```
 1  inventory        assembly_summary_genbank (no accession dedup for an
                     NCBI-sourced pool -- GCA only)
 2  all-vs-all ANI   sharded `skani dist`, never `skani triangle`
 3  threshold        quantile of this species' own distribution, above the
                     ~95% placement cliff. Record the value used.
 4  panels           `91_seed_panel_sampler.py`, mutual comparability
 5  Arm A            `-f 0` (default masking deletes the highest-copy
                     elements) -> mobility families
 6  Arm B            minigraph bubbles -> candidates; >5 kb parked
 7  Layer 2          -> structural catalogue                        LOCKED
 8  M join           -> mobility fields
 9  dedup            (insert_md5, target_context_md5)               TO BUILD
10  catalogue        S1 flag + all fields
11  Layer 5          annotation, optional, strictly last
```

**Two parameters are per-species and must not be inherited:** the ANI threshold
(99.0 sits on K. pneumoniae's mode; median neighbours 2,383 against tens in
E. coli) and the 5 kb scope cap (if a species' dominant elements are larger this
silently truncates the main signal; the `armB_large_events.tsv` parking mechanism
stays).

---

## 10. Where the numbers stand

| | E. coli (`10_`) | E. coli (`91_` sweep) | K. pneumoniae (`91_`) |
|---|---:|---:|---:|
| Arm B events <5 kb | — | — | **3,743** |
| parked >5 kb | — | — | 1,275 |
| Layer 2 resolved | 234 | 519 | **3,486** |
| decomposable | 133 (56.8%) | 373 (71.8%) | **2,940 (84.3%)** |
| per-panel median | — | 79 | 80.5 |
| S1 clean | 117 | — | pending |
| mobility positive | 83 | — | pending |
| usable | 77 | — | pending |

The `10_` column is single-linkage and **not comparable across species**. Use the
`91_` columns, **per panel, never as totals** — panel counts differ.

The three-way decomposable comparison separates two effects: **panel
construction** accounts for 56.8 → 71.8, **species difference** for 71.8 → 84.3.

**Interim array snapshots are unusable for quantiles.** Large panels finish last,
so a partly-complete array is a biased sample. Measured: at 18/40 the median read
74; the completed 40 gave 84.5, and taking the first 18 *by directory name* gives
81 — the early finishers were the small panels. Report only total and n until an
array completes.

---

## 11. Immediate queue

1. ~~Read the 14 `deletion_candidate` rows' `target_bases_lost` and overlap~~ —
   **DONE 2026-09-02: S1 already stops 8, Layer 3 uniquely protects 6.**
2. Delete the Layer 3 chain and `92_`
3. Rewrite `86_` as the catalogue writer; apply S1 as a **flag, not a filter**
4. Run the M join on K. pneumoniae — the nonzero-match assertion's first real
   test, since no "IS1 is obviously multi-copy" intuition exists for this species
5. Build the event-level dedup key — does not exist; overlapping neighbourhoods
   will produce ~10⁶ redundant rows without it
6. Update `PIPELINE_RECTIFICATION.md` and `NEXT_STEPS.md` to point here

**Not on the queue, and not to be re-proposed:** k=48; Pool O (below the
placement cliff, 30.8%); more E. coli for density (94.1% already have ≥5
partners); loosening the 0.95 identity gate (the 93 length-polymorphism loci
align end-to-end at ~85%; admitting them manufactures insertions); loosening
`--max-tsd 15`; the shared-target-site hypothesis at within-cluster E. coli
divergence (1 clean locus of 758).

---

## 12. One kept record

**`cl0002.B000006`** — 1,321 bp, 12 copies at 12 distinct loci, flank homology
0.0000, overlap 4 bp, unmatched by ISEScan. rRNA operons and REP elements do not
look like this. It is the one thing a library-first method could not have
returned, and it is the reason the annotation-free ordering was worth paying for
even after everything else was cut.

---

## 13. PRE-REGISTERED — mobility false-positive measurement (written 2026-09-07, before the data)

`mobility_positive` currently means "≥3 dispersed copies, or the same insert at
another event". rRNA operons, prophage, REP elements and segmental duplications
also satisfy that. The only measurement so far is **1.2% (1 of 83 M1 loci)** on
E. coli's 15-genome set.

**Criterion, fixed before looking:**

| outcome | reading |
|---|---|
| full-scale E. coli FP rate **≈ 1.2%** | the criterion is stable at scale; the K. pneumoniae run is confirmation only, and 75.4% usable stands |
| **> 5%** | the 1.2% was a small-sample artefact — the same phenomenon as the decomposable gap (12.5 → 4.5 points) and the Layer-1 bubble FP rate (3.4% → 1.6%). **75.4% must be revised down by the new rate** |
| 1.2–5% | intermediate; report the rate, revise, do not claim stability |

**Prerequisite, per standing practice 6 — report matched/total BEFORE the rate.**
`novelty_per_insert.tsv` is ISEScan output for the original **20-genome** batch;
the full run has **3,413** genomes. If the join covers only a small fraction, the
measured rate belongs to that fraction, not to the full set — and that fraction
is the earliest, smallest genome set, so the selection bias runs in an unknown
direction. **Coverage is reported first, and if it is low the measurement is
deferred until ISEScan is run on a random sample of full-run carriers.**

---

## 14. PRE-REGISTERED — does adding E. coli genomes still pay? (2026-09-09, before the data)

Retrospective rarefaction on the 400 panels already on disk. No download, no
re-run. E. coli currently uses **3,413 of 8,961 usable** genomes (38%), so
~5,500 are untouched.

### Criterion, fixed before looking

Marginal yield = new distinct events per panel. Compare the **300→400** band
against the **100→200** band:

| ratio | conclusion |
|---|---|
| **> 70%** | still climbing — more genomes are worth it |
| **< 40%** | saturating — go to a new species instead |
| 40–70% | marginal; decide on what the goal is, do not re-read the curve |

### Three curves, not one — because the event key can fragment

`event = (insert_md5, target_context_md5)`. A new genome brings new flanking
context, which can split one event into two purely because the context differs.
Some of any rise is therefore **key fragmentation, not discovery**. So plot:

| curve | saturating means |
|---|---|
| distinct `insert_key` | the **element repertoire** is complete |
| distinct `context_key` | the **target-site** repertoire is complete |
| distinct `event` | what is currently reported |

**The informative outcome is the DECOMPOSITION, not the total.** If insert-level
saturates while context-level keeps climbing, the reading is: *the element
repertoire is found, but those same elements keep turning up at new target
sites.* That is direct material for target-site preference — the founding
question — and it is worth more than any headline count.

### Second test, nearly free: are the unused genomes new territory or clones?

For every unused genome, its maximum ANI to the used set:

| result | conclusion |
|---|---|
| mostly ≥99.9 | clonal near-relatives of what is already sampled; adding them raises redundancy, not events |
| a substantial share <99.5 | unsampled lineages exist; more genomes is genuinely new ground |

E. coli's median reported ANI is 97.35 (K. pneumoniae 99.02), so a diffuse pool
is the prior — but that is a guess and the measurement decides.

---

## 15. PRE-REGISTERED — the E. coli targeted top-up (2026-09-11, before the run)

The rarefaction ANI test split the 4,310 unused E. coli genomes:

| max ANI to the used set | n | % |
|---|---:|---:|
| >=99.9 (clonal near-relative) | 2,088 | 48.4% |
| 99.5–99.9 | 1,067 | 24.8% |
| **<99.5 (candidate new lineage)** | **1,155** | **26.8%** |

Adding only the 1,155 is a **data-selected subset**, not a bulk expansion: the
2,088 clonal genomes would push redundancy above 1.98x while adding little.
All 1,155 are already on disk, so this costs no download.

### The falsifiable expectation, fixed before the data

> **If ANI distance predicts yield, panels seeded on these 1,155 must produce
> MORE NEW EVENTS PER PANEL than 34.8** — the observed 300->400 marginal on the
> existing 3,413 genomes. "New" means events absent from the 21,051 already
> catalogued, keyed on (insert_md5, target_context_md5).

| outcome | reading |
|---|---|
| **> 34.8 new events/panel** | ANI distance predicts yield; targeted top-up is the right expansion strategy for every later species |
| **<= 34.8** | **ANI distance does NOT predict yield.** A useful negative about the method: stop pre-selecting by ANI, and either expand in bulk or not at all |

The second outcome is the more valuable one, because it retires a selection step
that would otherwise be repeated on every species. It is reachable: these
genomes are more distant from the sampled set, but distance from what has been
sampled is not the same as carrying unseen elements, and a low yield would show
that directly.

Expected scale if positive: ~1,155 genomes -> ~48 panels -> ~1,700 new events,
about 8% on top of 21,051.

---

## 16. COMPUTE POLICY — downloads never go through SLURM (2026-09-11)

### The rule

> **Downloads run on the transfer node, by hand. Never `sbatch`. Never a
> compute node.**
>
> ```
> ssh lrc-xfer.lbl.gov
> cd /global/home/users/kh36969/fna_based_mgefinder_project
> ./tools/download_on_transfer_node.sh <species>    # or no arg for all queued
> ```
>
> `slurm/download_species.sh` is **RETIRED**. It must not be resubmitted.

### Why — measured, not assumed

`cf1` is **`OverSubscribe=EXCLUSIVE`**: a job is charged the **whole 64-core
node** regardless of `--cpus-per-task`. `--cpus-per-task=8` bills 64.

Downloads are pure network I/O and use ~0 CPU, so every download job billed a
64-core node to move bytes:

| download | CPU-hours billed | CPU actually used |
|---|---:|---|
| K. pneumoniae (5,916 genomes) | **225** | ~0 |
| A. baumannii (1,956 genomes) | **40** | ~0 |

The five-species queue (9,863 genomes, 35 GB) ran on the transfer node instead:
**0 SU, 0 failures, all five COMPLETE.**

### The same policy fixes a much larger waste: PACK THE NODE

Since the node is exclusive, requesting *fewer* CPUs saves nothing. The only way
to stop wasting is to **use all 64 cores**. Measured on this project since
2026-09-01: **73,054 CPU-hours over 660 tasks**, roughly two-thirds idle cores.

| stage | threads used | charged | CPU-hours |
|---|---:|---:|---:|
| `full_armB` | 16 | 64 | **60,571** |
| `kp_ani` | 8 | 64 | 6,632 |
| `full_armA` | 8 | 64 | 1,963 |
| downloads | ~0 | 64 | 265 |

Replacements, both written and to be used from S. aureus onward:

| script | packing | array |
|---|---|---|
| `run_packed_armB.sh` | **4 panels/node** at `-t 16` | 200 tasks -> 10 |
| `run_packed_armA.sh` | **8 genomes/node** at `-t 8` | 200 tasks -> 10 |

Expected: Arm B **~60,000 -> ~15,000 CPU-hours per species**. Across the five
queued species, ~300,000 -> ~75,000.

Single-threaded analysis steps (`100_event_dedup`, `86_catalogue`) burn a
64-core node to use one core. Batch them into ONE job rather than chaining four.

### Checklist for every new species

1. Manifest -> `download_queue/<slug>.tsv` (login node, awk only, seconds)
2. **Download on the transfer node** — not SLURM
3. ANI, Arm A, Arm B, Layer 2 via the **packed** scripts
4. dedup + M join + catalogue batched into one analysis job
5. Login node is for `ls`, `wc`, `head`, `squeue` — nothing else

---

## 17. PRE-REGISTERED — does k=24 earn its cost? (2026-09-12, before the data)

**Why this outranks any data-source question.** A k=24 panel needs enough
mutually-comparable genomes to form a 24-clique. That requirement, not the total
genome count, is what sets the addressable species universe:

| usable genomes >= | species |
|---|---:|
| 200 (comfortable for k=24) | **49** |
| 50 (feasible only at small k) | **194** |

If small k yields comparably per genome, the universe is **194 species, not 49**
— a 4x expansion that no amount of extra sequence data can buy.

### Measured on data already on disk

The saturation sweep ran 5 seeds x k=6/12/24 with NESTED panels (k6 subset of
k12 subset of k24), so the comparison is within-seed and paired.

**The cost unit is the GENOME, not the panel.** Genomes are what must exist,
be sequenced and be downloaded; panels are cheap. So the statistic is
**distinct events per unique genome**, deduplicated on
(insert_md5, target_context_md5) exactly as the species catalogues are.

### Criterion, fixed before looking

events-per-genome at k=12 as a fraction of k=24:

| ratio | conclusion |
|---|---|
| **>= 80%** | k=12 is preferable — same yield per genome at a quarter the panel compute, and the species universe opens to 194 |
| **50-80%** | k=24 earns its cost on yield but k=12 remains the right choice for species that cannot reach 24 |
| **< 50%** | k=24 is genuinely necessary; the universe stays at ~49 species and small-genome species are out of scope |

Also report k=6 on the same axis: if k=6 is within 80% of k=24, the floor for a
usable species drops to roughly 50 genomes.

**Caveat to carry:** the sweep predates the `decompose()` ambiguity fix and the
500 bp floor, so its absolute counts are not comparable to the current
catalogues. The k comparison is internally consistent because all three k values
ran the same code, and only the RATIO is used.

### 17a. CORRECTION before the S. aureus run (2026-09-12)

**The endpoint in §17 is wrong. It must be USABLE events per genome, not raw
events per genome.**

S1 requires **>=2 carriers on each side** of a locus. At k=6 that means >=2 with
and >=2 without out of six; at k=24 it is >=2 of twenty-four. Small k is
therefore materially harder to satisfy, so **per-genome yield can rise while the
usable rate falls**, and the net effect is unknown. The table in §17
(1.55 / 1.19 / 0.90 events per genome) is RAW yield and has not passed S1 — the
172.5% advantage at k=6 could be partly or wholly eaten by it.

The M line is unaffected: M1 is within-genome and M2 is cross-event, neither
depends on k.

### 17b. The S. aureus design — both k values, non-overlapping seeds

k=12 needs half the genomes per panel, so running both configurations costs
about what one k=24 run was going to cost. The choice sets the default for every
later species, so it is worth buying outright.

- **200 panels at k=24** — strictly comparable to E. coli / K. pneumoniae /
  A. baumannii. This is the **primary test**: does the method cross into
  Firmicutes? It must not be confounded by changing k at the same time.
- **200 panels at k=12** — same species, same ANI threshold, post-fix code,
  **non-overlapping seeds**, so the panels are INDEPENDENT rather than nested.
  That removes the §17 caveat that k=6 events are a subset of k=24 events.

Report for both: events, usable events, genomes used, and **usable events per
genome**. If k=12 holds up, the species universe is 194 rather than 49 — an
expansion no dataset can supply but one parameter can.

Per-species recalibration still applies: ANI threshold from S. aureus's own
distribution (2.8 Mb genome, GC 33%, different population structure), 5 kb cap
retained. SCCmec-class elements will park in `armB_large_events.tsv`, and
Stage 32 recovers only ~0.2% of the parked pile, so those are effectively out of
scope rather than deferred.

---

## 18. SETTLED — the ANI threshold is a COVERAGE parameter, nothing else

Two species, same result, so this is no longer a coincidence:

| species | threshold dropped | panel min pw ANI moved | eligible seeds gained |
|---|---|---:|---:|
| A. baumannii | 99.72 -> 99.00 (0.72) | **0.01** (99.730 -> 99.720) | +258 |
| S. aureus | 99.78 -> 99.00 (0.78) | **0.02** (99.800 -> 99.780) | +233 |

**Panel tightness is set by the greedy mutual-comparability clique growth, not
by the edge threshold.** A 24-clique lands inside one tight lineage whatever
looser edges are on offer.

> **Therefore: choose the ANI threshold on COVERAGE alone. Quality does not
> enter.** §9b of NEXT_STEPS (gap/floor calibration) is **formally retired**,
> not merely de-emphasised.

### And quantile transfer is a starting value, never a calibration — 4 cases

| species | distribution shape | quantile of ANI 99.0 | floor at q=0.9266 |
|---|---|---:|---:|
| E. coli | diffuse, median 97.35 | 0.927 | 99.13 |
| K. pneumoniae | single mode AT 99.0 | 0.386 | 99.59 |
| A. baumannii | bimodal, step at p75->p90 | 0.815 | 99.72 |
| S. aureus | step at p90->p95 (CC structure) | 0.791 | 99.78 |

Four species, four shapes, the quantile landing somewhere different each time.
**The actual calibration is: build 40 probe panels and read panel min pairwise
ANI.** That procedure has now run twice and settled the threshold both times.

## 19. TWO TIGHTNESS CONFOUNDS, and the one analysis that resolves both

Panel tightness is now known to vary systematically, and it is not controlled in
either headline comparison.

**Confound A — the k experiment.** k=12 panels come out TIGHTER than k=24
(S. aureus probe: 99.845–99.865 vs 99.780–99.800), because a smaller clique is
easier to find among close relatives. So the k=12 arm changes two things at
once: smaller panels AND more homogeneous ones. Homogeneity means fewer
polymorphic sites per panel — the OPPOSITE direction to the raw per-genome
advantage. The §17 figure of 172.5% at k=6 could be partly cancelled by it, and
the final numbers alone cannot separate the two effects.

**Confound B — the cross-phylum test.** S. aureus has the tightest panels of any
species (99.78–99.87 vs E. coli 99.62–100.00). If its per-panel yield comes in
low, tightness and biology both explain it — and distinguishing them is exactly
what the phylum-crossing test is for.

### The resolution, on data already produced

**Regress per-panel yield on panel min pairwise ANI, across all four species.**
E. coli has the widest panel spread (99.62–100.00), so it alone carries enough
internal variation to estimate the slope. With the slope in hand:

- S. aureus's yield can be tightness-corrected before the phylum comparison
- the k=12 arm can be stratified by tightness: if tight k=12 panels yield less
  than loose ones, the tightness effect is visible and separable; if yield is
  flat in tightness, the k effect is clean

Both joins already exist: `panels.tsv` carries `min_pairwise_ANI` per panel and
the catalogue carries a `panel` column. **No new experiment, one analysis.**

---

## 20. M1's limits, measured — and why the gate must NOT be changed

Three independent tests, all pointing the same way.

**Finding 1 — M1 under-reports OUTSIDE A SPECIES' MODAL SIZE BAND, and the
profile is NON-MONOTONIC.** Aligned rate by insert length, four species:

| band | E. coli | K. pneu | A. baum | S. aureus |
|---|---:|---:|---:|---:|
| <800 | 80.6% | 66.4% | 55.0% | **24.9%** |
| 1200-1599 | 85.9% | 86.4% | 81.4% | 87.2% |
| 1600-2399 | **28.2%** | 51.4% | **20.3%** | 54.1% |
| >=3500 | **9.5%** | **8.6%** | 21.7% | **0.8%** |

A dip at 1600-2399 with partial recovery at 2400-3499 cannot come from a
monotone coverage threshold. Say "outside the modal size band", never "long
elements".

**Finding 2 — `overlap > 15` is a STRUCTURALLY UNCOVERABLE class, not a
negative.** S. aureus: 15.3% of decomposable loci (vs 1.7-5.5% elsewhere).
Stratifying short inserts by repeat context:

| <800 bp | n | M1 aligned |
|---|---:|---:|
| overlap <= 15 | 876 | 34.6% |
| **overlap > 15** | 359 | **1.4%** |

25x. An insert that IS a fragment of a repeat array has no discrete element for
Arm A to match, **by construction**. Emit it with its own flag: not usable, not
discarded, and not counted as a mobility negative.

**Finding 3 — the no-hit failures are TRUE NEGATIVES. Do not relax the gate.**
76-91% of M1 failures produce no Arm A hit of any kind. Counting copies directly
in the carrier genome, bypassing Arm A, on modal-band failures:

| direct copies | % |
|---|---:|
| **0-1** | **84.7%** |
| 2 | 7.0% |
| >=3 (a real Arm A gap) | 8.3% |

**84.7% are genuinely single-copy.** The hypothesis that drove this check — "the
modal band holds dominant elements, which must be high-copy" — was wrong: the
modal band is where EVENTS cluster, not where high-copy FAMILIES cluster. A
common insertion size does not make each instance one of many copies in its own
genome.

> **Relaxing `--min-cov` or switching to containment would convert true
> negatives into false positives across four catalogues. The gate stays.**

### How to state the usable rate

**The usable rate is a CONSERVATIVE estimate: the true usable fraction is not
below it.** M1 is strict outside the modal size band and in repeat context.

But Finding 3 bounds how strict: 84.7% of the exclusions are correct, so M1 is
not killing much wrongly. What is genuinely excluded is the `overlap > 15`
class, and those are **structurally uncoverable by the mobility line, not
misjudged**.

So the cross-species spread of **59-87% is mostly real biological difference**,
with one annotation: S. aureus's low end includes 15.3% structurally uncoverable
loci. Do NOT present the whole range as unreliable.

## 21. PRE-REGISTERED — M. tuberculosis failure modes (before the run)

MTB can fail at three places and they mean different things. Naming them first
so the outcome is not rationalised afterwards:

| where it fails | meaning |
|---|---|
| **cannot form k=12 cliques** | a PANEL-level limit; the 194-species universe must be discounted |
| **cliques form, Arm B returns ~0 bubbles** | a real SPECIES property — MTB is famously low-diversity with little HGT. NOT a method failure |
| **bubbles exist, Layer 2 cannot decompose** | a METHOD problem; investigate |

The second outcome would add a dimension the census cannot see: **enough genomes
does not mean enough insertion polymorphism.** A species can clear the 50-genome
bar and still have nothing to find. If MTB lands there, the 194 figure needs a
second filter that no genome count can supply.

### 18a. QUALIFIER — §18 holds in DENSE pools only (2026-09-14)

§18 said the ANI threshold is a coverage parameter and does not set panel
tightness. Four species now, and the effect scales with pool density:

| species | genomes | tightness moved when the threshold dropped |
|---|---:|---:|
| S. aureus | 4,209 | 0.02 |
| A. baumannii | 1,956 | 0.01 |
| **E. faecium** | **1,038** | **0.125** |
| M. tuberculosis | 895 | 0.07 |

Mechanism: in a dense pool the greedy clique fills up near the threshold, so the
threshold never binds. In a sparse pool the clique must reach further down and
the threshold starts to bind.

> **Operational rule: below roughly 2,000 genomes, the ANI threshold is BOTH a
> coverage parameter AND a tightness parameter, and must be chosen on both.
> Above that, coverage alone.**

**This matters most exactly where it hurts.** The 194-species universe is by
definition the sparse end (50-200 genomes), so "choose on coverage alone" is
WRONG for the bulk of the species the k=12 result was meant to unlock.

MTB (895 genomes, 0.07) breaks strict monotonicity with E. faecium (1,038,
0.125), so density is not the only term: MTB's ANI 99.0 sits at quantile
**0.0000** -- nearly every pair is already above 99.0 -- so lowering the
threshold adds almost no edges. **Binding depends on distribution shape as well
as genome count.** Record the tightness cost in `<species>_params.tsv` either
way; it is needed later to explain yield differences.

### 21a. How to READ the MTB result when it lands

MTB has cleared failure mode 1 outright: **895/895 genomes can seed a k=12
panel, 0/40 panels short.** So the 194-species figure takes no discount from
panel formation -- clonality makes cliques EASIER.

Panel min pairwise ANI is **99.910**, far above any previous species (E. coli
99.62-100.00, S. aureus 99.78-99.87, A. baumannii 99.720), and median pairwise
ANI across the species is 99.920. Only mode 2 remains live.

> **Do not read Arm B's bubble count on its own. Read it against panel min
> pairwise ANI, on the §19 cross-species tightness/yield regression.**

- **MTB falls ON the curve** -> it is not a special species, just an extreme
  point on a relationship that holds everywhere. **That is the more useful
  outcome**, because it means a species' ANI distribution PREDICTS its yield --
  and that can be computed at census time, before any download.
- **MTB falls BELOW the curve** -> something specific to MTB beyond tightness.

The first outcome supplies the second filter the census cannot otherwise see:
enough genomes does not mean enough insertion polymorphism, and the ANI
distribution is what tells them apart in advance.

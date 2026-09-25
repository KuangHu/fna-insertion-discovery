# Pipeline Rectification — from gated pipeline to catalogue + verification

**Drafted 2026-09-02.** Supersedes the tier and gating design in `NEXT_STEPS.md`
§1–§3. Everything in §9 (census) is unchanged and remains the main line.
`PROJECT_SUMMARY.md` remains the record of how the pipeline got here.

> **Numbers were re-derived from the outputs on 2026-09-02.** A first
> verification pass proposed four corrections; **three of them were wrong and
> have been withdrawn** — the draft and `PROJECT_SUMMARY.md` were right. See §11,
> which records the failed audit rather than deleting it.

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
>
> **The same applies to K. pneumoniae, and more strongly.** It has ~37 usable
> panels against E. coli's 5. Comparing totals would compare 37 panels with 5.
> `k` is the same (24) but genome diversity per panel is not. **The comparable
> quantity is the per-panel median and its distribution, never the total.**


## 0. What forced the rectification

Four things changed during the 2026-08-31/09-02 review session, and together
they invalidate the gating design rather than any individual module:

1. **§6 retired.** Panel size is the wrong variable; the focal clade does not
   grow with the panel. The genealogy channel was reassigned from yield
   production to verification.
2. **§7 closed.** The shared-target hypothesis is not supported at
   within-cluster divergence in *E. coli*. 1 clean locus of 758.
3. **The TSD null is U-shaped.** 47% of chance calls land at 4 bp. IS3, IS5,
   IS21 and IS1380 are not detectable at any realistic *n*. Combined with the
   no-TSD families (IS110, IS200/IS605), the structure channel has **two
   systematic blind spots**, not one.
4. **Layer 3's P1 is 79.6% clade-constant.** For those loci the local genealogy
   contributed nothing; the answer came from outgroup state alone.

Neither polarity channel can serve as a mandatory gate without writing its own
family bias into the output. The old design made both mandatory in sequence,
which is why 234 Layer-2 loci became **20** downstream-ready targets. (The
Layer-4 catalogue holds 34 rows: **20 insertion targets + 14
`deletion_candidate`**. The 14 are Layer 3 calling the long allele ancestral and
Layer 4 correctly refusing to emit a target — design working, not data loss.)

---

## 1. The core reframing

> **Layer 2's output is the product. Everything downstream assigns confidence
> to it, and nothing downstream may suppress it.**

Layer 2 already emits what the project was founded to find: for each candidate
locus, the alternative alleles between homologous anchors, and the junction
offset inside the pre-event allele. `pre[:offset] + insert + pre[offset:]`
regenerates the derived allele at **median identity 1.0000, min 0.9915**, over
the 20 insertion targets. The metric is undefined — and blank in the file — for
the 14 `deletion_candidate` rows, which have no reconstructed insertion. That is
an insertion and its target site.

What Layer 3 and the TSD detector add is **direction** and **junction
verification**. Those are annotations on a record, not preconditions for the
record existing.

### Two catalogues, not one gated pipeline

| output | contents | gating |
|---|---|---|
| `structural_catalogue.tsv` | every Layer-2 resolved locus, with all verification fields populated or NA | none — emit always |
| `target_catalogue.tsv` | filtered view: polarity determined **and** long-allele-derived | filtered, not gated |

The refuse-rather-than-repair logic in Layer 4 was correct and stays — if
polarity says the long allele is ancestral, no target is emitted. The error was
that "no polarity" also meant "no output." A locus with a solid structure and
undetermined direction is still a defined length polymorphism at a defined site
with a defined boundary, which is most of what the founding question asked for.

---

## 2. Two independent axes, never combined into one gate

Conflating structural confidence with polarity confidence is what produced the
91.5% loss. They are separate measurements with separate failure modes.

### Axis 1 — structural confidence

| code | condition |
|---|---|
| **S1** | decomposed, `target_bases_lost == 0`, `junction_overlap_bp <= 15`, placement not ambiguous, junction ≥50 kb from nearest contig end, ≥2 carriers on each side |
| **S2** | decomposed (exact or tolerant), one or more S1 conditions unmet — recorded which |
| **S3** | length polymorphism, not decomposable (the 93-locus class: aligns end-to-end at ~85% identity, so short-plus-insert is the wrong model) |

> **S1 no longer requires `junction_ambiguity == 0`, and the old condition was
> backwards.** The overlap field IS the TSD: `overlap > 0` means the junction
> slides by that many bp because the element duplicated its target site.
> Requiring 0 therefore excluded **every TSD-bearing element** — the class with
> the strongest mechanistic evidence. Overlap is a *measurement* read against
> the U-shaped null (§9c of NEXT_STEPS), not a defect and not a certificate:
> `> 15` means repeat array or segmental duplication (9 loci) and is excluded;
> `4–15` is the usable TSD band (85 loci); `1–3` is below the chance floor
> (21 loci); `0` is a clean point insertion (18 loci).
>
> The implementation already did this — the 117 "reconstruction clean" loci were
> selected on `lost == 0 AND overlap <= 15 AND placement not ambiguous`. Only
> the written definition was stale.

### Axis 2 — polarity confidence

| code | condition |
|---|---|
| **D1** | two independent lines agree |
| **D2** | one line |
| **D0** | undetermined |
| **DX** | lines conflict — retain, do not resolve |

Independent lines, in descending scalability:

- **M — mobility (Arm A).** The insert sequence occurs at ≥3 dispersed loci with
  high internal identity and divergent flanks. Family-agnostic, no panel, no
  TSD dependence, already built (**1,610 families / 9,847 copies** over 200
  genomes, `-f 0` run; pre-fix 1,648 / 8,563). This establishes that the long allele contains a mobile
  element rather than an ordinary indel. **Never previously used as polarity
  evidence — the largest unexploited signal in the pipeline.**
- **T — TSD.** One copy becoming two is irreversible. Strong per locus,
  coverage **32.5% (76/234)**, with the two blind spots above.
- **X — target motif (external).** Family-specific published target
  specificity, e.g. the IS110/IS621 CT core. Validates the junction at
  nucleotide resolution. See §4.
- **G — genealogy (Layer 3).** Strongest per locus, not scalable. Reserved for
  targeted verification, see §5.

**DX is a retention class, not a failure class.** The two known DX cases had
Layer 3 calling the long allele ancestral against TSD calling it derived, and
both are legitimately ancient-insertion-then-loss. Layer 4 already emitted them
as `deletion_candidate`.

### The honest consequence, to be carried in the document

> The bulk of output will be S1/S2 + D2, and D2 by mobility evidence alone is
> per-locus unverified for direction. Accuracy rests on distribution-level
> evidence plus a small D1 anchor. **D1 cannot be produced at volume**, because
> two of the four lines (G, X) do not scale.

---

## 3. TSD, repositioned

TSD was doing three jobs. Only two of them are jobs it can do.

| job | verdict |
|---|---|
| boundary confirmation | **keep** — this is what it is good at |
| de novo validation (spike test) | **keep** — distributional, and it survived a proper null |
| polarity engine for the whole structure channel | **remove** — 32.5% coverage with systematic family bias |

Two corrections that follow, both documentation-only and **both already applied
on 2026-09-02**:

- **The *E. coli* 4–5 bp shoulder is downgraded.** It was read as IS3 alongside
  the IS1 peak at 9 bp. Under the U-shaped null, 5 bp is indistinguishable from
  chance (z = +0.35) and 4 bp is significantly **depleted** (z = −5.94) — the
  detector calls *fewer* 4 bp repeats than random placement produces. The peak
  is IS1; the shoulder cannot be claimed as IS3. 8 bp (z = +3.71) is suggestive
  and is not claimed either.
- **The 9 bp peak gets *stronger*.** Null 4.2%, not the 9.1% a uniform null
  implied; z = +10.18, the only bin clearing the ~3.9 critical value. 27.6% in
  the modal bin against a measured null is the strongest validation result the
  project has.

---

## 4. New: external target-motif validation (the CT core test)

The one route to a nucleotide-resolution junction claim that does not require
the read gold.

**Background.** Durrant, Perry et al. (Nature, 2024) showed IS110 elements
express a structured bridge RNA whose two loops base-pair with sequences
adjacent to a **central CT dinucleotide core** in both target and donor DNA.
IS621 is the characterised member and is native to some *E. coli* strains. IS110
excises scarlessly and produces no TSD — it is exactly the family the structure
channel is blind to.

**The test inverts the paper's ordering rather than borrowing its ordering.**
Never search for the CT core to find elements. Instead: take junctions this
pipeline placed structurally, with no motif knowledge, and ask whether they land
on CT.

**Procedure.**
1. Select the 1279 bp insert class from the existing *E. coli* Layer 4 output —
   the IS110 size, measured independently by Arm B and by Layer 2 decomposition.
2. Read the two bases at `offset` in the pre-event allele.
3. **Null, written before running:** CT frequency at randomly drawn offsets
   within the same pre-event allele, ≥1,000 draws per locus, true junction
   excluded ±5 bp. Chance CT is ~1/16, so without this the test succeeds by
   construction.
4. Statistic: excess of CT at the true junction over the per-locus null.

**Pre-registered outcomes.**

| outcome | reading |
|---|---|
| significant CT enrichment at `offset` | junction is nucleotide-accurate for this class; `NOT CLAIMED` in §5 becomes `CLAIMED for IS110-class` |
| no enrichment, but enrichment at `offset ± 1–2` | systematic junction offset — a real and fixable finding |
| no enrichment anywhere | either the junction is not nucleotide-accurate, or CT is IS621-specific rather than IS110-general |

**Confound to resolve first, from the paper's supplement:** whether the CT core
is general across the IS110 family or documented only for IS621. If only IS621,
the eligible set shrinks and the power calculation changes.

**Power caveat, to compute before running:** the eligible set is the 1279 bp
class within 34 catalogued events, so *n* is small. Compute the detectable
effect size before spending the test, per §10 rule 1.

**Why this matters for resource allocation.** The read gold on `igi.biotite` is
unreachable and has been the sole route to a precision claim. This is a positive
control from published biology, costs nothing, and validates precisely where TSD
is blind. It raises the priority of the CT test and lowers the urgency of the
read gold by one further notch.

**Distinct from the methodological comparison.** Durrant's MGEfinder is
read-based. Comparing this pipeline against it is a separate exercise requiring
the same samples through both, and is not scheduled here.

---

## 4b. LINE M — RUN 2026-09-02, and it is the largest single gain

`scripts/99_mobility_join.py`. Pure join on data already on disk, no new compute.

```
234  Layer-2 loci
133  decomposable (an insert exists to test)
117  reconstruction clean      88.0%   lost==0, overlap<=15, placement not ambiguous
 83  mobility positive         62.4%
 77  USABLE = clean AND mobile 57.9%   <- vs 20/234 under the old gated design
```

### M2 collapsed, and the lesson generalises

M2 was first written as "a byte-identical insert in **another genome**". At
>=99% ANI two genomes share almost everything byte-identically, so that tests
almost nothing:

| | naive (another genome) | corrected (another LOCUS) |
|---|---:|---:|
| M2 positive | 28 | **2** |
| M2 unique rescues | 19 | **0** |
| USABLE | 93 | **77** |

**26 of the 28 had the identical insert only at their own locus** — vertical
inheritance, not mobility — and **all 19 of M2's unique rescues were same-locus**.
Cross-genome identity carries essentially no information between close
relatives; only the *locus* constraint makes it evidence. **M1 carries the
entire line.**

### The false-positive test, which the module demanded and which now has an answer

`mobility_positive` means "≥3 dispersed copies", which rRNA operons, prophage,
REP elements and segmental duplications also satisfy. Measured against ISEScan
on the 83 M1 loci:

| | n |
|---|---:|
| **KNOWN IS** | **82** |
| UNMATCHED | 1 |

**False-positive rate 1/83 = 1.2%.** Families: IS1 23, IS110 19, IS4 18, IS3 12,
IS5 6, IS30 2, IS21 1, IS66 1.

**The control group is what makes this mean something.** The 50 decomposable
loci that are NOT M1 are 20/50 KNOWN. If both groups were near 100% the result
would only say "ISEScan recognises most things in E. coli". M1 is selecting.

**The three feared classes are excluded by mechanism, not by luck.** M1 insert
length is median 1279, range **776–2450**: zero above 4 kb (rRNA operons are
~5 kb), zero below 500 bp (REP/BIME), and prophage at tens of kb is held out by
the 5 kb cap. The length distribution does the work.

*Caveat:* ISEScan is an IS library, so this is a post-hoc false-positive
measurement against a known reference — exactly what Layer 5 is for. It is not
circular: Arm A never saw an IS library.

### `cl0002.B000006` — the one UNMATCHED locus, and it is a candidate not an error

```
len 1321 bp | 12 copies at 12 DISTINCT loci | flank_homology 0.0000 | overlap 4 bp
ISEScan: UNMATCHED
```

Twelve dispersed copies with completely divergent flanks is not what an rRNA
operon or a REP array looks like. **This is the only thing in the project so far
that the founding constraint — never start from a transposase — could have
produced and a library-first method could not.** It is recorded here so it does
not stay in a conversation log.

---

## 4c. Stage 31 is NOT on this path — and what that means for precision

**`31_armB_refine.py` cannot run here: `svim-asm` and `nucdiff` are both absent
from the environment.** This is an environment fact, not a defect, but it must
be recorded, because the script's existence implies it is on the chain and it is
not. Anyone reading the stage list will otherwise assume nucleotide-level
refinement happened.

**It is not needed, and the reason is a deliberate design choice.** The Layer-2
seed loci are derived directly from `armB_events.tsv`:

```
seed_genome = empty_ref_genome
contig      = empty_ref_contig
start       = junction_L - 200
end         = junction_R + 200
```

Verified against the E. coli `loci_armB_cl0000.tsv`: **100/100 loci reproduce
exactly.** The 200 bp pad is the point — §3.4 seeded loci with a pad
specifically so **the junction had to be relocated from sequence by Layer 2
rather than handed over**. Graph bubble coordinates locate the locus; they are
not a biochemical junction. So stage 31's nucleotide-level refinement was never
on the critical path: Layer 2 re-finds the junction itself.

> ### THE PRECISION CLAIM DOES NOT TRANSFER TO K. PNEUMONIAE
>
> K. pneumoniae junctions are produced by **exactly the same machinery** as
> E. coli's, with **no additional nucleotide-level evidence** of any kind. The
> §5 boundary — *nucleotide-resolution junction: NOT TESTED, NOT CLAIMED* —
> applies to K. pneumoniae verbatim, and more so: **K. pneumoniae has not even
> had the post-hoc ISEScan comparison** that gave E. coli its 1.2%
> false-positive figure.
>
> "40/40 panels, zero failures, comparable magnitude" is a statement about
> **enumeration**, not about accuracy. Nothing in the K. pneumoniae run so far
> is evidence that its junctions are placed correctly.

---

## 5. The genealogy channel's new job

Reassigned in the previous session and unchanged here: **produce verification
overlap, not yield.**

Current overlap between TSD calls and Layer 3 P1 is **18 loci**, of which only 5
are clade-polymorphic. That is the entire independent-reference base for 76 TSD
calls. The cheapest way to strengthen any accuracy claim is not more k — it is
**panels constructed specifically around loci that already have a structural
call**, pushing the overlap toward 60–70.

Requirement carried forward: overlap concordance must be **stratified by
clade-polymorphic vs clade-constant**. An invariance or concordance test run on
a mostly clade-constant set passes almost automatically and certifies nothing.

Consequence: §6's open question — should the focal clade grow with the panel —
affects only how fast the verification set grows. It blocks nothing and is not
scheduled.

---

## 6. Stage order, per species

```
 1  inventory            assembly_summary_genbank; NO accession dedup needed
                         for an NCBI-sourced pool (see §9a of NEXT_STEPS)
 2  all-vs-all ANI       sharded `skani dist`, never `skani triangle`
 3  threshold calibration  gap/floor test, not quantile match  (§9b)
 4  panel construction   clique at calibrated floor — TIME ONE PANEL FIRST
 5  Arm A                self-align with `-f 0` → mobility families  [line M]
 6  Arm B                minigraph bubbles → candidate loci
 7  Layer 2              anchors, tolerant decompose, place_locus
                         → structural_catalogue.tsv                 LOCKED
 8  verification         TSD [T] · mobility join [M] · target motif [X]
                         → fields on the record, never filters
 9  dedup                event-level key (§7)
10  tier assignment      S-code × D-code
11  target_catalogue     filtered view of 8+10
12  genealogy channel    OPTIONAL, targeted subset, for verification overlap
13  Layer 5 annotation   strictly last
```

Layer 3 moves from step 7-and-mandatory to step 12-and-optional. That is the
whole of the scalability change.

**Per-species parameters that must be recalibrated, never inherited:**

| parameter | why it cannot be a constant |
|---|---|
| ANI threshold | ANI 99.0 sits at *K. pneumoniae*'s mode — quantile 0.3864 there vs 0.9266 in *E. coli*; median neighbours 2,383 vs tens |
| TSD null shape and *n* floor | the U-shaped null was computed from *E. coli* allele-length distributions; a different distribution moves the critical value and the power floor |
| 5 kb scope cap | if a species' dominant elements are larger, this silently truncates the main signal — the `armB_large_events.tsv` parking mechanism stays |

Record every chosen value with its achieved floor/gap in
`census/<species>_params.tsv`, or inter-species yield differences will be
uninterpretable.

---

## 7. Scale-out requirements that only appear at volume

1. **Event-level dedup.** The same element at the same site will be captured by
   thousands of overlapping neighbourhoods. Key: `(insert md5 or
   length+TSD, target-context md5)`. Without it the output is ~10⁶ redundant
   rows. Does not exist yet.
2. **Single-carrier is the majority case at scale.** All 20 current
   downstream-ready targets have exactly one derived carrier. The collapse
   hypothesis was retracted on mechanism (9 bp ≪ assembler *k*; long-read
   assemblies; junctions 45 kb–2.2 Mb from contig ends), so this is not a
   correctness alarm — but the ≥2-carriers-per-side condition belongs in S1
   rather than being left to downstream judgment.
3. **Key consistency.** Item 11's ANI lookup failed because dedup preferred GCF
   while carriers were named GCA. Manifest, panel builder and Layer 2 must
   derive from the same accession key.

---

## 8. Work queue

### Phase 0 — documentation only, no compute
- [ ] Rewrite tier definitions on the S×D two-axis model (§2)
- [x] Downgrade the 4–5 bp shoulder claim — **done 2026-09-02** in
      `NEXT_STEPS.md` §4 and §9c. `PROJECT_SUMMARY.md` was checked and does
      **not** carry the claim, so no edit was needed there.
- [x] Record the structure channel's **two** systematic blind spots alongside
      "the bulk of output is per-locus unverified" — **done 2026-09-02**,
      `NEXT_STEPS.md` §3.
- [ ] Restate the §5 capability boundary: Layer 2 emits structure; direction is
      a separate annotated field with its own confidence
- [ ] **Scripted re-export audit of every table number in
      `PROJECT_SUMMARY.md`** — as a script with a reproduction check, never ad
      hoc shell (§11 is what ad hoc shell produced). Change numbers only, never
      conclusions, and record the original value at each edit.

### Phase 1 — cheap tests on existing *E. coli* output
- [ ] **CT core test** with its null, power calculation and pre-registered
      outcomes (§4) — highest value, and the only unblocked route to a
      precision claim
- [ ] **Arm A mobility join** across all 234 loci → populate line M. This is
      the field that converts D0 loci into D2 and is the single largest
      expected gain in usable output
- [ ] Queue item 12: TSD re-search in a ±20 bp window around middle boundaries
      (the earlier 0/8 negative probably searched the wrong window)
- [ ] Queue item 13: coverage on the max-0.683 revcomp pair

### Phase 2 — census (running; unchanged)
- [x] Last ANI shard → calibration on the complete set — **done 2026-09-02**:
      34,535,778 pairs, chosen floor **ANI 99.570** (q = 0.9266, retains 7.34%)
- [ ] Time one panel before launching the full set
- [ ] Panels → Arm B → Layer 2
- [ ] **Recompute the TSD null for *K. pneumoniae*** — do not inherit *n* ≥ 30
- [ ] Spike test with pre-registered O1–O5 and the canonical-length list

### Phase 3 — scale-out
- [ ] Event-level dedup key
- [ ] `structural_catalogue.tsv` / `target_catalogue.tsv` split
- [ ] Tier assignment as code, not judgment

### Phase 4 — unblocked but unscheduled
- [ ] **`cl0000` returned 0 bubbles from 9 genomes, 248 from 3.** Undiagnosed,
      upstream of every locked module. First look: rGFA `SN`/`SR` tags, to see
      whether all 9 genomes were integrated or minigraph dropped divergent
      members during construction.
- [ ] Stage 32 — nested sub-5 kb search inside the parked oversized bubbles.
      Designed, not written. Also the only place §7's closed hypothesis could
      still be visible.

### Branch — read gold on `igi.biotite`
- **Available** → tier calibration; "tier calibratable" becomes "tier
  calibrated"; per-locus precision claims permitted.
- **Unavailable** → all precision wording stays at distribution level
  permanently. Recording this as the answer is a legitimate close, not an open
  item. The CT test now covers part of what the read gold was for.

---

## 9. Not on the queue — rejected on evidence

Every item here was ruled out by a test that could have gone the other way.
Nothing underpowered or confirm-only belongs in this list.

| item | evidence |
|---|---|
| k = 48 | k is the wrong variable; nonmonophyly 17 → 30 → 162 on the collapsing seed |
| Pool O (sister species, 90–93% ANI) | below the ~95% anchor-placement cliff — 30.8% anchor placement at <95% ANI |
| more *E. coli* for density | 94.1% already have ≥5 partners |
| loosening the 0.95 identity gate | the 93 loci align end-to-end at ~85%; admitting them manufactures insertions |
| loosening `--max-tsd 15` | 67% of random hits pile up at the ceiling before censoring |
| shared target site, within-cluster *E. coli* | 1 clean locus of 758 across 7 backgrounds (§7, CLOSED and scoped) |

**Scoped, not universal:** the shared-target negative applies to within-cluster
divergence in *E. coli* at ≥99% ANI, where little time has elapsed for a second
element to reach an occupied site. It is not a claim that the phenomenon does
not occur.

---

## 10. Standing practice

The durable output of the review session, and the reason the rectification is
trustworthy. Four rules, each earned by a retraction:

1. **Before spending a criterion, name the result that would falsify it and
   confirm that result is reachable.** Five failures shared one shape: a design
   that could only move one way — confirm-only, an uncalibrated cutoff, a
   control defined to pass, an inclusion condition that guaranteed its own
   conclusion.
2. **Score relatedness against a measured, composition-matched null**, never a
   hand-picked cutoff. The revcomp floor (median 0.488) and the U-shaped TSD
   null both overturned conclusions drawn against intuition.
3. **No aligned columns is not 0% identity.** A hard zero from a windowed
   comparison means no alignment, which is equally consistent with "unrelated"
   and "reverse-complemented."
4. **Codify, reproduce, then extend.** No claim from ad hoc bash. A script, a
   check that it reproduces known numbers, and only then new data. The
   reproduction failures are what found the entry-condition bug that survived
   five rounds of review.

A fifth, from this session's largest error: **a quantity that does not move when
you change *n*, tighten coverage, or swap the control may not be sensitive to
what you are measuring.** Invariance was read as strength; it was a warning.

6. **A silent negative result is a forbidden failure mode. Every join and every
   input must assert non-zero, or print matched/total and fail below a stated
   floor.** Four instances this session, all the same shape — a key or a file
   that fails to a *negative answer* rather than an error:

   | # | what | how it failed | caught by |
   |---|---|---|---|
   | 1 | item-11 ANI lookup | dedup preferred GCF, carriers were GCA | a lookup returning nothing |
   | 2 | `cat *.tsv \| awk 'NR>1'` | skips only the FIRST file's header; 199 headers counted as data | arithmetic not matching the summary |
   | 3 | `split(".")[-1]` on `GCA_x.1.A0000` | family key never matched → **M1 = 0** | "IS1 is obviously multi-copy" |
   | 4 | Arm B on 8 missing + 1 fake input | no graph, 0-byte BEDs, header-only TSVs, **exit 0** | a human reading a 9-vs-3 comparison |

   Every one was caught by human intuition about the domain, not by the code.
   At census scale that intuition is unavailable: a silently failed panel and a
   panel with genuinely no insertions are **identical in the output**, and
   nobody inspects hundreds of them one by one. `30_armB_graph.py` now enforces
   this with three assertions (manifest entries exist and exceed
   `--min-genome-bytes`; the GFA is non-empty with ≥1 S line; zero bubbles is
   FATAL unless `--allow-zero-bubbles` is passed).

   > **The rule paid for itself the day it was added, on a fault introduced by
   > the person adding it.** The first K. pneumoniae Arm B submission put all 40
   > panels through a sbatch whose `PATH` carried the project env (minigraph,
   > gfatools) but not the conda env (samtools, python3). Every panel died with
   > `FileNotFoundError: 'samtools'` and a full traceback, within minutes.
   > **Under the previous behaviour those 40 panels would have written
   > header-only TSVs and exited 0**, and a samtools-less panel would have been
   > indistinguishable from a panel with no insertions. The fix was one PATH
   > line plus a `command -v` loop over the four required tools.
   >
   > **All four historical instances were caught by human domain intuition, none
   > by code.** That intuition does not transfer: nobody on this project can
   > look at a K. pneumoniae element spectrum and think "IS1 is obviously
   > multi-copy". `99_mobility_join.py`'s non-zero-match assertion was written
   > after the fact on E. coli, where the answer was already known —
   > **K. pneumoniae is its first real test.**

---

## 11. The verification pass that failed — recorded, not deleted

A first pass proposed four corrections to the draft. **Three were wrong.** The
draft and `PROJECT_SUMMARY.md` were correct; the audit introduced the errors.

| # | audit claimed | truth | what went wrong |
|---|---|---|---|
| 1 | Arm A is 1,809 / 10,046, not 1,610 / 9,847 | **1,610 / 9,847 is correct** | `cat *_families.tsv \| awk 'NR>1'` skips only the FIRST file's header. 199 of the 200 headers were counted as data: 1,610 + 199 = 1,809; 9,847 + 199 = 10,046 |
| 2 | regeneration median is 0.9941 over 34, not 1.0000 | **1.0000, min 0.9915, over 20 is correct** | the 14 `deletion_candidate` rows have a BLANK identity field; sorting numerically treated blank as 0 and moved the median. Blank is not zero |
| 3 | 34 = 20 non-truncated + 14 context-truncated | **34 = 20 insertion + 14 `deletion_candidate`** | `context_truncated` is blank for deletion candidates because it does not apply; blank was misread as a distinct truncation state |
| 4 | "§7a of NEXT_STEPS" should be §9a | **correct** | genuine cross-reference error |

**All three failures are the same failure**: a field that is *absent* was read as
a field that is *present with a low value*. That is standing-practice rule 3
(§10) — "no aligned columns is not 0% identity" — in two new costumes, plus an
awk stream-vs-file error. Rule 3 is hereby generalised:

> **3. An absent value is not a low value.** Blank, NA, "no alignment" and "not
> applicable" must never enter a numeric comparison as 0. Check for the empty
> case explicitly before sorting, averaging or thresholding.

**The premise that motivated the audit is also withdrawn.** The claim that "3 of
4 numbers sampled from `PROJECT_SUMMARY.md` are wrong" was an artifact of the
audit's own bugs: 0 of 4 were wrong. `PROJECT_SUMMARY.md` §3.1 reproduces
exactly from `armA_pre_f0` (1,648 / 8,563 / ≥20:18 / ≥50:0 / max 37) and
`armA_f0` (1,610 / 9,847 / ≥20:57 / ≥50:8 / max 109) — all ten values match.

**A full re-export audit of `PROJECT_SUMMARY.md` is still worth doing** — the
baseline for every later species is quoted from it — but it must be run as a
script with a reproduction check, not as ad hoc shell. Ad hoc shell is what
produced this section. Added to Phase 0.

---

## 12. Verified numbers, re-derived from disk 2026-09-02

| claim | source | status |
|---|---|---|
| Arm A pre-fix 1,648 / 8,563 / ≥20:18 / ≥50:0 / max 37 | `armA_pre_f0/*_families.tsv` | ✔ exact |
| Arm A post-fix 1,610 / 9,847 / ≥20:57 / ≥50:8 / max 109 | `armA_f0/*_families.tsv` | ✔ exact |
| 234 Layer-2 loci | 94 + 128 + 12 | ✔ |
| TSD coverage 32.5% | 76 / 234 | ✔ |
| TSD × Layer-3 P1 overlap = 18 | `structpol_calls.tsv` | ✔ |
| Layer-4: 20 insertion + 14 deletion_candidate, all P1 | `pre_insertion_targets.tsv` | ✔ |
| regeneration median 1.0000, min 0.9915 (n=20) | field 28 | ✔ |
| nonmonophyly 17 → 30 → 162 | saturation sweep | ✔ |
| 30.8% anchor placement below 95% ANI | outgroup distance test | ✔ |
| 94.1% with ≥5 partners | density report | ✔ |
| K. pneumoniae ANI floor 99.570 at q = 0.9266 | `census/kpneumoniae_params.tsv` | ✔ |

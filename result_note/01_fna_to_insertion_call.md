# Workflow 1 — from `.fna` to an insertion call

**Input:** assembled genome FASTA. Nothing else. No reads, no annotation, no IS
library at any stage.

**Output:** one row per distinct insertion event, carrying the empty target
sequence verbatim and the offset at which the element sits inside it.

---

## The whole path

```
 1  inventory        assembly_summary_genbank -> species census
 2  download         transfer node, by hand -- never SLURM
 3  all-vs-all ANI   sharded `skani dist`
 4  threshold        calibrate on this species' own distribution
 5  panels           mutually-comparable cliques, k = 12
 6  Arm A            self-alignment -> within-genome multi-copy families
 7  Arm B            minigraph bubbles -> candidate loci
 8  Layer 2          anchors + decomposition -> target, offset, insert
 9  dedup            collapse panel redundancy into distinct events
10  M join           mobility fields onto each event
11  catalogue        every event, with flags; nothing filtered
```

---

## 1. Inventory — `95_species_census.py`

Parses `assembly_summary_genbank_bacteria.txt`, keyed on `species_taxid` (never
`organism_name`, which carries strain suffixes and shatters a species).

**Rank by USABLE = Complete Genome + Chromosome, never by total.** Arm B compares
whole assemblies; a contig-level pool contributes fragmentation, not alleles.
Salmonella is the trap: largest species in GenBank at 626,212 assemblies, only
0.57% usable.

194 species have ≥50 usable genomes; 99,274 usable genomes exist in total.

---

## 2. Download — `tools/download_on_transfer_node.sh`

Run **by hand on the transfer node**, never through `sbatch`.

`cf1` is `OverSubscribe=EXCLUSIVE`, so a job is billed a whole 64-core node
regardless of `--cpus-per-task`. Downloads are pure network I/O and used ~0 CPU:
the K. pneumoniae download cost **225 CPU-hours**, A. baumannii **40**. The
five-species queue (9,863 genomes, 35 GB) ran on the transfer node at **0 SU,
0 failures**.

Concurrency 4 with exponential backoff. 128 parallel connections drew NCBI
rate-limiting and a 1.2% transient failure rate.

A partial download is invisible downstream — no stage checks pool completeness —
so the script reports `got N / want N` and says INCOMPLETE explicitly.

**Plasmid-only entries must be removed.** 65 GenBank records flagged "Complete
Genome" contain only a plasmid (55 in S. enterica, 7 in P. aeruginosa). Filter
at 50% of the species median, never an absolute byte floor: M. pneumoniae's
genome is 0.83 Mb and an absolute 1 MB floor rejects the entire species.

---

## 3. All-vs-all ANI

Sharded `skani dist`, **never `skani triangle`** — triangle buffers all output,
so a wall-clock kill loses 100% of the work (measured: 4,400 of 7,723 queries
after 10h47m, needing ~8h more).

---

## 4. Threshold — `98_species_ani_calibration.py`

**Quantile transfer is a starting value, never a calibration.** Four species,
four distribution shapes, the quantile landing somewhere different each time:

| species | shape | quantile of ANI 99.0 |
|---|---|---:|
| E. coli | diffuse, median 97.35 | 0.927 |
| K. pneumoniae | single mode AT 99.0 | 0.386 |
| A. baumannii | bimodal, step at p75→p90 | 0.815 |
| S. aureus | step at p90→p95 (CC structure) | 0.791 |

**The actual procedure: build 40 probe panels at two thresholds and read panel
min pairwise ANI.**

**Choose on coverage — with one qualifier.** Panel tightness is set by the
clique constraint, not the edge threshold, so the threshold buys coverage almost
for free *in a dense pool*. In a sparse pool it binds both:

| species | genomes | tightness moved when threshold dropped |
|---|---:|---:|
| S. aureus | 4,209 | 0.02 |
| A. baumannii | 1,956 | 0.01 |
| E. faecium | 1,038 | **0.125** |
| M. tuberculosis | 895 | 0.07 |

Below roughly 2,000 genomes, choose on both. This matters most where it hurts:
the 194-species universe *is* the sparse end.

---

## 5. Panels — `91_seed_panel_sampler.py`

A panel is a **clique**: a genome joins only if it clears the threshold against
every member already in, not merely against the seed. Single-linkage collapses
into one giant component because chaining is rampant — "degree 300" is not
"300 comparable genomes".

Among admissible candidates it maximises lineage diversity. Twenty resequenced
isolates of one outbreak supply one observation.

**k = 12, not 24.** Measured on S. aureus with disjoint seed pools, post-fix
code, 200 panels each: k=12 delivers **95.9% of k=24's usable events per genome
using 27% fewer genomes**. The endpoint must be *usable* events per genome — raw
yield reads 101.7% and overstates small k, because S1's ≥2-carriers-per-side
condition bites harder there.

This opens the addressable universe from ~49 species to ~194.

---

## 6. Arm A — `20_armA_multicopy.py`

Self-alignment per genome → families of ≥3 copies at distinct loci.
Family-agnostic; no library.

**`-f 0` is mandatory and not optional.** Default minimizer masking deletes the
highest-copy elements outright. Measured on one genome with 60 IS1 copies: the
default emitted 138 self-alignments and **zero** covering any IS1 locus; `-f 0`
emitted 3,442, with 57 covering one locus.

---

## 7. Arm B — `30_armB_graph.py`

`minigraph -cxggs` over a panel, then `gfatools bubble`. A bubble whose paths
differ by ≥ `--min-insert` is the empty/filled structure. Oversized events are
**parked** in `armB_large_events.tsv`, never discarded.

Four assertions, all of which have fired on real input:

| assertion | exit |
|---|---|
| every manifest entry exists | 2 |
| every entry above 50% of the panel median size | 2 |
| `graph.gfa` non-empty with ≥1 S line | 2 |
| bubble count 0 (unless `--allow-zero-bubbles`) | 3 |

These exist because a run with 8 missing inputs and one 2,015-byte file (an
insert sequence mis-saved as a genome) produced no graph, 0-byte call files,
header-only TSVs — **and exited 0**. At census scale a silently failed panel and
a genuinely invariant one are indistinguishable.

---

## 8. Layer 2 — `70_allele_reconstructor.py` — LOCKED

The product. Cut a 500 bp anchor each side, map both into every genome, extract
the interval **between** them verbatim, decompose the long allele against the
short one.

```
target_seq[:offset] + insert_seq + target_seq[offset:]  ==  observed derived allele
```

Loci are seeded with a **200 bp pad** deliberately, so the junction must be
relocated from sequence rather than handed over. Stage 31 (`svim-asm`) is
therefore not on the critical path — verified 100/100 against the E. coli seed
file.

**`--min-insert 500`, enforced on both decomposition paths.** Set from the
measured mobility cliff, not a guess: within-genome multi-copy rate is 1.3% at
400–599 bp and 71.8% at 600–799 bp. Cost, stated: MITEs (100–400 bp, genuinely
mobile) are out of scope — on structure alone they are not separable from REP
arrays and short indels.

### The signed overlap

`overlap = lcp + lcs − len(short)` is one signed number:

| | meaning |
|---|---|
| `> 0` | direct repeat of that length — a TSD when it is one |
| `= 0` | point insertion, target fully retained |
| `< 0` | target bases destroyed |

**TSD is a field, not a stage.** The overlap field reproduces the independent
k-sweep detector exactly — 76/76 calls, identical lengths — so `92_` was retired.

---

## 9. Dedup — `100_event_dedup.py` — mandatory before any count

Key: **(canonical insert, canonical ±100 bp target context)**, both
orientation-canonical via `min(seq, revcomp)`.

The context must be a **window centred on the offset**, not the whole short
allele: the short allele is essentially the 200 bp pad (it equals the seed
interval width in 95.5% of loci) and its boundaries shift between panels, so
hashing it under-merges.

Redundancy is real and species-dependent: E. coli 2.01×, K. pneumoniae 2.66×,
A. baumannii 3.40×, M. tuberculosis 4.55×, **M. pneumoniae 67.4×**; 2.54×
overall.

**The insert key must be junction-invariant.** A direct repeat at the junction
admits several equally valid boundaries, all reconstructing the identical
derived allele, and the alignment's choice is not symmetric under reverse
complementation — measured, re-deciding the junction in the revcomp frame keeps
the key for 79.7% of loci at overlap 0 and **0.0% at every overlap ≥ 1**. So
`lib_insert.canonical_insert_key()` slides the insert across the whole repeat
and takes the smallest canonical key over the equivalent placements. Before
this, 351 of 21,051 E. coli events (1.67%) were one event counted twice.

---

## 10. M join — `99_mobility_join.py` — requires `--events`

| line | question |
|---|---|
| **M1** | is the insert at ≥3 dispersed loci inside its own carrier genome? |
| **M2** | does a byte-identical insert occur at a **different event**? |

**M2 must key on event, not locus.** `locus_id` is panel-scoped, so the same
insertion seen by two panels reads as "the same insert elsewhere". Measured:
**85.7% by locus_id vs 35.2% by event_id** — a 50-point inflation. Dedup is
therefore a hard prerequisite.

M1's limits, measured and not to be patched:

- it under-reports **outside a species' modal size band**, and the profile is
  non-monotonic, so it is not a length effect
- of its no-hit failures, **84.7% are genuinely single-copy** — true negatives.
  Relaxing the gate would convert them into false positives
- inserts in repeat context (`overlap > 15`) are **1.4% M1 by construction** —
  a fragment of a repeat array has no discrete element to match

---

## 11. Catalogue — `86_catalogue.py`

**Nothing filters.** S1, mobility and length are columns.

```
S1 = target_bases_lost == 0
 AND overlap <= 15          (above 15 is a repeat array, not a junction feature)
 AND placement unambiguous
```

`overlap == 0` is **not** required. The old condition `ambiguity == 0` was
backwards — it excluded every element that produces a direct repeat, which is
most of them.

---

## What a row asserts, and what it does not

**Asserts:** at this anchored site these genomes differ by an inserted segment;
`target_seq` is the empty allele verbatim; the insert goes in at `offset`.

**Does not assert:** which allele is ancestral (no polarity in this pipeline);
what family the insert belongs to (no annotation); that the junction is correct
to the base (assembly tier, never read-verified); that the insert is a mobile
element.

---

## Result: 18 species, 8 phyla, 74,067 events

**Rebuilt 2026-09-24** after a code review found that the insert sequence was
being sliced out of the long allele with a short-allele offset. See
`docs/REVIEW_0262ce3.md`. Two things changed in this table: three more species
landed (P. aeruginosa, H. pylori, C. jejuni), and **every species lost events**
because the corrected key merges duplicates the old one split.

| species | phylum | events | was | usable |
|---|---|---:|---:|---:|
| E. coli | Gammaproteobacteria | 20,734 | 21,051 | 77.5% |
| K. pneumoniae | Gammaproteobacteria | 11,507 | 11,794 | 76.5% |
| E. faecium | Firmicutes | 9,115 | 9,281 | 82.4% |
| A. baumannii | Gammaproteobacteria | 7,248 | 7,459 | 88.1% |
| S. aureus | Firmicutes | 6,235 | 6,375 | 57.9% |
| P. aeruginosa | Gammaproteobacteria | 3,565 | — | 53.9% |
| M. tuberculosis | Actinobacteria | 3,278 | 3,399 | 41.8% |
| S. enterica | Gammaproteobacteria | 2,906 | 2,944 | 61.3% |
| E. faecalis | Firmicutes | 2,401 | 2,440 | 65.1% |
| S. pneumoniae | Firmicutes | 1,893 | 1,944 | 53.4% |
| H. pylori | Campylobacterota | 1,645 | — | 22.8% |
| B. subtilis | Firmicutes | 1,338 | 1,367 | 37.8% |
| L. monocytogenes | Firmicutes | 1,040 | 1,084 | 27.4% |
| N. gonorrhoeae | Betaproteobacteria | 716 | 755 | 9.9% |
| C. jejuni | Campylobacterota | 355 | — | 23.7% |
| T. pallidum | Spirochaetes | 54 | 54 | 1.9% |
| C. trachomatis | Chlamydiae | 26 | 28 | 42.3% |
| M. pneumoniae | Mollicutes | 11 | 11 | 18.2% |

`usable` = S1 **and** mobility-positive, per event. B. pertussis is still
running. `sa_k12` is the k=12 arm of the panel-size experiment and is excluded
here to avoid double-counting S. aureus; it holds 4,661 events (was 4,701).

Across the 18 species the correction removed **1,666 duplicate events (-2.2%)**
and no species gained one, which is the direction a de-duplication fix must
move. It also raised the mobility-positive rate from **71.5% to 78.9%**, because
M2 now recognises the same element on opposite strands.

**Every row records how its insert was obtained.** `insert_frame_status` is
`coord` when Layer 2 wrote the long-allele coordinate directly, `repaired` when
it was recovered by matching `inserted_md5`, and `unresolved` when it could not
be verified — in which case no sequence is exported at all. In this rebuild:
22,972 repaired, **0 unresolved**, and `md5(insert_seq) == insert_md5` for all
192,386 rows.

### Genome count does not predict insertion polymorphism

M. pneumoniae had **385 usable genomes, 200 panels, 0 panels short of k, 0
zero-bubble panels, Layer 2 decomposing at 72.9% — and 11 distinct events.**
Every stage succeeded. T. pallidum (54) and C. trachomatis (28) repeat it.

This is a filter the census cannot see: a species can clear any genome threshold
and still have almost nothing to find, and only running the pipeline reveals it.

### Layer 1 precision, self-contained

`monomorphic + substitution_only` — bubbles Arm B called where the two alleles
are equal-length or identical — run **1.2–1.6%** across species at full scale.
That is a reportable Arm B false-positive rate needing no external truth.

---

## Files

```
database/insertions.tsv            74,067 rows, 26 columns
database/insertions_targets.fna    empty target sites
database/insertions_inserts.fna    inserted sequences
catalogue/<species>_{loci,events}.tsv + FASTAs
full_<species>/                    panels, Arm A/B, Layer 2, dedup, mobility
```

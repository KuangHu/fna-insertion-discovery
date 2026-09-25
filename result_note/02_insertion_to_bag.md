# Workflow 2 — from an insertion call to a ready bag

**Input:** `database/insertions.tsv` — 69,986 (empty target site, insertion)
pairs from 15 species across 7 phyla.

**Output:** `bags/insertions_sites.jsonl` — 63,266 sites in 4,450 bags,
conforming to `CANONICAL_BAG_SPEC.md` of the DL project with departures stated.

---

## The path

```
1  merge      per-species catalogues -> one database, ids namespaced
2  CDS        prodigal -p meta, ab initio          -> ORFs
3  cluster    mmseqs on the DOMINANT protein       -> cds_cluster_id
4  bag        bag = CDS cluster
5  flank      120 bp of the EMPTY site, 60 | 60, junction at index 60
6  noncoding  insert minus its CDS spans
```

---

## 1. Merge — `110_build_database.py`

**`db_id = <species>.<event_id>`.** `event_id` is panel-scoped *within* a
species and collides across them. Namespacing is done once, here, at the point
where the collision would occur.

`sa_k12` is excluded so S. aureus is not double-counted; it remains available as
the separate k-experiment arm.

```
database/insertions.tsv         165 MB   69,986 rows, 25 columns
database/insertions_targets.fna  38 MB
database/insertions_inserts.fna 108 MB
```

---

## 2. CDS — prodigal, ab initio

```
120,661 ORFs over 69,986 inserts     1.72 per insert
 69,448 inserts carry >=1 ORF         99.2%
```

No HMM, no IS library. prodigal predicts ORFs from sequence alone.

---

## 3. Cluster — `111_cds_cluster.py`

`mmseqs easy-cluster` on protein, `--min-seq-id 0.50 -c 0.80`.

```
  6,095 clusters
        singletons 3,370   >=10 members 445   >=100 71   max 6,143
  dominant-ORF -> cluster join: 69,448/69,448 = 100.0%
```

**Protein, not nucleotide.** `106_element_families.py` already clusters the
inserts by DNA; that answers "is this the same element", not "is this the same
transposase family". Elements of one family diverge far more in DNA than in
protein. Both groupings are kept and are expected to disagree.

**Named `cds_cluster_id`, deliberately not `transposase_id`.** Nothing in this
step identifies a transposase. A cluster is a cluster.

**Multi-CDS rule:** when an insert carries several ORFs, the cluster of the
**longest ORF** defines its `cds_cluster_id` — where an element has a
transposase it is normally the longest ORF, and passenger genes (resistance,
toxin-antitoxin) are shorter. This is a heuristic, so `dominant_orf_frac`
reports how much of the insert the chosen ORF covers and `n_orfs` reports how
often the question arises.

**Cross-species clusters, no library consulted:** 281 in ≥2 species, 79 in ≥3,
3 in ≥5.

---

## 4. Bag = CDS cluster

### The alternative was measured before it was rejected

The spec's §1 requires every site in a bag to share **identical**
`noncoding_regions`. Keying on CDS cluster departs from that, so the departure
was going to be excused by "a consumer can re-split on `nc_sequence_hash`".
That escape was measured first:

| | CDS-cluster bag | strict nc bag |
|---|---:|---:|
| bags | 4,607 | **50,379** |
| bags with ≥2 sites | 2,115 (45.9%) | **5,239 (10.4%)** |
| sites in those bags | 62,742 (96.2%) | **20,094 (30.8%)** |

Median strict bag size **1**. The re-split is available and yields a **~90%
singleton corpus** — technically an escape, practically close to worthless.
`nc_sequence_hash` is written per site so a consumer can still take it, with
that number attached.

**Bags span species by design.** Key by `(corpus, bag_id)`.

---

## 5. Flank — 120 bp of the empty site

```
flank = target_seq[offset-60 : offset+60]
        [0:60] | [60:120]        junction between index 59 and 60
```

**It comes from the EMPTY allele, not the filled one.** `target_seq` is the site
as it was before the element arrived, reconstructed from genomes that lack the
insertion, so the window holds no element sequence at either end. A flank carved
from the filled allele would leak the element's own termini into the input.

**Exact centring is enforced, never padded.** A site is emitted only when
`offset >= 60` and `len(target) - offset >= 60`. Padding a short window would
move the junction off 60|60 and the model would learn the padding.

---

## 6. Non-coding regions

The insert minus its CDS spans, runs of ≥20 bp kept. This is the ncRNA context
and a DEPLOY field.

---

## Three fields that were specified, and how each is grounded

### `empty_site_source` — evidence, not assumption

`observed` only when `carriers_empty` names at least one genome that actually
carries the empty allele. Sites without that evidence are **dropped, not
relabelled**. Nothing in this corpus is `tsd_reconstructed`; the enum exists so
a future source that does reconstruct can say so and be filtered on.

### `orient` — canonical, not `unknown`

The flank is stored in whichever orientation is lexicographically smaller of
`(flank, revcomp(flank))`; `orient` records which was applied. Deterministic, so
two records of the same site agree.

Result: **fwd 31,583 / rc 31,683** — the near-exact 50/50 expected of a
lexicographic rule on unbiased sequence. A lopsided split would have meant the
rule was picking up something real.

**Cost, stated:** revcomp swaps the two halves, so upstream/downstream identity
is not fixed afterwards. The junction stays at index 60 because the window is
symmetric, and `orient` makes the original recoverable. **It is not a strand
label and must not be used as one.**

### Drops — by reason × species, not a total

| reason | total |
|---|---:|
| `no_noncoding_region` | 3,839 |
| `window_would_not_fit` | 2,071 |
| `no_cds_cluster` | 538 |
| `too_many_N` | 272 |
| **total** | **6,720** |

`no_noncoding_region` is the largest: inserts whose CDS covers essentially the
whole element. That is a real property of compact IS elements, not a failure —
but it means **the bag corpus is depleted in the most compact elements**, which
must be known before anyone reads an nc-length distribution.

Full breakdown in `bags/insertions_drops.tsv`.

---

## Result

```
bags/insertions_sites.jsonl    86 MB    63,266 sites
bags/insertions_bags.tsv      174 KB     4,450 bags
bags/insertions_drops.tsv     1.5 KB
```

| | |
|---|---:|
| sites | 63,266 of 69,986 |
| bags | 4,450 |
| sites per bag | median 1, max 6,015, **719 bags with ≥5 sites** |
| cross-species bags | 244 in ≥2 species, 72 in ≥3, **21 spanning ≥2 phyla** |

---

## §6 check 8 — scaffold sharing

```
distinct nc_sequence_hash    48,578
bags / unique nc hashes       0.092
```

The spec flags a **large** ratio as scaffold-sharing. This is the opposite
extreme: nc is almost entirely per-event, so there is no shared-scaffold
degeneracy of the kind behind Durrant's `guide_start_in_nc = 49`. Good for
localisation statistics.

It is also the re-split result seen from the other side: nc varies per site, so
a bag is held together by its CDS cluster and nothing else. **If the loader
assumes any nc coherence within a bag, this corpus will not satisfy it.**

---

## Conformance report (spec §7)

| | |
|---|---|
| **DEPARTURE — §1 bag-level nc invariant** | bags key on CDS cluster; sites within a bag may differ in `noncoding_regions`. Strict re-split measured: 50,379 bags, 10.4% with ≥2 sites, covering 30.8% of sites |
| bags span species | by design; key by `(corpus, bag_id)` |
| flank | `joined`, empty spacer, 60+60, junction at index 60, **no padding ever applied** |
| `orient` | canonical lexicographic, **not gold**; must not be read as a strand label |
| `labels.*` | **ABSENT** — unlabelled real observations, not generator positives. Must not be loaded as a labelled corpus |
| `empty_site_source` | `observed` for every emitted site, evidence-backed; sites without evidence dropped |

---

## Rebuild

Built from the 15 species complete at the time. P. aeruginosa, H. pylori,
B. pertussis and C. jejuni are still running; when they land, re-run
`110 → 111 → 112`. The scripts are idempotent.

# Workflow 3 — the rule a bag is built by

The parameters that decide what ends up in the same bag, stated once, with what
each was chosen for and what it costs.

---

## The clustering threshold

```
mmseqs easy-cluster  --min-seq-id 0.50  -c 0.80  --cov-mode 0
```

**50% amino-acid identity, 80% bidirectional coverage.**

| | |
|---|---|
| input | protein, not nucleotide |
| identity | ≥ 0.50 |
| coverage | ≥ 0.80, `--cov-mode 0` = over **both** query and target |
| clustering | `easy-cluster` (greedy set-cover), not connected components |

`--cov-mode 0` is the part that matters most and is easy to miss. Coverage is
required on both sequences, so a short protein matching a fragment of a long one
is **not** clustered with it. The alternative modes accept one-sided coverage and
would let a 120 aa fragment join the family of a 900 aa transposase.

Greedy set-cover rather than connected components: a chain of partial matches
A–B, B–C, C–D cannot drag A and D into one cluster when A and D share nothing.

**50% is this corpus's operating point, not a family boundary.** It is not
claimed to be the identity at which two transposases stop being the same family.
Nothing here measures that. It is the threshold the bags were built at, and any
statement about bag membership is conditional on it.

---

## Protein, not DNA — and both groupings are kept

`106_element_families.py` clusters the same inserts by **nucleotide**. That
answers "is this the same element". This step answers "is this the same protein
family". Elements of one family diverge far more in DNA than in protein, so the
two groupings are **expected to disagree**, and neither supersedes the other.

---

## What defines a bag's identity when an insert has several ORFs

**The longest ORF.** Its cluster becomes the insert's `cds_cluster_id`.

Rationale: where an element carries a transposase it is normally the longest
ORF, and passenger genes (resistance, toxin–antitoxin) are shorter.

**RATIFIED 2026-09-26 — do not change this.** The alternatives (cluster on all
ORFs, or on a concatenation, or split multi-ORF inserts across bags) were
raised together with the measurements below and the longest-ORF rule was kept
deliberately. Treat it as settled, not as an open heuristic awaiting a fix.

What it costs, stated so nobody has to rediscover it:

- for the 29.3% of inserts with two ORFs, the shorter ORF is predicted and
  recorded but **never clustered** — bag identity rests on one of the two genes
- the dominant ORF covers a median **85%** of a one-ORF insert but only **58%**
  of a two-ORF one, so for two-ORF inserts the bag describes a bit over half
  the element
- if the passenger gene is the longer one, the bag is defined by the passenger.
  Nothing in the pipeline can detect this, because no annotation is used
- the sort has no tie-break and Python's sort is stable, so on equal lengths
  the leftmost ORF wins

It remains instrumented rather than trusted:

- `dominant_orf_frac` — how much of the insert the chosen ORF covers, so a
  caller can see where the choice is weak
- `n_orfs` — how often the question arises at all

Measured on v3: 183,347 ORFs across 105,048 inserts, **1.75 per insert**, so the
question arises for most inserts and the field should be read.

---

## No annotation anywhere in this step

`prodigal -p meta` is ab initio. `mmseqs` is sequence-only. Nothing is matched
against a transposase HMM or an IS library.

The field is therefore named **`cds_cluster_id`**, deliberately not
`transposase_id`. Nothing in this step identifies a transposase. A cluster is a
cluster.

---

## The flank, which is a separate rule

```
flank = target_seq[offset-60 : offset+60]
        [0:60] | [60:120]        junction between index 59 and 60
```

**From the EMPTY allele, never the filled one.** `target_seq` is the site as it
was before the element arrived. A flank carved from the filled allele would leak
the element's own termini into the input.

**Exact centring, never padded.** A site is emitted only when `offset >= 60` and
`len(target) - offset >= 60`. Padding a short window would move the junction off
60|60 and the model would learn the padding.

### Multi-allele pairs are NOT all eligible for a flank

The multi-allele analysis produces two classes of pair and only one may supply a
bag flank:

| class | target | eligible |
|---|---|---|
| `empty_target` | the pair's shorter allele IS the locus's shortest | **yes** |
| `nested_target` | the shorter allele is an intermediate, so the target already carries an element | **no** |

Measured on v3: 101,508 novel pairs, of which 87,831 are `nested_target`.
Feeding those into a bag would put element sequence into a window that is
*defined* as pre-insertion. They are kept in a separate table.

---

## Orientation

The flank is stored in whichever orientation is lexicographically smaller of
`(flank, revcomp(flank))`; `orient` records which was applied. Deterministic, so
two records of the same site agree.

**`orient` is not a strand label and must not be used as one.** revcomp swaps
the two halves, so upstream/downstream identity is not preserved. The junction
stays at index 60 because the window is symmetric.

The near-exact 50/50 split is the check that the rule is doing nothing but
sorting: a lopsided split would mean it had picked up something real.

---

## v3 result at these settings

```
inserts                     105,048
ORFs                        183,347      1.75 per insert
inserts with >=1 ORF        104,159      99.2%
CDS clusters                  8,672
sites emitted                93,020
bags with >=5 sites           1,118
cross-species bags              308 in >=2 species, 91 in >=3, 18 spanning >=2 phyla
```

---

## What changing the threshold would do

Not measured. If 0.50 / 0.80 is ever changed, the bag corpus is not comparable
across the change and the cluster count is the first thing that moves. A
sensitivity sweep over `--min-seq-id` has not been run, so no claim is made that
0.50 is optimal for anything — only that it is what these bags were built with.

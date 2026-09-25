# Code review of 0262ce3 — what was verified, what was fixed

Five findings were reported against the first commit. Each was reproduced
against real output before anything was changed, because this project has a
history of both real bugs and wrong corrections. Two of my own first readings
here were wrong and are marked as such.

| # | finding | verdict | scale |
|---|---|---|---|
| 1 | tolerant offset used to slice the long allele | **CONFIRMED** | 22,972 of 199,766 catalogue rows (11.5%); 31.1% of tolerant-path rows |
| 2 | dedup counts one event twice | **CONFIRMED** | 351 of 21,051 E. coli events (1.67%) |
| 3 | `Snakefile` broken on fresh checkout | **CONFIRMED** | cannot build a DAG at all |
| 4 | M2 keys on an orientation-dependent hash | **CONFIRMED** | by the repo's own comment in `100_event_dedup.py` |
| 5 | empty catalogue divides by zero | **CONFIRMED** | 2 call sites, before any output is written |

---

## 1. The insert sequence — the one that shipped bad data

`lcp_bp` is an offset in the **short** allele. That is its correct meaning: it
is the insertion point in the empty target, the offset that makes

```
target_seq[:offset] + insert_seq + target_seq[offset:] == derived allele
```

true. But `86_catalogue.py`, `99_mobility_join.py` and `100_event_dedup.py`
each independently did `long_allele[lcp : lcp + inserted_len]`.

On the **exact** path the two frames coincide — a common prefix counts the same
bases in both alleles — so the bug was invisible there. On the **tolerant**
path `tolerant_decompose()` returns `"lcp": sp` (short frame) while the
sequence it hashed into `inserted_md5` was `long_[bl0:bl1]` (long frame). `sp`
and `bl0` differ by exactly the upstream indels.

**Synthetic dose–response** (`121_offset_frame_repro.py`) — a *d* bp deletion
upstream of the junction gives shift −*d*, a *d* bp insertion gives +*d*, and
`ok=True` throughout, so the quality gate does not catch it:

```
control (no upstream indel)   AGREE
del  3 bp   reslice_agrees=False  shift=-3    ins 3 bp   shift=+3
del 12 bp   reslice_agrees=False  shift=-12   ins 8 bp   shift=+8
```

**Real data**, 19 species, `md5(insert_seq)` vs the row's own `insert_md5`:

| method | rows | wrong | |
|---|---:|---:|---:|
| exact | 125,947 | 0 | 0.00% |
| tolerant_alignment | 73,819 | 22,972 | **31.12%** |

The exact path being 0 is the control that makes this attributable: if the
fault were anywhere but the frame, exact rows would be wrong too. Worst
species: H. pylori 52.7%, L. monocytogenes 23.9%, C. jejuni 22.1%.

### Repaired, not re-run

`inserted_md5` was hashed from the **correct** sequence, so it is ground truth
even where the re-slice is wrong. `lib_insert.extract_insert()` searches
offsets nearest-first for the one reproducing it. On 41,607 E. coli loci:

```
re-slice wrong    3,676
repaired          3,512 with a +/-300 bp window, all of them once the search
                  covers the whole allele
ambiguous         0      <- no locus had two offsets matching
```

Zero ambiguity is what makes the repair safe: the md5 anchor is decisive and
"nearest" never has to break a real tie. So Layer 2 does **not** have to be
re-run across 19 species. Post-fix runs of `70_` also write
`insert_start_in_long_bp` and skip the search entirely.

---

## 2. Duplicate events — and a wrong first answer of mine

My first probe reported 99.99% revcomp invariance and I wrote that finding #2
did not reproduce. **That was wrong.** The probe mirrored coordinates
arithmetically, which holds the junction choice fixed — and the junction choice
is the whole mechanism. It could not have found the bug it was testing for.

Re-testing by actually reverse-complementing the alleles and re-deciding the
junction:

| junction overlap | loci | invariant |
|---|---:|---:|
| 0 | 158 | 79.7% |
| 1–3 | 575 | **0.0%** |
| 4–15 | 2,162 | **0.0%** |
| >15 | 66 | **0.0%** |

Any microhomology at the junction — which is every locus carrying a TSD —
and the key depends on which strand the assembly happened to be deposited on.

**Consequence in shipped output:** holding the context key fixed, 351 of 21,051
E. coli events (1.67%) collapse under a junction-jitter-immune key. These are
under-merges by dedup's own definition of an event, not a judgement call.

Fixed losslessly rather than by trimming: `canonical_insert_key()` slides the
insert across the whole direct repeat and takes the smallest canonical key over
the equivalent placements. Trimming the ends would also work but would discard
the ability to separate elements differing only at their termini.

---

## 3. The Snakefile

Two rules call scripts that moved to `scripts/_retired/` during the
simplification work, and `rule annotate` requires
`stage40/mmseqs_rep_seq.fasta`, which **no rule produces** — so the DAG cannot
be built at all, independent of the missing scripts. Retired to
`scripts/_retired/Snakefile.legacy`; README now states the real entry points.

---

## 4 and 5

`inserted_md5` hashes the insert as stored, so the same element on opposite
strands hashes to two values and M2 answers False for both. `100_event_dedup.py`
already said so in a comment; `99_mobility_join.py` used it anyway. M2 now uses
`canonical_insert_key`, which fixes orientation and junction jitter together.

The two `100.0 * n / len(dec)` calls in `86_catalogue.py` run before anything is
written, so a species with nothing decomposable lost its (legitimately empty)
outputs. M. pneumoniae produced 11 events from 385 genomes, so this is a real
outcome and not an error state.

---

## What the review also turned up: the frozen gates had been failing

Re-running Layer 2 on the three bench clusters failed 8 of the then-16 gates —
**on pre-edit code as well**, byte for byte. Cause, established by sweeping
`--min-insert`: at `--min-insert 0` seven values return exactly
(85/21/9/57/58, `resolved_loci` 234, `agreement_regressions` 0) and the other
two land one locus away. The baseline was taken from a run in which the 500 bp
floor was not in force on both paths — the floor was added the same day, and
the harness was never re-run after it.

The eleven loci that move are all sub-500 bp (1, 220, 245, 300, 300, 300, 328,
390, 392, 403, 453). The floor is working; the baseline was stale. Re-frozen
against `--min-insert 500` on post-fix code, with the old values kept in the
comment.

**Today's change is provably not the cause:** pre-edit and post-edit `70_`
differ in **0 of 11,700 cells** on those clusters. The only column that changed
is the appended `insert_start_in_long_bp`.

### Gate 10

Gate 9 exists because the eight gates before it were bit-identical across a fix
that changed 99 loci. Gate 10 exists for the same reason one step further out:
every field gate 9 reads was **correct** here. The corrupted value was produced
downstream, by re-deriving a sequence `70_` had already computed. So gate 10
asserts the re-derivation agrees:

```
post-fix 70_   insert_frame: {coord: 124, lcp: 0, repaired: 0, unresolved: 0}   PASS
pre-fix  70_   insert_frame: {coord:   0, lcp: 111, repaired: 13, unresolved: 0}  FAIL
```

It is deliberately **not** a count of decomposed loci — that is exactly the
quantity that drifted — so it cannot go stale the way gate 9 did.

---

## Standing lesson

A gate that is never re-run is not a gate. Both the stale baseline and finding
#1 survived because nothing re-derived a value it already had and compared the
two. That comparison is now gate 10, and `lib_insert` refuses to return a
sequence it cannot verify against `inserted_md5` — an unresolved insert is
dropped, never exported, because a wrong insert contaminates CDS calling,
clustering and every downstream bag while looking perfectly well-formed.

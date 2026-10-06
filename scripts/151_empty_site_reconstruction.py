#!/usr/bin/env python3
"""RETRACTED 2026-10-02 -- THIS TEST IS CIRCULAR. Superseded by 152_.

It built the derived allele as `target[:off] + insert + target[off:]` and then
removed `insert` from it. Concatenation and removal are inverses, so its
100.0%% exact result was guaranteed before any data was read. Do not quote it.

Original docstring follows.

Can the empty site be RECONSTRUCTED from a carrier, without observing it?

A copy found by homology tells us an element sits somewhere. It does not give
the pre-insertion sequence, so the copy cannot become a bag site. But the
pipeline's defining identity runs both ways:

    target_seq[:off] + insert_seq + target_seq[off:]  ==  derived allele

so excising the element from a carrier should rejoin the flanks into the empty
site. Where the junction carries a TSD the inserted unit already includes one
repeat copy, so removing it leaves exactly one behind -- which is what the
pre-insertion site had. Nothing is inferred, PROVIDED no target bases were lost.

THE TEST. Every decomposed event already has BOTH alleles, so reconstruction
can be scored against ground truth rather than argued:

  1 take the derived (long) allele and the insert, as a copy-finder would see
    them -- locate the element by ALIGNMENT, not by the stored offset, because
    a copy-finder only ever has homology to go on
  2 excise it and rejoin the flanks
  3 compare to the OBSERVED target_seq

Scored in two ways, because they answer different questions:

  exact        reconstruction == observed empty allele, byte for byte
  flank120     the 120 bp window a bag actually uses is identical -- a
               mismatch far from the junction does not affect the corpus

Stratified by target_bases_lost and by junction overlap, since those are the
two conditions the method's validity rests on: bases lost means the empty site
is NOT recoverable by excision, and a long overlap means the boundary is
ambiguous by exactly that many bases.
"""
import collections
import csv
import sys as _sys
print('RETRACTED: 151_ is circular, use 152_excision_validation.py', file=_sys.stderr)
_sys.exit(2)
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)
CAT = sys.argv[1]
HALF = 60
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 200000


def find_excise(derived, insert):
    """Locate the insert in the derived allele by exact match, then excise.

    A copy-finder has only homology, so the stored offset is deliberately NOT
    used. Exact match is the conservative case: if even an exact-present
    element cannot be excised correctly, an inexact homology hit certainly
    cannot, so this is an upper bound on the method.
    """
    i = derived.find(insert)
    if i < 0:
        return None, "insert_not_found_in_derived"
    if derived.find(insert, i + 1) >= 0:
        return derived[:i] + derived[i + len(insert):], "ambiguous_multiple_hits"
    return derived[:i] + derived[i + len(insert):], "ok"


tot = collections.Counter()
by_lost = collections.defaultdict(collections.Counter)
by_ov = collections.defaultdict(collections.Counter)

for p in sorted(glob.glob(os.path.join(CAT, "*_loci.tsv"))):
    for r in csv.DictReader(open(p), delimiter="\t"):
        if tot["seen"] >= LIMIT:
            break
        tgt, ins = r.get("target_seq", ""), r.get("insert_seq", "")
        if not tgt or not ins:
            continue
        try:
            off = int(r["insertion_point_offset"])
            lost = int(r.get("target_bases_lost") or 0)
            ov = int(r.get("junction_overlap_bp") or 0)
        except (ValueError, KeyError):
            continue
        derived = tgt[:off] + ins + tgt[off:]
        tot["seen"] += 1
        rec, status = find_excise(derived, ins)
        lb = "lost=0" if lost == 0 else ("lost=1-10" if lost <= 10 else "lost>10")
        ob = ("ov=0" if ov == 0 else "ov=1-3" if ov <= 3
              else "ov=4-15" if ov <= 15 else "ov>15")
        if rec is None:
            tot[status] += 1
            by_lost[lb][status] += 1
            by_ov[ob][status] += 1
            continue
        exact = (rec == tgt)
        # the window a bag actually uses
        a = tgt[off - HALF:off + HALF] if off >= HALF and len(tgt) - off >= HALF else None
        b = rec[off - HALF:off + HALF] if off >= HALF and len(rec) - off >= HALF else None
        f120 = (a is not None and a == b)
        for d, k in ((tot, None), (by_lost[lb], None), (by_ov[ob], None)):
            d["scored"] += 1
            d["exact"] += exact
            d["flank120"] += f120
            if status == "ambiguous_multiple_hits":
                d["ambiguous"] += 1

print("events scored: %d of %d seen" % (tot["scored"], tot["seen"]))
print("  insert not found in derived : %d" % tot["insert_not_found_in_derived"])
print()
print("%-12s %9s %10s %10s %11s" % ("", "n", "exact", "flank120", "ambiguous"))


def row(lab, c):
    n = max(1, c["scored"])
    print("%-12s %9d %9.1f%% %9.1f%% %10d"
          % (lab, c["scored"], 100.0 * c["exact"] / n,
             100.0 * c["flank120"] / n, c["ambiguous"]))


row("ALL", tot)
print()
print("by target_bases_lost -- excision CANNOT be right when bases were lost")
for k in ("lost=0", "lost=1-10", "lost>10"):
    if by_lost[k]["scored"]:
        row("  " + k, by_lost[k])
print()
print("by junction overlap -- a TSD makes the boundary ambiguous by its length")
for k in ("ov=0", "ov=1-3", "ov=4-15", "ov>15"):
    if by_ov[k]["scored"]:
        row("  " + k, by_ov[k])
print()
print("flank120 is the operative number: it is the window a bag site uses.")

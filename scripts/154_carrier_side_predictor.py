#!/usr/bin/env python3
"""Can excision correctness be predicted from the CARRIER ALONE?

153_ found that every feature it tested requires both alleles -- they are all
products of the decomposition -- and the one rule that reached 100% (the exact
decomposition path) is tautological: `decompose()` calls a locus exact exactly
when the short allele IS the long allele minus one contiguous block, so
excising that block returns it by definition.

The real setting has no second allele. A copy-finder sees one genome holding
an element and must decide, from that alone, whether excising it yields a
trustworthy pre-insertion site. So the feature set here is restricted to what
is visible WITHOUT the empty allele:

  tsd_len        a direct repeat immediately flanking the element in the
                 carrier, found by extending outward from both ends while the
                 bases agree. Detectable from one genome -- the same signal
                 structural polarity already uses at 88.9% concordance. This
                 is the load-bearing candidate: excision is only well defined
                 if the boundary is, and a clean TSD says where the boundary is.
  end_clean      the element's terminal 20 bp match the reference element's
                 terminal 20 bp -- a ragged end means an uncertain cut point
  dist_to_edge   bases from the element to the nearest contig end
  n_hits         times the element occurs in this carrier

GROUND TRUTH is still the observed empty allele, and the set is still selected
for sites demonstrably once empty, so any precision here remains an UPPER
BOUND on the population reconstruction is actually for.

Split-sample by panel hash, as before: derive proposes, TEST decides.
"""
import collections
import csv
import glob
import hashlib
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI

HALF = 60
dirs = []
for pat in sys.argv[1].split(","):
    dirs.extend(sorted(glob.glob(pat)))


def tsd_len(seq, a, b, cap=40):
    """Direct repeat flanking [a,b) in `seq`, measured from the carrier only.

    Walks outward: seq[a-1-k] vs seq[b-1-k] -- the base just left of the
    element against the last base of the element. A TSD of length d makes the
    first d comparisons agree.
    """
    n = 0
    while n < cap and a - 1 - n >= 0 and b - 1 - n >= 0:
        if seq[a - 1 - n] != seq[b - 1 - n]:
            break
        n += 1
    return n


rows = []
for d in dirs:
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)):
        continue
    panel = LI.panel_label(d)
    half = int(hashlib.md5(panel.encode()).hexdigest()[:8], 16) % 2
    aseq = LI.read_alleles(ap)
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        av = aseq.get(r["locus_id"], {})
        lng, sht = av.get(r["longest_allele"], ""), av.get(r["shortest_allele"], "")
        if not lng or not sht:
            continue
        i0, _ = LI.locate_insert(r, lng)
        if i0 < 0:
            continue
        ilen = int(r["inserted_len"])
        ins = lng[i0:i0 + ilen]
        j = lng.find(ins)
        if j < 0:
            continue
        nhits = lng.count(ins)
        rec = lng[:j] + lng[j + ilen:]
        if j < HALF or len(sht) - j < HALF or len(rec) - j < HALF:
            continue
        ok = sht[j - HALF:j + HALF] == rec[j - HALF:j + HALF]
        rows.append({
            "half": half, "ok": ok,
            "tsd": tsd_len(lng, j, j + ilen),
            "dist_edge": min(j, len(lng) - (j + ilen)),
            "nhits": nhits,
            "ilen": ilen})

print("loci scored %d   derive %d / test %d"
      % (len(rows), sum(1 for x in rows if x["half"] == 0),
         sum(1 for x in rows if x["half"] == 1)))
print("baseline flank120 accuracy: %.1f%%"
      % (100.0 * sum(1 for x in rows if x["ok"]) / max(1, len(rows))))
print()


def ev(rule, label):
    out = []
    for h in (0, 1):
        sub = [x for x in rows if x["half"] == h]
        kept = [x for x in sub if rule(x)]
        if not kept:
            out.append("kept 0")
            continue
        out.append("prec %5.1f%%  cov %5.1f%%" %
                   (100.0 * sum(1 for x in kept if x["ok"]) / len(kept),
                    100.0 * len(kept) / max(1, len(sub))))
    print("  %-40s derive %-26s TEST %s" % (label, out[0], out[1]))


print("carrier-side TSD length (the load-bearing candidate)")
for lo, hi in ((0, 0), (1, 3), (4, 8), (9, 15), (16, 40)):
    ev(lambda x, a=lo, b=hi: a <= x["tsd"] <= b, "tsd %d-%d" % (lo, hi))
print()
print("other carrier-side features")
ev(lambda x: x["nhits"] == 1, "element occurs once in carrier")
ev(lambda x: x["dist_edge"] >= 500, "at least 500bp from contig edge")
ev(lambda x: x["dist_edge"] >= 2000, "at least 2000bp from contig edge")
ev(lambda x: x["ilen"] >= 1000, "element >= 1kb")
print()
print("combinations")
ev(lambda x: x["tsd"] >= 4 and x["nhits"] == 1, "tsd>=4 AND single hit")
ev(lambda x: 4 <= x["tsd"] <= 15 and x["nhits"] == 1, "tsd 4-15 AND single hit")
ev(lambda x: 4 <= x["tsd"] <= 15 and x["nhits"] == 1 and x["dist_edge"] >= 500,
   "  + >=500bp from edge")
print()
print("Ground truth is still the observed empty allele, so these are upper")
print("bounds on the population where no empty allele exists.")

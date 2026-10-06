#!/usr/bin/env python3
"""Does excising an element from a REAL carrier reproduce the REAL empty site?

RETRACTS 151_empty_site_reconstruction.py, which was circular. That script
built the derived allele as `target[:off] + insert + target[off:]` and then
removed `insert` from it. Concatenation and removal are inverses, so its
100.0% exact was guaranteed before any data was read and it measured nothing.

This uses the OBSERVED alleles from alleles.fna on both sides:

    long  = the real filled allele as assembled in a carrier genome
    short = the real empty allele as assembled in a different genome

and asks whether locating the element in `long` BY ALIGNMENT and excising it
reproduces `short`. That is not an identity. The two alleles come from
different genomes, so they differ by whatever SNPs and indels separate those
strains, and the tolerant decomposition path exists precisely because they do.
A mismatch is therefore informative rather than definitional.

WHY ALIGNMENT AND NOT THE STORED OFFSET. A copy-finder has only homology: it
sees an element somewhere in a contig and must decide where its boundaries
are. Using the recorded lcp would hand the answer over. So the element is
re-located here the way a copy-finder would have to.

WHAT THIS STILL CANNOT TEST, and no amount of data on disk can:

  * whether the site was EVER empty. Excision yields the sequence that would
    be there if the element left precisely. For an element that is ancestral
    or fixed, that sequence may never have existed in any organism. Every
    locus scored here HAS an observed empty allele, so the set is selected for
    sites that were demonstrably once empty -- the very population where
    reconstruction is easiest. Accuracy here is an UPPER BOUND on accuracy
    where no empty allele is known, which is the case reconstruction is for.
"""
import collections
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI

RECON = sys.argv[1]
HALF = 60
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 60000

tot = collections.Counter()
by_lost = collections.defaultdict(collections.Counter)
by_ov = collections.defaultdict(collections.Counter)
dist = []

dirs = []
for pat in RECON.split(","):
    dirs.extend(sorted(glob.glob(pat)))

for d in dirs:
    if tot["scored"] >= LIMIT:
        break
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)):
        continue
    aseq = LI.read_alleles(ap)
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if tot["scored"] >= LIMIT:
            break
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        av = aseq.get(r["locus_id"], {})
        lng = av.get(r["longest_allele"], "")
        sht = av.get(r["shortest_allele"], "")
        if not lng or not sht:
            continue
        i0, st = LI.locate_insert(r, lng)
        if i0 < 0:
            tot["insert_unresolved"] += 1
            continue
        ilen = int(r["inserted_len"])
        ins = lng[i0:i0 + ilen]
        tot["seen"] += 1
        # how a copy-finder would do it: find the element by homology, excise,
        # rejoin the flanks. Exact search is the generous case.
        j = lng.find(ins)
        if j < 0:
            tot["not_found"] += 1
            continue
        if lng.find(ins, j + 1) >= 0:
            tot["ambiguous"] += 1
        rec = lng[:j] + lng[j + ilen:]
        try:
            lost = int(r.get("target_bases_lost") or 0)
            ov = int(r.get("junction_ambiguity_bp") or 0)
        except ValueError:
            continue
        exact = (rec == sht)
        # FRAME. `j` is a LONG-allele coordinate. The junction in the SHORT
        # allele sits at lcp_bp, and on the tolerant path the two differ by the
        # upstream indels -- the same long-vs-short confusion as the catalogue
        # bug. Comparing both windows at `j` misaligns them and reports a
        # mismatch rate that is an artefact of the index, not of the excision.
        # `rec` is the long allele minus the element, so ITS junction is at j.
        try:
            soff = int(r["lcp_bp"])
        except (KeyError, ValueError):
            continue
        a = sht[soff - HALF:soff + HALF] if soff >= HALF and len(sht) - soff >= HALF else None
        b = rec[j - HALF:j + HALF] if j >= HALF and len(rec) - j >= HALF else None
        f120 = (a is not None and a == b)
        if a is not None and b is not None and not f120:
            dist.append(sum(1 for x, y in zip(a, b) if x != y))
        lb = "lost=0" if lost == 0 else ("lost=1-10" if lost <= 10 else "lost>10")
        ob = ("ov=0" if ov == 0 else "ov=1-3" if ov <= 3
              else "ov=4-15" if ov <= 15 else "ov>15")
        for c in (tot, by_lost[lb], by_ov[ob]):
            c["scored"] += 1
            c["exact"] += exact
            c["f120"] += f120
            c["f120_eval"] += (a is not None and b is not None)

print("scored %d  (insert unresolved %d, not found %d, multi-hit %d)"
      % (tot["scored"], tot["insert_unresolved"], tot["not_found"], tot["ambiguous"]))
print()
print("%-12s %9s %11s %11s" % ("", "n", "exact", "flank120"))


def row(lab, c):
    n = max(1, c["scored"])
    e = max(1, c["f120_eval"])
    print("%-12s %9d %10.1f%% %10.1f%%"
          % (lab, c["scored"], 100.0 * c["exact"] / n, 100.0 * c["f120"] / e))


row("ALL", tot)
print()
for k in ("lost=0", "lost=1-10", "lost>10"):
    if by_lost[k]["scored"]:
        row("  " + k, by_lost[k])
print()
for k in ("ov=0", "ov=1-3", "ov=4-15", "ov>15"):
    if by_ov[k]["scored"]:
        row("  " + k, by_ov[k])
if dist:
    dist.sort()
    print()
    print("when the 120bp window differs: median %d mismatched bases, p90 %d, max %d"
          % (dist[len(dist) // 2], dist[int(len(dist) * .9)], dist[-1]))
print()
print("UPPER BOUND: every locus here HAS an observed empty allele, so the set")
print("is selected for sites demonstrably once empty -- the easy case.")

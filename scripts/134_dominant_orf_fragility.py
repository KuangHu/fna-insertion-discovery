#!/usr/bin/env python3
"""When an insert has two ORFs, how safe is "take the longest"?

Bag identity comes from ONE ORF per insert, chosen by nucleotide length with no
tie-break: `lst.sort(key=lambda x: -x[3]); dom = lst[0]`. 29.3% of inserts carry
exactly two ORFs, so for those the bag rests on which of two genes happened to
be longer.

That is safe when the winner is clearly longer and fragile when the two are
close. A near-tie means a few codons of prodigal boundary difference would send
the insert into a different bag entirely -- and since `sort` is stable, today's
winner on an exact tie is just whichever prodigal numbered first, i.e. the
leftmost.

  ratio = len(second ORF) / len(longest ORF), per insert, nucleotide span

  ratio > 0.9   near-tie: the choice is close to arbitrary
  ratio > 0.8   fragile
  ratio < 0.5   the dominant ORF is at least twice the other -- a clear winner

This measures FRAGILITY, not error. Nothing here says which ORF is the
transposase; no annotation is used anywhere in this pipeline. It says how much
of the corpus has a bag assignment that a small change in ORF calling could
flip.
"""
import collections
import os
import sys

gff = sys.argv[1]

# prodigal GFF: seqid is the insert record; each CDS line is one ORF
lens = collections.defaultdict(list)
for line in open(gff):
    if not line or line[0] == "#":
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 5 or f[2] != "CDS":
        continue
    try:
        st, en = int(f[3]), int(f[4])
    except ValueError:
        continue
    lens[f[0]].append(en - st + 1)

hist = collections.Counter()
ratios = []
for rec, L in lens.items():
    hist[len(L)] += 1
    if len(L) == 2:
        L.sort(reverse=True)
        ratios.append(L[1] / float(L[0]))

print("inserts with >=1 ORF : %d" % len(lens))
print("inserts with exactly 2 ORFs : %d" % hist[2])
print()
if not ratios:
    print("no two-ORF inserts found")
    sys.exit(0)
ratios.sort()
n = len(ratios)


def q(p):
    return ratios[min(n - 1, int(n * p))]


print("second/longest ORF length ratio, two-ORF inserts")
print("  p10 %.3f   p25 %.3f   median %.3f   p75 %.3f   p90 %.3f"
      % (q(.10), q(.25), q(.50), q(.75), q(.90)))
print()
for thr, lab in ((0.9, "near-tie   ratio > 0.90"),
                 (0.8, "fragile    ratio > 0.80"),
                 (0.5, "clear win  ratio < 0.50")):
    if thr == 0.5:
        c = sum(1 for r in ratios if r < thr)
    else:
        c = sum(1 for r in ratios if r > thr)
    print("  %-26s %7d  (%.1f%% of two-ORF inserts, %.1f%% of all)"
          % (lab, c, 100.0 * c / n, 100.0 * c / max(1, len(lens))))
print()
print("Fragility, not error: nothing here identifies which ORF is the")
print("transposase. It bounds how much of the corpus could change bag on a")
print("small difference in ORF boundaries.")

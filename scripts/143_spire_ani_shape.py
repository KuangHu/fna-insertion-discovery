#!/usr/bin/env python3
"""What shape is a MAG species' ANI distribution, and where could a panel
threshold sit?

The project's standing finding is that an ANI threshold does NOT transfer: 99.0
is E. coli's q=0.927 but K. pneumoniae's q=0.386, and the right procedure is to
read the FLOOR and the GAP off each species' own distribution rather than
carrying a quantile across. Those were two well-sampled isolate species. These
are uncultured MAGs inside a 95% ANI cluster boundary, so there is even less
reason to expect 99.0 to mean anything here.

This prints the distribution rather than a verdict:

  ANI quantiles           where the mass actually sits
  AF quantiles            alignment fraction is the parameter most likely to
                          behave differently for MAGs -- an incomplete genome
                          has genuinely lower AF against a complete one for
                          reasons unrelated to relatedness
  histogram               a mode at the cluster boundary (95%) with nothing
                          above it means no fine substructure and no panels
  neighbours per genome   at a ladder of candidate thresholds, with the MEDIAN
                          and the share of genomes having >= 11 neighbours,
                          since a k=12 clique needs 11 partners as a necessary
                          (not sufficient) condition

The last column is the one that matters. It is an UPPER BOUND on panel
formability: having 11 neighbours does not make them mutually comparable, which
is why single-linkage degree was already shown to be a bad predictor -- only
~10% of a neighbourhood is mutually comparable. If a species cannot even clear
the upper bound, it certainly cannot form panels.
"""
import collections
import sys

path, label = sys.argv[1], sys.argv[2]
AF_MIN = float(sys.argv[3]) if len(sys.argv) > 3 else 85.0

ani, af = [], []
deg = collections.defaultdict(lambda: collections.Counter())
THR = [95.0, 97.0, 98.0, 98.5, 99.0, 99.5, 99.9]
genomes = set()

with open(path) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    cA = ix.get("ANI", 2)
    cQ, cR = ix.get("Query_file", 0), ix.get("Ref_file", 1)
    cAF = ix.get("Align_fraction_query", 3)
    for line in fh:
        r = line.rstrip("\n").split("\t")
        if len(r) <= max(cA, cQ, cR, cAF):
            continue
        q, t = r[cQ], r[cR]
        if q == t:
            continue
        try:
            a = float(r[cA])
            fq = float(r[cAF])
        except ValueError:
            continue
        genomes.add(q)
        ani.append(a)
        af.append(fq)
        if fq < AF_MIN:
            continue
        for th in THR:
            if a >= th:
                deg[q][th] += 1

if not ani:
    print("%s: no usable rows" % label)
    sys.exit(0)
ani.sort()
af.sort()
n = len(ani)


def q(v, p):
    return v[min(len(v) - 1, int(len(v) * p))]


print("=" * 72)
print("%s   %d genomes, %d non-self pairs" % (label, len(genomes), n))
print("=" * 72)
print("ANI  p05 %.2f  p25 %.2f  median %.2f  p75 %.2f  p95 %.2f  max %.2f"
      % (q(ani, .05), q(ani, .25), q(ani, .50), q(ani, .75), q(ani, .95), ani[-1]))
print("AF   p05 %.1f  p25 %.1f  median %.1f  p75 %.1f  p95 %.1f"
      % (q(af, .05), q(af, .25), q(af, .50), q(af, .75), q(af, .95)))
print()
print("ANI histogram (all non-self pairs)")
h = collections.Counter()
for a in ani:
    h[min(100, int(a))] += 1
for b in sorted(h):
    if h[b]:
        print("  %3d-%3d  %9d  %s" % (b, b + 1, h[b], "#" * min(44, h[b] * 44 // n + 1)))
print()
print("neighbours per genome at AF >= %.0f  (a k=12 clique needs >= 11)" % AF_MIN)
print("  %-8s %10s %10s %12s" % ("ANI", "median", "p90", "genomes>=11"))
ng = len(genomes)
for th in THR:
    d = sorted(deg[g][th] for g in genomes)
    if not d:
        d = [0]
    ok = sum(1 for x in d if x >= 11)
    print("  %-8.1f %10d %10d %7d (%4.1f%%)"
          % (th, d[len(d) // 2], d[int(len(d) * .9)], ok, 100.0 * ok / max(1, ng)))
print()
print("UPPER BOUND only: 11 neighbours does not make them mutually comparable.")

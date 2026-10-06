#!/usr/bin/env python3
"""How many proteins does an insert actually encode?

The bag's identity comes from ONE ORF -- the longest -- on the reasoning that
where an element carries a transposase it is normally the longest ORF and
passenger genes are shorter. That heuristic is safe when inserts carry one or
two ORFs and increasingly shaky as the count rises: with five ORFs, "the
longest" is a much weaker claim about what the element IS.

The corpus mean is 1.75 ORFs per insert, but a mean cannot answer this. What
matters is the distribution, and specifically the share at n_orfs <= 2.

Also reported, because it is the instrument that shows where the heuristic is
weak:

  dominant_orf_frac   how much of the insert the chosen ORF covers. A dominant
                      ORF spanning 90% of the element is a strong identity; one
                      spanning 25% is a label taken from a fragment of a
                      multi-gene element.

Broken out per species AND for the sites that actually reach a bag, since an
insert with no ORF never becomes a bag at all and would otherwise dilute the
denominator.
"""
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)

path = sys.argv[1]
hist = collections.Counter()
per_sp = collections.defaultdict(collections.Counter)
frac_by_n = collections.defaultdict(list)
n = 0
for r in csv.DictReader(open(path), delimiter="\t"):
    try:
        k = int(r["n_orfs"])
    except (ValueError, KeyError):
        continue
    n += 1
    hist[k] += 1
    per_sp[r.get("species", "?")][k] += 1
    try:
        frac_by_n[min(k, 6)].append(float(r["dominant_orf_frac"]))
    except (ValueError, KeyError):
        pass

print("inserts with a CDS record: %d" % n)
print()
print("%-10s %10s %8s %10s" % ("n_orfs", "inserts", "pct", "cumulative"))
cum = 0
for k in sorted(hist):
    cum += hist[k]
    print("%-10d %10d %7.1f%% %9.1f%%"
          % (k, hist[k], 100.0 * hist[k] / n, 100.0 * cum / n))
le2 = sum(v for k, v in hist.items() if k <= 2)
print()
print("n_orfs <= 2 : %d of %d  (%.1f%%)" % (le2, n, 100.0 * le2 / n))
print("n_orfs >= 4 : %d of %d  (%.1f%%)"
      % (sum(v for k, v in hist.items() if k >= 4), n,
         100.0 * sum(v for k, v in hist.items() if k >= 4) / n))
print()
print("dominant_orf_frac by ORF count -- how much of the insert the chosen ORF covers")
print("%-10s %8s %9s %9s %9s" % ("n_orfs", "n", "median", "p25", "p10"))
for k in sorted(frac_by_n):
    v = sorted(frac_by_n[k])
    if not v:
        continue
    lab = "%d" % k if k < 6 else "6+"
    print("%-10s %8d %9.3f %9.3f %9.3f"
          % (lab, len(v), v[len(v) // 2], v[len(v) // 4], v[len(v) // 10]))
print()
print("%-16s %8s %9s %9s %9s" % ("species", "inserts", "n=1", "n=2", "<=2 pct"))
for sp in sorted(per_sp, key=lambda x: -sum(per_sp[x].values())):
    c = per_sp[sp]
    t = sum(c.values())
    print("%-16s %8d %9d %9d %8.1f%%"
          % (sp, t, c[1], c[2], 100.0 * (c[1] + c[2]) / max(1, t)))

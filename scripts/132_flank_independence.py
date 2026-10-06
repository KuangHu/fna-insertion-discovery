#!/usr/bin/env python3
"""Are the "sites" in a bag independent, or the same site more than once?

The step-5 metric is "bags reaching >=5 / >=10 INDEPENDENT sites". 112_ emits
one site per database row and never dedups on the flank, so two events that
share a target site but carry different inserts each become a site. They are
two alleles of one locus, not two independent observations, and counting both
inflates the metric that decides whether a bag is usable.

This measures the gap on the corpus as built, before anything is added to it:

  sites            rows in insertions_sites.jsonl
  distinct flanks  distinct canonical 120 bp windows
  bags >=5 / >=10  counted BOTH ways -- by raw site rows, and by distinct flanks

If the two counts agree the concern is theoretical. If they diverge, the
published bag counts are overstated and the fix belongs in 112_, not in a
caveat.
"""
import collections
import json
import sys

path = sys.argv[1]
by_bag_sites = collections.Counter()
by_bag_flanks = collections.defaultdict(set)
n = 0
dup_flank = collections.Counter()
for line in open(path):
    line = line.strip()
    if not line:
        continue
    try:
        d = json.loads(line)
    except ValueError:
        continue
    n += 1
    bag = d.get("bag_id") or d.get("bag") or "?"
    fl = d.get("flank") or ""
    by_bag_sites[bag] += 1
    by_bag_flanks[bag].add(fl)
    dup_flank[fl] += 1

tot_flanks = len(dup_flank)
print("sites                     %8d" % n)
print("distinct flanks           %8d  (%.1f%% of sites)"
      % (tot_flanks, 100.0 * tot_flanks / max(1, n)))
rep = sum(v - 1 for v in dup_flank.values() if v > 1)
print("repeat site rows          %8d  (same 120bp window seen again)" % rep)
print("flanks seen >=2 times     %8d" % sum(1 for v in dup_flank.values() if v > 1))
print("max times one flank seen  %8d" % (max(dup_flank.values()) if dup_flank else 0))
print()
print("%-26s %10s %10s" % ("", "by rows", "by flanks"))
for k in (2, 5, 10, 20):
    a = sum(1 for b, c in by_bag_sites.items() if c >= k)
    b_ = sum(1 for b, s in by_bag_flanks.items() if len(s) >= k)
    print("%-26s %10d %10d" % ("bags with >=%d sites" % k, a, b_))
print()
print("'by flanks' is the honest denominator for the >=5 / >=10 metric:")
print("two alleles at one locus are one observation of that site.")

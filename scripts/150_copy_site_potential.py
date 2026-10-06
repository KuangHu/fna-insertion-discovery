#!/usr/bin/env python3
"""If every detected COPY became a bag site, how many sites would that be?

Today one event contributes exactly one site. But copy-finding already runs:
M1 counts copies inside the carrier genome and records how many sit at
distinct loci, and M3 counts genomes of the species carrying the insert. Each
of those locations is a place the element sits -- a potential (empty site,
insert) pair -- and none of them reaches the corpus, because a bag site needs
the PRE-INSERTION sequence and copy-finding establishes only presence.

This bounds the prize before anything is built.

  sites now              events that produced a bag site
  copies within genome   sum of insert_distinct_loci - 1 per event
                         (minus 1: the called locus is already a site)
  extra from M3          genomes carrying the insert, beyond the carriers
                         already represented

CEILING, NOT FORECAST. Three deductions apply before any of this is real:
  * a copy is only a site if its empty allele can be recruited or
    reconstructed -- unmeasured here
  * copies of ONE element at many loci are many sites, but copies of one
    element at the SAME locus across genomes are one site; distinct_loci is
    the right field for the first, M3 counts the second and overcounts
  * the bag key is (bag, flank), so copies landing in identical sequence
    context collapse
So read the M1 column as the defensible bound and the M3 column as the
absolute ceiling.
"""
import collections
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)
CAT = sys.argv[1]

print("%-16s %9s %10s %12s %12s %11s"
      % ("species", "events", "sites_now", "M1_extra", "copies_med", "M1_pos"))
T = collections.Counter()
for p in sorted(glob.glob(os.path.join(CAT, "*_events.tsv"))):
    sp = os.path.basename(p).replace("_events.tsv", "")
    n = extra = m1pos = 0
    cops = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        n += 1
        try:
            dl = int(r.get("insert_distinct_loci") or 0)
        except ValueError:
            dl = 0
        try:
            cp = int(r.get("insert_copies_in_genome") or 0)
        except ValueError:
            cp = 0
        if cp > 1:
            cops.append(cp)
        if r.get("M1_within_genome") == "True":
            m1pos += 1
        if dl > 1:
            extra += dl - 1
    cops.sort()
    med = cops[len(cops) // 2] if cops else 0
    T["events"] += n
    T["extra"] += extra
    T["m1pos"] += m1pos
    print("%-16s %9d %10d %12d %12d %11d" % (sp, n, n, extra, med, m1pos))
print("%-16s %9d %10d %12d %12s %11d"
      % ("TOTAL", T["events"], T["events"], T["extra"], "-", T["m1pos"]))
print()
print("potential sites if every distinct-locus copy were recruited: %d -> %d  (%.1fx)"
      % (T["events"], T["events"] + T["extra"],
         (T["events"] + T["extra"]) / max(1, T["events"])))
print()
print("CEILING only. Each copy needs an empty allele before it is a site.")

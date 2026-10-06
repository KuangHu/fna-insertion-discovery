#!/usr/bin/env python3
"""Is a novel multi-allele PAIR a novel SITE, or a second allele at a known one?

129_ stage 3 dedups on (insert_key, context_key). A pair A->B carrying insert X
and the catalogue event A->C carrying X+Y have DIFFERENT insert keys, so both
survive -- yet they describe the SAME target site in the same genomes. Counting
both as discoveries double-counts one composite difference expressed two ways,
which is precisely the thing that must not happen.

The context key is the +/-100bp window around the insertion point in the empty
allele, so it identifies the SITE independently of what was inserted:

  known_site_new_allele   the context key already exists in the catalogue. The
                          site was already found; this is another allele at it.
                          Real information about allele structure, but NOT a
                          new insertion-target pair to add to a site count.

  new_site                the context key does not occur in the catalogue at
                          all. A target site nothing has reported.

Only `new_site` may be added to a site count. `known_site_new_allele` belongs
in an allele-structure table attached to the site it refines.
"""
import collections
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)

W = sys.argv[1]
print("%-16s %10s %10s %12s %10s %8s"
      % ("species", "pairs", "new_site", "known_site", "v3_events", "new%"))
T = collections.Counter()
for p in sorted(glob.glob(os.path.join(W, "multiallele",
                                       "*_empty_target_pairs.tsv"))):
    sp = os.path.basename(p).replace("_empty_target_pairs.tsv", "")
    if sp.startswith("v2_"):
        continue
    ev = os.path.join(W, "full_%s" % sp, "v3", "dedup_loci_to_event.tsv")
    if not os.path.exists(ev):
        continue
    known_ctx = {r["context_key"]
                 for r in csv.DictReader(open(ev), delimiter="\t")}
    n = ns = ks = 0
    for r in csv.DictReader(open(p), delimiter="\t"):
        n += 1
        if r["context_key"] in known_ctx:
            ks += 1
        else:
            ns += 1
    cat = os.path.join(W, "catalogue_v3", "%s_events.tsv" % sp)
    ne = (sum(1 for _ in open(cat)) - 1) if os.path.exists(cat) else 0
    T["pairs"] += n
    T["new"] += ns
    T["known"] += ks
    T["events"] += ne
    print("%-16s %10d %10d %12d %10d %7.1f%%"
          % (sp, n, ns, ks, ne, 100.0 * ns / max(1, ne)))
print("%-16s %10d %10d %12d %10d %7.1f%%"
      % ("TOTAL", T["pairs"], T["new"], T["known"], T["events"],
         100.0 * T["new"] / max(1, T["events"])))
print()
print("new%  = new SITES as a share of the existing event count.")
print("known_site rows are additional alleles at sites already reported and")
print("must not be added to a site count.")

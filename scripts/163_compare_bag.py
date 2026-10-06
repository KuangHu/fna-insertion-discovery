#!/usr/bin/env python3
"""Compare one bag between two corpora, field by field.

v5 carries v4's observed sites through unchanged and adds reconstructed ones
that INHERIT their non-coding payload from the source event's own v4 site. So
for any given source event the nc content should be identical in both, and a
difference in the bag's nc-length distribution should come only from which
events contribute sites -- not from the same event reporting different nc.

This separates those two cases:

  per-site      for sites present in BOTH corpora (same site_id), is the nc
                payload byte-identical? A difference here is a bug.
  per-bag       how the bag's nc-length distribution shifts, and whether the
                shift is explained by new source events entering the bag.
"""
import collections
import json
import sys

V4, V5, BAG = sys.argv[1], sys.argv[2], sys.argv[3]


def load(path, bag):
    out = []
    for line in open(path):
        line = line.strip()
        if not line or bag not in line:
            continue
        d = json.loads(line)
        if d.get("bag_id") != bag and d.get("cds_cluster_id") != bag:
            continue
        out.append(d)
    return out


a = load(V4, BAG)
b = load(V5, BAG)
print("bag %s" % BAG)
print("  v4 sites %d   v5 sites %d" % (len(a), len(b)))

src = collections.Counter(x.get("empty_site_source", "?") for x in b)
print("  v5 by source: %s" % dict(src))


def ncstats(rows, lab):
    L = sorted(int(x.get("nc_total_len") or 0) for x in rows)
    C = sorted(int(x.get("nc_region_count") or 0) for x in rows)
    if not L:
        print("  %-26s no sites" % lab)
        return
    print("  %-26s n=%-7d nc_total_len median %-6d min %-6d max %-6d | "
          "regions median %d"
          % (lab, len(L), L[len(L) // 2], L[0], L[-1], C[len(C) // 2]))


ncstats(a, "v4 (all observed)")
ncstats(b, "v5 (all)")
ncstats([x for x in b if x.get("empty_site_source") == "observed"],
        "v5 observed only")
ncstats([x for x in b if x.get("empty_site_source") == "tsd_reconstructed"],
        "v5 reconstructed only")

# per-site identity for sites in both
ai = {x.get("site_id"): x for x in a}
same = diff = 0
examples = []
for x in b:
    sid = x.get("site_id")
    if sid not in ai:
        continue
    y = ai[sid]
    if (y.get("nc_total_len") == x.get("nc_total_len")
            and y.get("nc_sequence_hash") == x.get("nc_sequence_hash")
            and y.get("noncoding_regions") == x.get("noncoding_regions")):
        same += 1
    else:
        diff += 1
        if len(examples) < 3:
            examples.append((sid, y.get("nc_total_len"), x.get("nc_total_len"),
                             y.get("nc_sequence_hash"), x.get("nc_sequence_hash")))
print()
print("  sites present in BOTH: %d identical nc, %d DIFFERENT" % (same, diff))
for e in examples:
    print("    %s  v4 len=%s hash=%s   v5 len=%s hash=%s"
          % (e[0], e[1], e[3], e[2], e[4]))

# do reconstructed sites match their source event's observed nc?
print()
mismatch = 0
checked = 0
for x in b:
    if x.get("empty_site_source") != "tsd_reconstructed":
        continue
    s = x.get("source_event")
    y = ai.get("%s.site" % s)
    if not y:
        continue
    checked += 1
    if y.get("nc_sequence_hash") != x.get("nc_sequence_hash"):
        mismatch += 1
print("  reconstructed sites whose source has a v4 site: %d checked, "
      "%d nc mismatches" % (checked, mismatch))

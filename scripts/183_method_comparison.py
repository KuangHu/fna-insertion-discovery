#!/usr/bin/env python3
"""Which mining method is worth the compute? One species, one metric, three methods.

Raw output size is not the metric. A method that returns 83,593 rows of the
same element at the same locus in 400 clonal strains has produced one piece of
evidence. What a candidate system needs is the same element at MANY INDEPENDENT
TARGET SITES, because that is what makes guide-target coherence testable.

So every method is scored the same way:

  element families     distinct CDS clusters represented
  independent sites    distinct (genome, contig, position/50bp) -- two copies
                       in one genome at different loci are two backgrounds;
                       the same locus in 400 strains is one
  families with >=2    families having two or more independent sites. This is
                       the number that matters: a family at one site can never
                       support a coherence test regardless of how many rows
                       carry it
  per node-hour        the same counts divided by measured wall-clock x cores

WHY THE SITE KEY IS COARSE (position // 50). Assemblies differ by a few bases
in where they place a junction, so an exact coordinate would count one locus
as several. 50bp is wide enough to absorb that and narrow enough that two real
loci are not merged -- elements are >=500bp apart by construction.

WHAT THIS DOES NOT MEASURE. Whether the sites are correct. The copy-based
method's sites rest on excision accuracy measured at 84.7% on the easiest
population, and the cross-species method's on a homology threshold calibrated
for within-species comparison. A method can win on count and lose on truth;
this ranks yield per compute, not quality.
"""
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)


def load(path, cluster_col, genome_col, contig_col, pos_col, node_hours, label,
         cluster_map=None):
    fam = collections.defaultdict(set)
    n = 0
    try:
        for r in csv.DictReader(open(path), delimiter="\t"):
            n += 1
            c = r.get(cluster_col, "")
            if cluster_map is not None:
                c = cluster_map.get(c, "")
            if not c or c == ".":
                continue
            g = r.get(genome_col, ".")
            ct = r.get(contig_col, ".")
            try:
                p = int(r.get(pos_col) or 0) // 50
            except ValueError:
                p = 0
            fam[c].add((g, ct, p))
    except OSError as e:
        print("  %-26s UNREADABLE: %s" % (label, e))
        return None
    sites = sum(len(v) for v in fam.values())
    multi = sum(1 for v in fam.values() if len(v) >= 2)
    return {"label": label, "rows": n, "families": len(fam), "sites": sites,
            "multi": multi, "nh": node_hours}


W = "/global/scratch/users/kh36969/fna_ins_discovery"
SP = "spneumoniae"

# db_id -> cds_cluster, needed by the methods that key on the event
cmap = {}
for r in csv.DictReader(open(W + "/cds_v3/insert_cds.tsv"), delimiter="\t"):
    cmap[r["db_id"]] = r.get("cds_cluster_id", ".")

res = []

# 1. within-species: the full pipeline. node-hours measured from the species
#    driver -- Arm A + Arm B + Layer 2 for one species on a 64-core node.
# The catalogue carries no per-locus genomic coordinate, so it cannot go
# through load(): an event IS a distinct (insert, context) pair, i.e. one site
# by construction. Computed directly instead -- an earlier version ran it
# through load() with a degenerate site key and reported 0 families with >=2
# sites, which was an artefact of the key, not a property of the method.
fam_w = collections.defaultdict(set)
nrow_w = 0
for x in csv.DictReader(open(W + "/catalogue_v3/%s_events.tsv" % SP),
                        delimiter="\t"):
    nrow_w += 1
    c = cmap.get("%s.%s" % (SP, x["event_id"]), ".")
    if c and c != ".":
        fam_w[c].add(x["event_id"])          # each event is its own site
res.append({"label": "within-species panels", "rows": nrow_w,
            "families": len(fam_w),
            "sites": sum(len(v) for v in fam_w.values()),
            "multi": sum(1 for v in fam_w.values() if len(v) >= 2),
            "nh": 64 * 1.5})

# 2. copy-based reconstruction
res.append(load(W + "/copysites/%s_sites.tsv" % SP, "source_event", "carrier",
                "contig", "elem_start", 64 * 0.40, "copy-based excision",
                cluster_map=cmap))

# 3. cross-species expansion
res.append(load(W + "/xspecies/%s_hits.tsv" % SP, "cds_cluster", "genome",
                "contig", "elem_start", 64 * 0.13, "cross-species expansion"))

res = [r for r in res if r]
print("ONE SPECIES (%s) -- three methods, same metric" % SP)
print("node-hours are ESTIMATES from job wall-clock x cores, not instrumented.")
print()
print("%-26s %9s %9s %9s %9s %10s" %
      ("method", "rows", "families", "sites", "fam>=2", "node-hours"))
for r in res:
    print("%-26s %9d %9d %9d %9d %10.1f"
          % (r["label"], r["rows"], r["families"], r["sites"], r["multi"], r["nh"]))
print()
print("%-26s %14s %14s" % ("per node-hour", "sites", "families>=2"))
for r in res:
    print("%-26s %14.1f %14.2f"
          % (r["label"], r["sites"] / max(0.1, r["nh"]),
             r["multi"] / max(0.1, r["nh"])))
print()
print("`fam>=2` is the one that decides usefulness: a family at a single site")
print("cannot support a coherence test no matter how many rows carry it.")
print("Yield per compute only -- this says nothing about whether the sites")
print("are correct.")

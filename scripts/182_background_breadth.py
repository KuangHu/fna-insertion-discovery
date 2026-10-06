#!/usr/bin/env python3
"""How many INDEPENDENT target backgrounds does each element family occupy?

This is the selection metric: a candidate system is worth pursuing when the
same element appears at many independent target sites, because that is what
makes guide-target coherence testable. Total event count does not measure it --
one element recurring across clonal strains and overlapping panels inflates
events without adding a single independent background.

ONLY CROSS-SPECIES HITS ARE COUNTED AS NEW. A hit in the element's own source
species is largely what Arm A already reports: within-genome multi-copy
detection covers it. A hit in a DIFFERENT species is an independent background
that no amount of within-species mining produces, so same-species hits are
tallied separately and excluded from the breadth score.

INDEPENDENCE IS COUNTED AT THREE LEVELS, because they answer different things:

  species      how many different species carry this family
  genomes      how many different assemblies -- NOT assumed to be distinct
               samples: one BioSample re-deposited under two accessions is one
               background, and this script cannot see BioSample, so the genome
               count is an UPPER bound on independent hosts
  sites        distinct (genome, contig, position/50bp) -- the quantity that
               actually matters, since two copies in one genome at different
               loci are two target backgrounds

ONLY FLANK-COMPLETE LOCI COUNT. A hit running into a contig end has unknown
extent and its flanking sequence is not usable as a target context, so it is
excluded rather than counted with a caveat.
"""
import collections
import csv
import glob
import sys

csv.field_size_limit(sys.maxsize)

rows = []
for p in sorted(glob.glob(sys.argv[1])):
    rows.extend(csv.DictReader(open(p), delimiter="\t"))
if not rows:
    print("no hits found")
    sys.exit(0)

MIN_SP = int(sys.argv[2]) if len(sys.argv) > 2 else 2

fam = collections.defaultdict(lambda: {
    "xsp_species": set(), "xsp_genomes": set(), "xsp_sites": set(),
    "same_sites": set(), "src": None, "trunc": 0})
n_all = n_complete = 0
for r in rows:
    n_all += 1
    if r.get("flank_complete") != "yes":
        fam[r["cds_cluster"]]["trunc"] += 1
        continue
    n_complete += 1
    f = fam[r["cds_cluster"]]
    f["src"] = r.get("source_species")
    site = (r["genome"], r["contig"], int(r["elem_start"]) // 50)
    if r.get("cross_species") == "yes":
        f["xsp_species"].add(r["target_species"])
        f["xsp_genomes"].add(r["genome"])
        f["xsp_sites"].add(site)
    else:
        f["same_sites"].add(site)

print("hits read %d, flank-complete %d (%.1f%%)"
      % (n_all, n_complete, 100.0 * n_complete / max(1, n_all)))
print("element families seen: %d" % len(fam))
print()

breadth = collections.Counter()
for c, f in fam.items():
    breadth[len(f["xsp_species"])] += 1
print("families by number of DIFFERENT species they occupy (cross-species only)")
print("  %-10s %8s" % ("species", "families"))
for k in sorted(breadth):
    print("  %-10d %8d" % (k, breadth[k]))
print()

sel = [(c, f) for c, f in fam.items() if len(f["xsp_species"]) >= MIN_SP]
sel.sort(key=lambda kv: (-len(kv[1]["xsp_species"]), -len(kv[1]["xsp_sites"])))
print("families in >= %d species: %d" % (MIN_SP, len(sel)))
print()
print("  %-12s %-16s %7s %8s %8s %9s" %
      ("family", "source_species", "species", "genomes", "sites", "same_sp"))
for c, f in sel[:30]:
    print("  %-12s %-16s %7d %8d %8d %9d"
          % (c[:12], (f["src"] or ".")[:16], len(f["xsp_species"]),
             len(f["xsp_genomes"]), len(f["xsp_sites"]), len(f["same_sites"])))
print()
tot_sites = sum(len(f["xsp_sites"]) for _, f in sel)
print("cross-species target sites in those families: %d" % tot_sites)
print()
print("`genomes` is an UPPER bound on independent hosts: one BioSample")
print("re-deposited under two accessions reads as two genomes here.")

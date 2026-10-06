#!/usr/bin/env python3
"""Four-tier background counting. Retracts the site counts in 183_.

183_ claimed to count independent target backgrounds and did not. Two errors,
both of which inflated the expansion method:

  * its site key was (genome, contig, position). One ancestral insertion
    inherited by 50 strains gives 50 distinct keys and was counted as 50
    backgrounds. Coordinate dedup removes duplicate records inside one
    assembly; it cannot merge a homologous locus across strains.
  * within-species events are ALREADY collapsed on (canonical insert,
    canonical +/-100bp context), so one ancestral locus across 50 strains is
    ONE event there. Comparing a context-deduplicated count against a
    non-deduplicated one biases the comparison entirely one way.

So the tiers are computed identically for every method:

  T1 hit records     raw alignment rows
  T2 physical sites  hits merged at one position. 25% of expansion rows are
                     several element families matching the SAME position --
                     overlapping or nested hits that were never merged.
  T3 host loci       the same locus across different strains merged by FLANK
                     HOMOLOGY. Two hits whose flanking sequence matches are
                     one ancestral locus, however many genomes carry it.
  T4 informative     host loci that can actually inform a coherence test:
                     flank-complete, and belonging to a family present at >=2
                     DISTINCT host loci. A family at one locus cannot support
                     the test regardless of row count.

SAME-SPECIES LOCI ARE NOT EXCLUDED. 183_ dropped them as "already covered by
Arm A". A new locus in the source species is a real new background, and
conversely a cross-species hit may just be the same host region that travelled
with the element. Cross-species breadth is reported as a SEPARATE column, not
folded into the independence count.

COST IS ALLOCATED CORE-HOURS, not node-hours: 64 cores for 1.5h is 96
core-hours. And the expansion's cost is the INCREMENTAL cost of deepening an
existing catalogue -- it cannot discover a family that panel mining did not
already find, so it is not comparable to from-scratch discovery.
"""
import collections
import csv
import hashlib
import sys

csv.field_size_limit(sys.maxsize)
W = "/global/scratch/users/kh36969/fna_ins_discovery"
SP = "spneumoniae"
FLANKKEY = 60     # bp of each flank used to decide "same host locus"


def flank_key(lf, rf):
    """Identity of a host locus from its flanking sequence.

    Exact on a trimmed window rather than an alignment: strains of one species
    differ by ~1 base per 100, so 60bp either side is short enough that most
    copies of one ancestral locus match exactly, and long enough that two
    unrelated loci do not. This UNDER-merges (a SNP in the window splits one
    locus in two), so the host-locus count is an UPPER bound on independence --
    the conservative direction for a method being argued for.
    """
    a = (lf or "")[-FLANKKEY:].upper()
    b = (rf or "")[:FLANKKEY].upper()
    if len(a) < FLANKKEY or len(b) < FLANKKEY:
        return None
    return hashlib.md5((a + "|" + b).encode()).hexdigest()[:16]


cmap = {r["db_id"]: r.get("cds_cluster_id", ".")
        for r in csv.DictReader(open(W + "/cds_v3/insert_cds.tsv"), delimiter="\t")}

out = []

# ---- cross-species expansion -------------------------------------------
rows = list(csv.DictReader(open(W + "/xspecies/%s_hits.tsv" % SP), delimiter="\t"))
t1 = len(rows)
pos = set()
loci = collections.defaultdict(set)      # family -> {host locus key}
xsp = collections.defaultdict(set)       # family -> {species}
nofl = 0
for r in rows:
    pos.add((r["genome"], r["contig"], int(r["elem_start"]) // 50))
    if r.get("flank_complete") != "yes":
        continue
    k = flank_key(r.get("left_flank"), r.get("right_flank"))
    if k is None:
        nofl += 1
        continue
    loci[r["cds_cluster"]].add(k)
    if r.get("cross_species") == "yes":
        xsp[r["cds_cluster"]].add(r["target_species"])
t3 = sum(len(v) for v in loci.values())
t4fam = {c for c, v in loci.items() if len(v) >= 2}
t4 = sum(len(loci[c]) for c in t4fam)
out.append(("cross-species expansion", t1, len(pos), t3, t4, len(t4fam),
            len(loci), sum(1 for c in t4fam if xsp[c]), 64 * 0.13))

# ---- copy-based excision ------------------------------------------------
rows = list(csv.DictReader(open(W + "/copysites/%s_sites.tsv" % SP), delimiter="\t"))
t1 = len(rows)
pos = set()
loci = collections.defaultdict(set)
for r in rows:
    pos.add((r["carrier"], r["contig"], int(r["elem_start"]) // 50))
    fam = cmap.get(r.get("source_event", ""), ".")
    if fam == ".":
        continue
    # the 120bp flank IS the locus context here; split it to match the
    # left/right convention used above
    fl = r.get("flank_fwd") or r.get("flank") or ""
    if len(fl) != 120:
        continue
    k = flank_key(fl[:60], fl[60:])
    if k:
        loci[fam].add(k)
t3 = sum(len(v) for v in loci.values())
t4fam = {c for c, v in loci.items() if len(v) >= 2}
t4 = sum(len(loci[c]) for c in t4fam)
out.append(("copy-based excision", t1, len(pos), t3, t4, len(t4fam),
            len(loci), 0, 64 * 0.40))

# ---- within-species panels ----------------------------------------------
# Events are already collapsed on (canonical insert, canonical context), so an
# event IS a host locus by construction. Re-keying on flanks here would merge
# nothing and would misrepresent the method.
rows = list(csv.DictReader(open(W + "/catalogue_v3/%s_events.tsv" % SP),
                           delimiter="\t"))
loci = collections.defaultdict(set)
for r in rows:
    fam = cmap.get("%s.%s" % (SP, r["event_id"]), ".")
    if fam != ".":
        loci[fam].add(r["event_id"])
t3 = sum(len(v) for v in loci.values())
t4fam = {c for c, v in loci.items() if len(v) >= 2}
t4 = sum(len(loci[c]) for c in t4fam)
out.append(("within-species panels", len(rows), len(rows), t3, t4, len(t4fam),
            len(loci), 0, 64 * 1.5))

print("ONE SPECIES (%s). Cost is ALLOCATED CORE-HOURS, estimated." % SP)
print()
print("%-26s %8s %8s %8s %9s %8s %8s" %
      ("method", "T1 hits", "T2 pos", "T3 loci", "T4 infor", "families", "fam>=2"))
for lab, t1, t2, t3, t4, nf4, nf, nx, ch in out:
    print("%-26s %8d %8d %8d %9d %8d %8d" % (lab, t1, t2, t3, t4, nf, nf4))
print()
print("%-26s %10s %12s %10s" % ("", "core-hours", "T4/core-hr", "fam>=2/core-hr"))
for lab, t1, t2, t3, t4, nf4, nf, nx, ch in out:
    print("%-26s %10.1f %12.1f %10.2f" % (lab, ch, t4 / ch, nf4 / ch))
print()
print("T3 collapses one ancestral locus across strains by flank identity.")
print("It UNDER-merges (a SNP in the 60bp window splits one locus), so T3 and")
print("T4 are UPPER bounds on independence -- conservative for the method")
print("being argued for.")
print()
print("Expansion cannot discover a family panel mining did not already find;")
print("its cost is incremental deepening, not from-scratch discovery.")

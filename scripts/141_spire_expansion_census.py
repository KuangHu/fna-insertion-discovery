#!/usr/bin/env python3
"""Which species become MINEABLE once SPIRE MAGs are added? Whole-database census.

The earlier profile only asked whether SPIRE tops up the 18 species already
running. That is the narrow question. The real one is whether SPIRE brings
species that were never analysable into range at all -- and answering it needs
a census over every SPIRE cluster, not a lookup of names we already use.

WHY THE JOIN IS ON CLUSTERS, NOT NAMES. GTDB and NCBI disagree on names
constantly (E. faecium is `Enterococcus_B faecium`; M. pneumoniae is
`Mycoplasmoides pneumoniae`), so joining our GenBank census to SPIRE by name
would silently drop exactly the species we are hunting. SPIRE's cluster table
already carries `size.pg3` -- the count of proGenomes3 ISOLATE genomes in the
same 95% ANI cluster -- so "existing isolates" and "new MAGs" can be compared
inside one clustering with no name matching at all. proGenomes3's inclusion
criteria are not identical to our GenBank census, so pg3 is a PROXY for how
much isolate material exists, not a restatement of our own pool.

THREE CLASSES, in the order they were asked for:

  1 crosses_50    pg3 < 50 and pg3 + qualified MAGs >= 50. Uses the existing
                  census threshold, so these are species that genuinely cross
                  the bar the project already set.
  2 mag_rich      pg3 < 50 and qualified MAGs >= 100. Enough material that
                  several k=12 panels might form -- worth an ANI check.
  3 mag_only      pg3 == 0 and qualified MAGs >= 50. Never had isolate
                  representation at all; this is where phylogenetic coverage
                  actually widens.

WHAT THIS CANNOT DECIDE, and must not be read as deciding:

  * PANELS. A SPIRE cluster is 95% ANI; our panels need mutual ~99% ANI AND
    alignment fraction. Members of one cluster are not all mutually
    comparable, so MAG count is an UPPER BOUND on panel material and a poor
    predictor of it. The panel column stays TBD until the genomes are
    downloaded and run through skani. Ranking here is provisional.
  * NET GAIN. A qualified MAG may be the same organism as an isolate we
    already hold, or two MAGs may come from re-assemblies of one sample.
    Counts here are CANDIDATES, before cross-database and sample dedup.

QUALITY BAR. completeness >= 90, contamination <= 5, gunc_pass. Deliberately
NOT filtered on n50: N50 is a global statistic and a locus only needs its own
flanks intact locally, so screening on it here would discard usable genomes for
the wrong reason. n50 tiers are reported alongside as information.
"""
import collections
import gzip
import sys

GEN, CLU = sys.argv[1], sys.argv[2]
OUT = sys.argv[3]
RUNNING = {s.strip().lower() for s in (sys.argv[4] if len(sys.argv) > 4 else "").split(":") if s.strip()}


def f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


# ---- cluster table: taxonomy + isolate (proGenomes3) representation --------
clu = {}
with gzip.open(CLU, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        if len(r) < len(hdr):
            continue
        clu[r[ix["spire_cluster"]]] = {
            "species": r[ix["species"]], "genus": r[ix["genus"]],
            "family": r[ix["family"]], "phylum": r[ix["phylum"]],
            "pg3": int(f(r[ix["size.pg3"]])),
            "spire": int(f(r[ix["size.spire"]]))}
print("clusters in table: %d" % len(clu))

# ---- genome table: qualified MAGs per cluster ------------------------------
q = collections.Counter()
n50_100 = collections.Counter()
n50_50 = collections.Counter()
samples = collections.defaultdict(set)
with gzip.open(GEN, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    has_smp = "derived_from_sample" in ix
    for line in fh:
        r = line.rstrip("\n").split("\t")
        if len(r) < len(hdr):
            continue
        if f(r[ix["completeness"]]) < 90 or f(r[ix["contamination"]], 100) > 5:
            continue
        g = r[ix["gunc_pass"]].strip().lower()
        if not (g in ("true", "t", "yes") or f(g, 0.0) == 1.0):
            continue
        c = r[ix["spire_cluster"]]
        q[c] += 1
        n = f(r[ix["n50"]])
        if n >= 50000:
            n50_50[c] += 1
        if n >= 100000:
            n50_100[c] += 1
        if has_smp:
            samples[c].add(r[ix["derived_from_sample"]])
print("clusters with >=1 qualified MAG: %d" % len(q))

rows = []
for c, n in q.items():
    m = clu.get(c)
    if not m:
        continue
    pg3 = m["pg3"]
    cls = []
    if pg3 < 50 and pg3 + n >= 50:
        cls.append("1_crosses_50")
    if pg3 < 50 and n >= 100:
        cls.append("2_mag_rich")
    if pg3 == 0 and n >= 50:
        cls.append("3_mag_only")
    if not cls:
        continue
    sp = m["species"] or c
    rows.append({
        "spire_cluster": c, "species": sp, "genus": m["genus"],
        "family": m["family"], "phylum": m["phylum"],
        "isolates_pg3": pg3, "qualified_mags": n,
        "mags_n50_50k": n50_50[c], "mags_n50_100k": n50_100[c],
        "distinct_samples": len(samples[c]) if samples else -1,
        "combined": pg3 + n,
        "classes": ";".join(cls),
        "already_running": "yes" if sp.strip().lower() in RUNNING else "no",
        "panels_formable": "TBD_needs_ANI"})

rows.sort(key=lambda r: -r["qualified_mags"])
cols = ["species", "spire_cluster", "phylum", "family", "genus",
        "isolates_pg3", "qualified_mags", "mags_n50_50k", "mags_n50_100k",
        "distinct_samples", "combined", "classes", "already_running",
        "panels_formable"]
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

print()
print("candidate clusters: %d" % len(rows))
for k in ("1_crosses_50", "2_mag_rich", "3_mag_only"):
    print("  %-14s %6d" % (k, sum(1 for r in rows if k in r["classes"])))
print("  already running %d" % sum(1 for r in rows if r["already_running"] == "yes"))
print()
print("TOP 30 by qualified MAGs (provisional -- panels not yet computed)")
print("%-44s %7s %8s %9s %8s  %s"
      % ("species", "pg3", "qualMAG", "n50>100k", "samples", "classes"))
for r in rows[:30]:
    print("%-44s %7d %8d %9d %8d  %s"
          % (r["species"][:44], r["isolates_pg3"], r["qualified_mags"],
             r["mags_n50_100k"], r["distinct_samples"], r["classes"]))
print()
print("wrote %s" % OUT)

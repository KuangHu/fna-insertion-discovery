#!/usr/bin/env python3
"""How many SPIRE MAGs could this pipeline actually use, per species?

SPIRE holds 1.16M medium+ quality MAGs. "Medium quality" is a CheckM2 statement
about completeness and contamination and says nothing about CONTIGUITY, which
is what this pipeline needs:

  Arm B   minigraph over whole assemblies; a fragmented pool contributes
          fragmentation, not alleles -- the same reason the GenBank census
          ranks on Complete+Chromosome rather than total
  Layer 2 a 500 bp anchor on EACH side of a locus, both inside ONE contig, so
          a locus needs >= ~1 kb + insert of contiguous sequence around it

So the filter that matters is n50 / max_contig_length, not completeness. This
counts MAGs per species at a ladder of thresholds, so the decision to download
is made on how many are USABLE rather than how many exist.

Tiers, cumulative and increasingly strict:

  all            every MAG of the species
  +qual          completeness >= 90, contamination <= 5
  +gunc          and gunc_pass (chimerism -- a binned chimera would create
                 false structural differences, the exact signal we call)
  +n50_50k       and n50 >= 50 kb
  +n50_100k      and n50 >= 100 kb
  +n50_250k      and n50 >= 250 kb

Also reported per species: how many distinct spire_clusters the MAGs fall into,
and how many distinct source samples, since 500 MAGs from one study is not the
same evidence as 500 from 500 studies.
"""
import collections
import gzip
import sys

META = sys.argv[1]
# ONE comma-separated token, underscores for spaces. slurm/analysis.sh
# word-splits ARGS on purpose (so multi-value flags survive), which turns
# 'Escherichia coli' into two tokens carrying literal quotes and matches
# nothing at all -- silently, since "no species matched" looks like a real
# answer about SPIRE rather than a quoting bug.
# COLON separated, not comma: `sbatch --export` uses commas as its OWN
# separator, so ARGS="a,b,c" reaches the job as ARGS="a" and b,c become stray
# export vars. The list silently truncates to its first element.
WANT = {w.strip().lower().replace("_", " ")
        for w in (sys.argv[2] if len(sys.argv) > 2 else "").split(":") if w.strip()}

TIERS = ["all", "+qual", "+gunc", "+n50_50k", "+n50_100k", "+n50_250k"]


def f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


counts = collections.defaultdict(lambda: collections.Counter())
clusters = collections.defaultdict(lambda: collections.defaultdict(set))
samples = collections.defaultdict(lambda: collections.defaultdict(set))
n50s = collections.defaultdict(list)

with gzip.open(META, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        if len(r) < len(hdr):
            continue
        sp = r[ix["species"]].strip()
        key = sp.lower()
        if WANT and not any(w in key for w in WANT):
            continue
        comp, cont = f(r[ix["completeness"]]), f(r[ix["contamination"]], 100)
        n50 = f(r[ix["n50"]])
        # gunc_pass is written as the FLOAT "1.0", while gunc_pass_5 in the
        # same row is the string "True". Testing for "1" made gunc_pass False
        # everywhere and drove every stricter tier to 0 -- which reads as
        # "no usable MAGs" rather than as a parsing bug.
        _gp = r[ix["gunc_pass"]].strip().lower()
        gp = _gp in ("true", "t", "yes") or f(_gp, 0.0) == 1.0
        cl = r[ix["spire_cluster"]]
        smp = r[ix["derived_from_sample"]] if "derived_from_sample" in ix else ""
        tiers = ["all"]
        if comp >= 90 and cont <= 5:
            tiers.append("+qual")
            if gp:
                tiers.append("+gunc")
                if n50 >= 50000:
                    tiers.append("+n50_50k")
                if n50 >= 100000:
                    tiers.append("+n50_100k")
                if n50 >= 250000:
                    tiers.append("+n50_250k")
        for t in tiers:
            counts[sp][t] += 1
            clusters[sp][t].add(cl)
            samples[sp][t].add(smp)
        n50s[sp].append(n50)

if not counts:
    print("no species matched %s" % sorted(WANT))
    sys.exit(0)

print("%-42s %9s %9s %9s %10s %11s %11s"
      % ("species", *TIERS))
for sp in sorted(counts, key=lambda s: -counts[s]["all"]):
    print("%-42s %9d %9d %9d %10d %11d %11d"
          % (sp, *[counts[sp][t] for t in TIERS]))

print()
print("at the +n50_100k tier -- distinct clusters and source samples")
print("%-42s %9s %9s %9s" % ("species", "MAGs", "clusters", "samples"))
for sp in sorted(counts, key=lambda s: -counts[s]["+n50_100k"]):
    if not counts[sp]["+n50_100k"]:
        continue
    print("%-42s %9d %9d %9d"
          % (sp, counts[sp]["+n50_100k"], len(clusters[sp]["+n50_100k"]),
             len(samples[sp]["+n50_100k"])))

print()
print("n50 distribution of ALL MAGs per species (kb)")
print("%-42s %8s %8s %8s %8s" % ("species", "p25", "median", "p75", "p90"))
for sp in sorted(counts, key=lambda s: -counts[s]["all"]):
    v = sorted(n50s[sp])
    if not v:
        continue
    q = lambda p: v[min(len(v) - 1, int(len(v) * p))] / 1000.0
    print("%-42s %8.1f %8.1f %8.1f %8.1f"
          % (sp, q(.25), q(.50), q(.75), q(.90)))

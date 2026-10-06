#!/usr/bin/env python3
"""Emit per-cluster MAG id lists for download, at the same quality bar the
census counted at.

The census counted MAGs passing completeness >= 90, contamination <= 5 and
gunc_pass. This must apply the IDENTICAL filter or the download will not match
the numbers the decision was made on -- a mismatch there is the kind of thing
that silently turns "704 candidates" into some other pool and makes every later
comparison unreadable.

NOT filtered on n50, deliberately: N50 is a global statistic and a locus only
needs its own flanks locally intact, so screening on it here would discard
usable genomes for the wrong reason. n50 is written into the manifest instead,
so a downstream step can stratify on it once panel formability is known.

One file per cluster, plus a combined manifest carrying the fields any later
step needs to stratify or deduplicate: n50, contigs, completeness, and the
source sample (two MAGs from one sample are not two observations).
"""
import gzip
import os
import sys

GEN = sys.argv[1]
OUTDIR = sys.argv[2]
CLUSTERS = [c.strip() for c in sys.argv[3].split(":") if c.strip()]


def f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


want = set(CLUSTERS)
rows = {c: [] for c in want}
with gzip.open(GEN, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        if len(r) < len(hdr):
            continue
        c = r[ix["spire_cluster"]]
        if c not in want:
            continue
        if f(r[ix["completeness"]]) < 90 or f(r[ix["contamination"]], 100) > 5:
            continue
        g = r[ix["gunc_pass"]].strip().lower()
        if not (g in ("true", "t", "yes") or f(g, 0.0) == 1.0):
            continue
        rows[c].append({
            "genome_id": r[ix["genome_id"]],
            "species": r[ix["species"]],
            "n50": int(f(r[ix["n50"]])),
            "n_contigs": int(f(r[ix["n_contigs"]])),
            "genome_size": int(f(r[ix["genome_size"]])),
            "completeness": r[ix["completeness"]],
            "contamination": r[ix["contamination"]],
            "sample": r[ix["derived_from_sample"]] if "derived_from_sample" in ix else ""})

os.makedirs(OUTDIR, exist_ok=True)
man = open(os.path.join(OUTDIR, "manifest.tsv"), "w")
man.write("genome_id\tspire_cluster\tspecies\tn50\tn_contigs\tgenome_size\t"
          "completeness\tcontamination\tsample\n")
tot = 0
for c in CLUSTERS:
    rs = rows[c]
    sp = rs[0]["species"].replace(" ", "_") if rs else c
    p = os.path.join(OUTDIR, "%s.ids" % sp)
    with open(p, "w") as fh:
        for r in rs:
            fh.write(r["genome_id"] + "\n")
    nsmp = len({r["sample"] for r in rs})
    print("%-34s %-16s %6d MAGs  %6d distinct samples  -> %s"
          % (sp, c, len(rs), nsmp, os.path.basename(p)))
    for r in rs:
        man.write("\t".join(str(r[k]) for k in
                            ("genome_id",)) + "\t" + c + "\t" +
                  "\t".join(str(r[k]) for k in
                            ("species", "n50", "n_contigs", "genome_size",
                             "completeness", "contamination", "sample")) + "\n")
    tot += len(rs)
man.close()
print()
print("total %d MAGs across %d clusters" % (tot, len(CLUSTERS)))
print("wrote %s/manifest.tsv" % OUTDIR)

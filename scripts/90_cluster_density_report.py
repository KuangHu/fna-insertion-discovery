#!/usr/bin/env python3
"""Database capacity measured as sampling topology, not genome count.

NOTE 2026-08-30: the union-find (single-linkage) path below is RETAINED ONLY FOR
REFERENCE and must not be used at this scale. Measured on the full 7,723-genome
all-vs-all, mean degree at ANI>=99 & AF>=85 is 382 (median 318, max 1231), and
still 95 at ANI>=99.7. Any single-linkage partition therefore collapses into one
giant component and the cluster-size table is degenerate.

Use `--mode neighbourhood` (the default). It reports, per genome, how many
partners clear the Layer 1 criteria. That needs no partition, does not chain,
and maps directly onto what a discovery run actually requires: a handful of
mutually comparable genomes at one locus, not a globally defined cluster.

The yield measurement says cluster density dominates: 6-genome clusters gave
5.0-6.8 P1 per genome, a 3-genome cluster gave 0.33. That is strongly
non-linear, so `N genomes` is the wrong number to describe a database with.
Three quantities predict pipeline yield instead:

    N>=3, N>=6, and the genomes actually sitting inside dense clusters

    F_productive        = genomes in ANI clusters of >=6 / all unique genomes
    N_productive        = number of ANI clusters with >=6 members

Singletons are NOT written off. A singleton is an unfilled cluster seed: if new
data lands around it, it becomes a new productive cluster. It contributes
nothing to Arm B *in the current pool*, which is a statement about sampling,
not about the genome.

PDS CROSS-CHECK, never a gate. Layer 1 clusters by ANI. NCBI Pathogen Detection
clusters the same isolates by SNP distance on a surveillance backbone. They
answer different questions, and 40.5% of this pool carries no PDS label at all,
so PDS is used only to ask whether an ANI cluster is many independent lineages
or one resequenced outbreak:

    D_lineage = n_PDS_clusters / n_genomes_with_a_PDS_label

An ANI cluster of 20 genomes spanning 8 PDS clusters is worth far more than one
of 20 genomes in a single PDS cluster, and this column is what tells them apart.

    90_cluster_density_report.py --ani ani_sparse.tsv --genomes ecoli_dedup.txt \\
        --pds cluster_list.tsv --min-ani 99 --min-af 85 --out density
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402

BANDS = [(1, 1, "1 (singleton)"), (2, 2, "2"), (3, 5, "3-5"),
         (6, 10, "6-10"), (11, 20, "11-20"), (21, 10 ** 9, ">20")]

# measured P1 yield per genome, from the 3 benchmark clusters
YIELD = {"1 (singleton)": 0.0, "2": None, "3-5": 0.33,
         "6-10": 5.9, "11-20": None, ">20": None}


def core(n):
    p = n.split("_")
    return p[1].split(".")[0] if len(p) > 1 else n


def sname(path):
    b = os.path.basename(path)
    for s in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(s):
            b = b[: -len(s)]
    return b



NB_BANDS = [(0, 0, "0 (isolated)"), (1, 2, "1-2"), (3, 5, "3-5"),
            (6, 10, "6-10"), (11, 50, "11-50"), (51, 10 ** 9, ">50")]


def neighbourhood(args):
    """Per-genome partner counts. No partition, so no chaining."""
    import glob as _glob
    paths = sorted(_glob.glob(args.ani_glob)) if args.ani_glob else [args.ani]
    log("%d ANI file(s)" % len(paths))
    allg = [sname(l.strip()) for l in open(args.genomes) if l.strip()]
    deg = collections.Counter()
    seen = set()
    n_pairs = 0
    for p in paths:
        with open(p) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                r_i, q_i = head.index("Ref_file"), head.index("Query_file")
                a_i = head.index("ANI")
                f1 = head.index("Align_fraction_ref")
                f2 = head.index("Align_fraction_query")
            except ValueError:
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(r_i, q_i, a_i, f1, f2):
                    continue
                q = sname(f[q_i])
                seen.add(q)
                n_pairs += 1
                if f[r_i] == f[q_i]:
                    continue
                try:
                    if float(f[a_i]) < args.min_ani:
                        continue
                    if min(float(f[f1]), float(f[f2])) < args.min_af:
                        continue
                except ValueError:
                    continue
                deg[q] += 1
    log("%d pairs read; %d genomes appeared as a query" % (n_pairs, len(seen)))

    pds_of = {}
    if args.pds and os.path.exists(args.pds):
        with open(args.pds) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = r.get("gencoll_acc", "")
                if g.startswith(("GCA_", "GCF_")):
                    pds_of[core(g)] = r["PDS_acc"]

    with open(args.out + "_neighbourhood.tsv", "w") as fh:
        fh.write("genome\tn_neighbours\tpds_acc\n")
        for g in sorted(seen):
            fh.write("%s\t%d\t%s\n" % (g, deg[g], pds_of.get(core(g), ".")))

    n = len(seen)
    log("=" * 84)
    log("NEIGHBOURHOOD DENSITY  (ANI>=%.1f, AF>=%.1f, %d genomes queried of %d)"
        % (args.min_ani, args.min_af, n, len(allg)))
    log("")
    log("  %-16s %9s %9s   %s" % ("neighbours", "genomes", "%", "meaning"))
    log("  " + "-" * 70)
    meaning = {"0 (isolated)": "cannot enter any comparison",
               "1-2": "pairwise only",
               "3-5": "marginal; 3-genome clusters yielded 0.33 P1/genome",
               "6-10": "productive; 6-genome clusters yielded 5.0-6.8",
               "11-50": "productive, room to subsample",
               ">50": "saturated; selection matters more than acquisition"}
    for lo, hi, name in NB_BANDS:
        c = sum(1 for g in seen if lo <= deg[g] <= hi)
        log("  %-16s %9d %8.1f%%   %s"
            % (name, c, 100.0 * c / n if n else 0, meaning[name]))
    ok = sum(1 for g in seen if deg[g] >= 5)
    log("")
    log("  genomes with >=5 partners : %d  (%.1f%%)  <- can sit in a 6-genome set"
        % (ok, 100.0 * ok / n if n else 0))
    if deg:
        v = sorted(deg.values())
        log("  degree: median %d, p90 %d, max %d"
            % (v[len(v) // 2], v[int(len(v) * 0.9)], v[-1]))
    log("")
    log("  Density is a property of the POOL, not of any sample drawn from it.")
    log("  The 20-genome benchmark gave clusters of 6/6/3 because 20 genomes")
    log("  drawn from 7,723 are mostly unrelated -- not because E. coli is")
    log("  sparsely sampled. Build discovery sets from dense neighbourhoods.")
    log("=" * 84)
    log("wrote %s_neighbourhood.tsv" % args.out)
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ani", default=None, help="single skani output file")
    ap.add_argument("--genomes", required=True, help="the deduplicated genome list")
    ap.add_argument("--pds", default=None, help="PDG cluster_list.tsv")
    ap.add_argument("--min-ani", type=float, default=99.0)
    ap.add_argument("--min-af", type=float, default=85.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["neighbourhood", "singlelinkage"],
                    default="neighbourhood")
    ap.add_argument("--ani-glob", default=None,
                    help="glob for sharded skani dist output; overrides --ani")
    args = ap.parse_args()
    if not args.ani and not args.ani_glob:
        ap.error("give --ani or --ani-glob")

    if args.mode == "neighbourhood":
        return neighbourhood(args)

    allg = [sname(l.strip()) for l in open(args.genomes) if l.strip()]
    idx = {g: i for i, g in enumerate(allg)}
    parent = list(range(len(allg)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    n_pairs = n_kept = 0
    with open(args.ani) as fh:
        head = fh.readline()
        cols = head.rstrip("\n").split("\t")
        r_i = cols.index("Ref_file") if "Ref_file" in cols else 0
        q_i = cols.index("Query_file") if "Query_file" in cols else 1
        a_i = cols.index("ANI") if "ANI" in cols else 2
        f1 = cols.index("Align_fraction_ref") if "Align_fraction_ref" in cols else 3
        f2 = cols.index("Align_fraction_query") if "Align_fraction_query" in cols else 4
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(r_i, q_i, a_i, f1, f2):
                continue
            n_pairs += 1
            try:
                ani = float(f[a_i])
                af = min(float(f[f1]), float(f[f2]))
            except ValueError:
                continue
            # Layer 1's criteria: identity is necessary, aligned fraction is the
            # real guard -- two genomes can be 99% identical over 30% of themselves
            if ani < args.min_ani or af < args.min_af:
                continue
            a, b = sname(f[r_i]), sname(f[q_i])
            if a in idx and b in idx:
                union(idx[a], idx[b])
                n_kept += 1
    log("%d pairs read, %d passed ANI>=%.1f AF>=%.1f"
        % (n_pairs, n_kept, args.min_ani, args.min_af))

    comp = collections.defaultdict(list)
    for g, i in idx.items():
        comp[find(i)].append(g)
    clusters = sorted(comp.values(), key=len, reverse=True)

    pds_of = {}
    if args.pds and os.path.exists(args.pds):
        with open(args.pds) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = r.get("gencoll_acc", "")
                if g.startswith(("GCA_", "GCF_")):
                    pds_of[core(g)] = r["PDS_acc"]
        log("%d genomes carry a PDS label" % len(pds_of))

    rows = []
    for k, mem in enumerate(clusters):
        labelled = [pds_of[core(g)] for g in mem if core(g) in pds_of]
        npds = len(set(labelled))
        big = (collections.Counter(labelled).most_common(1)[0][1] / len(labelled)
               if labelled else float("nan"))
        rows.append({"cluster_id": "ac%05d" % k, "n_genomes": len(mem),
                     "n_pds_clusters": npds,
                     "frac_pds_labelled": len(labelled) / len(mem),
                     "largest_pds_fraction": big,
                     "D_lineage": npds / len(labelled) if labelled else float("nan"),
                     "members": ",".join(sorted(mem)[:50])})

    with open(args.out + "_clusters.tsv", "w") as fh:
        cols = ["cluster_id", "n_genomes", "n_pds_clusters", "frac_pds_labelled",
                "largest_pds_fraction", "D_lineage", "members"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.3f" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    tot = len(allg)
    log("=" * 88)
    log("ANI CLUSTER DENSITY  (ANI>=%.1f, AF>=%.1f, %d unique genomes)"
        % (args.min_ani, args.min_af, tot))
    log("")
    log("  %-16s %9s %9s %10s %12s %14s"
        % ("cluster size", "clusters", "genomes", "% genomes", "med D_lineage",
           "P1/genome (obs)"))
    log("  " + "-" * 78)
    for lo, hi, name in BANDS:
        sel = [r for r in rows if lo <= r["n_genomes"] <= hi]
        g = sum(r["n_genomes"] for r in sel)
        dl = sorted(r["D_lineage"] for r in sel if r["D_lineage"] == r["D_lineage"])
        med = dl[len(dl) // 2] if dl else float("nan")
        y = YIELD.get(name)
        log("  %-16s %9d %9d %9.1f%% %12s %14s"
            % (name, len(sel), g, 100.0 * g / tot if tot else 0,
               "%.2f" % med if med == med else "n/a",
               "%.2f" % y if y is not None else "not measured"))
    prod = [r for r in rows if r["n_genomes"] >= 6]
    gp = sum(r["n_genomes"] for r in prod)
    p3 = [r for r in rows if r["n_genomes"] >= 3]
    log("")
    log("  N>=3 clusters      %6d   covering %6d genomes  (%.1f%%)"
        % (len(p3), sum(r["n_genomes"] for r in p3),
           100.0 * sum(r["n_genomes"] for r in p3) / tot))
    log("  N_productive (>=6) %6d   covering %6d genomes"
        % (len(prod), gp))
    log("  F_productive       %6.1f%%   <- fraction of the database that can yield P1"
        % (100.0 * gp / tot if tot else 0))
    log("")
    log("  Singletons are unfilled cluster seeds, not waste: adding relatives")
    log("  converts them into productive clusters. Download priority should")
    log("  maximise d(P1)/d(genome), which means filling size-2..5 clusters")
    log("  toward 6, not adding isolated accessions.")
    log("=" * 88)
    log("wrote %s_clusters.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

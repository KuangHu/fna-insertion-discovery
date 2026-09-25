#!/usr/bin/env python3
"""Build seed-centred comparison panels that are MUTUALLY comparable, not just seed-similar.

A high neighbour count does not mean a usable panel. A seed with 300 neighbours
at ANI>=99/AF>=85 tells you only that each of those 300 resembles the SEED; they
need not resemble each other. Single-linkage on this graph collapses into one
giant component precisely because that chaining is rampant, so "degree 300"
cannot be read as "300 genomes that can be compared together".

This module therefore grows a panel under an explicit mutual-comparability
constraint: a genome joins only if it clears the threshold against EVERY member
already in the panel, not merely against the seed. The panel is a clique in the
ANI graph (found greedily -- max clique is NP-hard and unnecessary here).

Among admissible candidates it maximises LINEAGE DIVERSITY, because redundancy
that is phylogenetically distributed is what supplies alternative alleles;
twenty resequenced isolates of one outbreak supply one observation. Priority:

  1. a PDS cluster not yet represented in the panel
  2. an unlabelled genome not near-identical (>= --dup-ani) to any member --
     a fallback independence proxy, since 38.7% of dense genomes carry no PDS
  3. otherwise the candidate with the highest minimum ANI to the panel

PDS is never a hard requirement. Gating on it would silently discard the
unlabelled 38.7% and bias every panel toward surveillance-sequenced lineages.

Reported per panel, so a coherent panel can be told from a merely large one:
  n_genomes, n_distinct_PDS, n_unlabelled, largest_PDS_fraction,
  median/min pairwise ANI, median/min AF, seed_degree,
  fraction_of_pairs_passing (1.000 by construction for the clique panels)

    91_seed_panel_sampler.py --edges edges_99_85.tsv --pds cluster_list.tsv \\
        --k 6 12 24 48 --n-seeds 20 --out panels
"""
import argparse
import collections
import csv
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def core(n):
    p = n.split("_")
    return p[1].split(".")[0] if len(p) > 1 else n


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edges", required=True,
                    help="compacted passing edges: a, b, ANI, minAF")
    ap.add_argument("--pds", default=None)
    ap.add_argument("--k", type=int, nargs="+", default=[6, 12, 24, 48])
    ap.add_argument("--n-seeds", type=int, default=20)
    ap.add_argument("--dup-ani", type=float, default=99.95,
                    help="unlabelled genomes at/above this to a panel member are "
                         "treated as the same lineage, not a new observation")
    ap.add_argument("--genome-dir", default=None,
                    help="if given, also write a manifest per panel")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--seed-pool", default=None,
                    help="file of genome names (one per line) to draw SEEDS "
                         "from. Used to make two panel sets INDEPENDENT rather "
                         "than nested: the S. aureus k=12 and k=24 sets draw "
                         "from disjoint halves, so neither set's events are a "
                         "subset of the other's. Panels may still recruit any "
                         "genome as a member -- only the seed is restricted.")
    ap.add_argument("--seed-pool-exclude", default=None,
                    help="file of genome names that must NOT be used as seeds")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    random.seed(args.seed)
    adj = collections.defaultdict(dict)
    with open(args.edges) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            a, b = f[0], f[1]
            try:
                ani, af = float(f[2]), float(f[3])
            except ValueError:
                continue
            adj[a][b] = (ani, af)
            adj[b][a] = (ani, af)
    log("%d genomes with >=1 passing edge" % len(adj))

    pds_of = {}
    if args.pds and os.path.exists(args.pds):
        with open(args.pds) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = r.get("gencoll_acc", "")
                if g.startswith(("GCA_", "GCF_")):
                    pds_of[core(g)] = r["PDS_acc"]
        log("%d PDS labels loaded" % len(pds_of))

    kmax = max(args.k)
    seeds = [g for g in adj if len(adj[g]) >= kmax]
    if args.seed_pool:
        keep = {l.strip() for l in open(args.seed_pool) if l.strip()}
        keep |= {os.path.basename(x)[:-4] for x in keep if x.endswith(".fna")}
        n0 = len(seeds)
        seeds = [g for g in seeds if g in keep]
        log("seed pool: %d of %d eligible seeds kept (%.1f%%)"
            % (len(seeds), n0, 100.0 * len(seeds) / n0 if n0 else 0))
        if n0 and not seeds:
            log("FATAL: the seed pool matched NO eligible seed -- a name-format "
                "mismatch between the pool file and the edge list")
            return 1
    if args.seed_pool_exclude:
        bad = {l.strip() for l in open(args.seed_pool_exclude) if l.strip()}
        bad |= {os.path.basename(x)[:-4] for x in bad if x.endswith(".fna")}
        n0 = len(seeds)
        seeds = [g for g in seeds if g not in bad]
        log("seed exclusion: %d of %d eligible seeds kept" % (len(seeds), n0))
    log("%d genomes have degree >= %d (eligible seeds)" % (len(seeds), kmax))
    if not seeds:
        log("no eligible seeds")
        return 1
    random.shuffle(seeds)
    seeds = seeds[: args.n_seeds]

    def grow(seed, k):
        """Greedy clique growth: a member must pass against ALL current members."""
        panel = [seed]
        seen_pds = {pds_of[core(seed)]} if core(seed) in pds_of else set()
        while len(panel) < k:
            best, best_key = None, None
            for c in adj[seed]:
                if c in panel:
                    continue
                # mutual comparability: c must be adjacent to every member
                if not all(c in adj[m] for m in panel):
                    continue
                mn = min(adj[m][c][0] for m in panel)
                p = pds_of.get(core(c))
                if p and p not in seen_pds:
                    tier = 0                      # new lineage
                elif not p and all(adj[m][c][0] < args.dup_ani for m in panel):
                    tier = 1                      # unlabelled, not a near-duplicate
                else:
                    tier = 2                      # redundant lineage
                key = (tier, -mn)
                if best_key is None or key < best_key:
                    best, best_key = c, key
            if best is None:
                break
            panel.append(best)
            p = pds_of.get(core(best))
            if p:
                seen_pds.add(p)
        return panel

    rows = []
    for seed in seeds:
        for k in sorted(args.k):
            panel = grow(seed, k)
            pairs = [(a, b) for i, a in enumerate(panel) for b in panel[i + 1:]]
            anis = [adj[a][b][0] for a, b in pairs if b in adj[a]]
            afs = [adj[a][b][1] for a, b in pairs if b in adj[a]]
            npass = len(anis)
            labels = [pds_of[core(g)] for g in panel if core(g) in pds_of]
            unl = len(panel) - len(labels)
            big = (collections.Counter(labels).most_common(1)[0][1] / len(labels)
                   if labels else float("nan"))
            rows.append({
                "seed": seed, "k_requested": k, "n_genomes": len(panel),
                "n_distinct_PDS": len(set(labels)), "n_unlabelled": unl,
                "largest_PDS_fraction": big,
                "median_pairwise_ANI": statistics.median(anis) if anis else float("nan"),
                "min_pairwise_ANI": min(anis) if anis else float("nan"),
                "median_AF": statistics.median(afs) if afs else float("nan"),
                "min_AF": min(afs) if afs else float("nan"),
                "seed_degree": len(adj[seed]),
                "fraction_of_pairs_passing": npass / len(pairs) if pairs else float("nan"),
                "members": ",".join(panel)})
            if args.genome_dir:
                d = os.path.join(args.out + "_manifests")
                os.makedirs(d, exist_ok=True)
                with open(os.path.join(d, "%s_k%02d.manifest" % (seed, k)), "w") as fh:
                    for g in panel:
                        fh.write(os.path.join(args.genome_dir, g + ".fna") + "\n")

    cols = ["seed", "k_requested", "n_genomes", "n_distinct_PDS", "n_unlabelled",
            "largest_PDS_fraction", "median_pairwise_ANI", "min_pairwise_ANI",
            "median_AF", "min_AF", "seed_degree", "fraction_of_pairs_passing",
            "members"]
    with open(args.out + "_panels.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.3f" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    log("=" * 90)
    log("SEED-CENTRED PANELS  (mutual comparability enforced, %d seeds)" % len(seeds))
    log("")
    log("  %-4s %9s %9s %10s %11s %11s %11s"
        % ("k", "achieved", "distinct", "unlabelled", "med pw ANI", "min pw ANI",
           "largest PDS"))
    log("  " + "-" * 74)
    for k in sorted(args.k):
        sel = [r for r in rows if r["k_requested"] == k]
        if not sel:
            continue
        ach = statistics.median(r["n_genomes"] for r in sel)
        dp = statistics.median(r["n_distinct_PDS"] for r in sel)
        un = statistics.median(r["n_unlabelled"] for r in sel)
        ma = statistics.median(r["median_pairwise_ANI"] for r in sel
                               if r["median_pairwise_ANI"] == r["median_pairwise_ANI"])
        mn = statistics.median(r["min_pairwise_ANI"] for r in sel
                               if r["min_pairwise_ANI"] == r["min_pairwise_ANI"])
        lp = [r["largest_PDS_fraction"] for r in sel
              if r["largest_PDS_fraction"] == r["largest_PDS_fraction"]]
        log("  %-4d %9.0f %9.0f %10.0f %11.3f %11.3f %11s"
            % (k, ach, dp, un, ma, mn,
               "%.2f" % statistics.median(lp) if lp else "n/a"))
    short = [r for r in rows if r["n_genomes"] < r["k_requested"]]
    log("")
    log("  panels that could NOT reach k: %d of %d" % (len(short), len(rows)))
    log("  (a shortfall means the seed's neighbours are not mutually comparable,")
    log("   which is exactly what raw degree hides)")
    log("=" * 90)
    log("wrote %s_panels.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

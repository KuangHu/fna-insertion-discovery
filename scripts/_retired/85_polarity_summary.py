#!/usr/bin/env python3
"""Layer 3 summary -- failure waterfall and P1 taxon-sampling stability.

Two questions, and the second is the one that decides whether Layer 3 works.

WATERFALL. P1/P2/P0 counts alone do not say where the layer loses loci, so the
loci are walked through the gates in order and the survivors counted at each
step. That shows whether the bottleneck is outgroup information, local
genealogy, ancestral ambiguity or support -- each of which has a different fix,
and only one of which (outgroup information) can be bought with more genomes.

P1 STABILITY. A P1 label is only worth having if it survives changing which
genomes were sampled. Re-running Layer 3c against outgroup sets built at
several ANI ceilings gives, per locus:

    P1 -> same allele        the label held
    P1 -> different allele   THE FAILURE THAT MATTERS: a confident call that
                             reverses under resampling. Majority voting scored
                             34.4% here; if the phylogenetic gates work this
                             should be near zero.
    P1 -> unresolved         acceptable: the gates withdrew the claim when the
                             evidence thinned, rather than keeping it.
    unresolved -> P1         rescue, reported for completeness.

A high P1 -> unresolved rate is not a defect. A non-trivial P1 -> different
allele rate is.

    85_polarity_summary.py --main polarity/all_cl0000 polarity/all_cl0001 \\
        --perturb 98.5:polarity_sens/98.5 98:polarity_sens/98 --out summary
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def load(dirs):
    out = {}
    for d in dirs:
        p = os.path.join(d, "polarity.tsv") if os.path.isdir(d) else d
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[r["locus_id"]] = r
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main", nargs="+", required=True)
    ap.add_argument("--perturb", nargs="*", default=[],
                    help="TAG:dir[,dir...] runs at other outgroup ANI ceilings")
    ap.add_argument("--refinements", nargs="*", default=[],
                    help="refinement.tsv files, for clade size vs success")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-lineages", type=int, default=2)
    args = ap.parse_args()

    main_rows = load(args.main)
    n = len(main_rows)
    log("=" * 78)
    log("LAYER 3 POLARITY SUMMARY  (%d Layer-2 resolved loci)" % n)

    # ------------------------------------------------------------ waterfall
    def num(r, k):
        try:
            return int(r.get(k, 0) or 0)
        except ValueError:
            return 0

    s1 = [r for r in main_rows.values()
          if num(r, "n_informative_outgroups") >= args.min_lineages]
    s2 = [r for r in s1 if num(r, "n_taxa") > 0]
    s3 = [r for r in s2 if r.get("monophyly_LR") == "yes"]
    s4 = [r for r in s3 if r.get("local_genealogy_conflict") != "yes"]
    s5 = [r for r in s4 if r.get("ancestral_allele") not in (".", "", None)]
    p1 = [r for r in s5 if r.get("polarity_tier") == "P1"]
    p2 = [r for r in main_rows.values() if r.get("polarity_tier") == "P2"]

    steps = [("Layer-2 resolved loci", n),
             (">=%d informative outgroups" % args.min_lineages, len(s1)),
             ("flanks placed in enough taxa", len(s2)),
             ("focal clade locally monophyletic", len(s3)),
             ("no L/R genealogy conflict", len(s4)),
             ("unique ancestral state", len(s5)),
             ("bootstrap support pass -> P1", len(p1))]
    log("")
    log("  FAILURE WATERFALL")
    prev = None
    for label, v in steps:
        drop = "" if prev is None else "  (-%d)" % (prev - v)
        log("    %-42s %5d  %5.1f%%%s"
            % (label, v, 100.0 * v / n if n else 0, drop))
        prev = v
    log("")
    log("    P2 (supported, sampling-limited)          %5d  %5.1f%%"
        % (len(p2), 100.0 * len(p2) / n if n else 0))
    tiers = Counter(r.get("polarity_tier") for r in main_rows.values())
    log("    P0 unresolved                             %5d  %5.1f%%"
        % (tiers["P0"], 100.0 * tiers["P0"] / n if n else 0))

    log("")
    log("  UNRESOLVED BY REASON")
    for k, v in Counter(r.get("note") for r in main_rows.values()
                        if r.get("polarity_tier") == "P0" and r.get("note")
                        ).most_common():
        log("    %-42s %5d" % (k, v))
    fc = Counter(r.get("monophyly_failure_class") for r in main_rows.values()
                 if r.get("monophyly_failure_class") not in ("", "none", None))
    if fc:
        log("")
        log("  MONOPHYLY FAILURE TAXONOMY")
        for k, v in fc.most_common():
            log("    %-42s %5d" % (k, v))

    # ------------------------------------------------------- P1 stability
    rows_out = []
    if args.perturb:
        log("")
        log("  P1 TAXON-SAMPLING STABILITY")
        agg = Counter()
        for spec in args.perturb:
            tag, _, ds = spec.partition(":")
            pr = load(ds.split(","))
            same = diff = unres = rescue = miss = 0
            for lid, r in main_rows.items():
                q = pr.get(lid)
                if r.get("polarity_tier") == "P1":
                    if q is None:
                        miss += 1
                    elif q.get("polarity_tier") in ("P1", "P2"):
                        if q.get("ancestral_allele") == r.get("ancestral_allele"):
                            same += 1
                        else:
                            diff += 1
                            rows_out.append([lid, tag, r.get("ancestral_allele"),
                                             q.get("ancestral_allele"),
                                             "P1_FLIP_TO_DIFFERENT_ALLELE"])
                    else:
                        unres += 1
                elif q is not None and q.get("polarity_tier") == "P1":
                    rescue += 1
            tot = same + diff + unres
            agg["same"] += same
            agg["diff"] += diff
            agg["unres"] += unres
            log("    ceiling %-6s P1 kept=%-4d flipped=%-3d withdrawn=%-4d "
                "rescued=%-4d  (not run=%d)"
                % (tag, same, diff, unres, rescue, miss))
            if tot:
                log("        %-14s same=%.1f%%  DIFFERENT=%.1f%%  unresolved=%.1f%%"
                    % ("", 100.0 * same / tot, 100.0 * diff / tot,
                       100.0 * unres / tot))
        t = sum(agg.values())
        if t:
            log("")
            log("    OVERALL  P1 -> same allele      %5.1f%%" % (100.0 * agg["same"] / t))
            log("    OVERALL  P1 -> DIFFERENT allele %5.1f%%   <-- the number that matters"
                % (100.0 * agg["diff"] / t))
            log("    OVERALL  P1 -> unresolved       %5.1f%%   (acceptable: claim withdrawn)"
                % (100.0 * agg["unres"] / t))

    # -------------------------------------------- clade size / gap vs P1
    if args.refinements:
        log("")
        log("  FOCAL CLADE / OUTGROUP GAP")
        for p in args.refinements:
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    log("    %-28s clade=%s/%s floor=%s gap=%s met=%s"
                        % (os.path.basename(os.path.dirname(p)),
                           r["focal_clade_size"], r["discovery_cluster_size"],
                           r["ingroup_ani_floor"], r["ani_gap"],
                           r["gap_requirement_met"]))

    with open(args.out + "_p1_flips.tsv", "w") as fh:
        fh.write("locus_id\tperturbation\tmain_allele\tperturbed_allele\tclass\n")
        for r in rows_out:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(args.out + "_waterfall.tsv", "w") as fh:
        fh.write("step\tn\tpct\n")
        for label, v in steps:
            fh.write("%s\t%d\t%.1f\n" % (label, v, 100.0 * v / n if n else 0))
    log("=" * 78)
    log("wrote %s_{waterfall,p1_flips}.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

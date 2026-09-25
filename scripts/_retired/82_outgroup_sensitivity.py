#!/usr/bin/env python3
"""Layer 3b-sensitivity -- does the inferred ancestral allele survive removing close relatives?

Every in-band outgroup for these clusters sits at 98.2-99.0 ANI, hard against
the --max-ani ceiling; the 95-98 band is nearly empty. Close relatives are the
genomes most likely to SHARE the event under test, so a majority among them
could report a derived allele as ancestral with high apparent confidence. This
grid tests exactly that, by rebuilding the outgroup set at progressively lower
ANI ceilings and asking whether the supported allele changes.

Everything except --max-ani is held fixed: aligned fraction, the independence
rule, place_locus(), and the allele-matching thresholds.

Three things are tracked per cutoff:
  callability   can anchors still place as outgroups get more distant
  agreement     unanimous / majority / split among the outgroups
  STABILITY     does the same locus keep the same supported allele

The third is the point. A locus supporting allele B at 99, 98.5 and 98 is
strong evidence. A locus that supports A at 99 and B at 98 is evidence that the
99% neighbours share the derived event -- and that simple majority voting over
close relatives is not ancestry.

Output is one row per locus with the supported allele, support count and
callable count at each cutoff, plus `polarity_stable_across_cutoffs`.

    82_outgroup_sensitivity.py --workdir W --clusters cl0000 cl0001 cl0002 \\
        --cutoffs 99.0 98.5 98.0 97.5 --out sens
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd):
    return subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode


def supported_allele(rows):
    """Plurality allele among CONFIDENT votes only.

    A vote is confident when the outgroup placed cleanly AND committed to one
    allele. `ambiguous_between_*` assignments are deliberately not counted:
    at 98% ANI ordinary SNPs make identity to two similar alleles nearly equal,
    and letting those through would manufacture agreement.
    """
    votes = [r["allele_assignment"] for r in rows
             if r["allele_match_status"] != "unresolved"
             and r["allele_assignment"] not in (".", "")
             and not r["allele_assignment"].startswith("ambiguous_between")]
    if not votes:
        return ".", 0, 0, "none"
    c = Counter(votes)
    top, n = c.most_common(1)[0]
    if n == len(votes):
        agree = "unanimous"
    elif n > len(votes) / 2:
        agree = "majority"
    else:
        agree = "split"
    return top, n, len(votes), agree


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--clusters", nargs="+", required=True)
    ap.add_argument("--cutoffs", nargs="+", type=float,
                    default=[99.0, 98.5, 98.0, 97.5])
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=80.0)
    ap.add_argument("-n", "--n-outgroups", type=int, default=5)
    ap.add_argument("--n-candidates", type=int, default=300)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--skip-run", action="store_true",
                    help="only aggregate existing per-cutoff output")
    args = ap.parse_args()

    W = args.workdir
    tags = ["%g" % c for c in args.cutoffs]

    for cut, tag in zip(args.cutoffs, tags):
        for c in args.clusters:
            od = os.path.join(W, "outgroup_sens", tag, c)
            mo = os.path.join(W, "og_sens", tag, c)
            if args.skip_run and os.path.exists(
                    os.path.join(mo, "outgroup_alleles.tsv")):
                continue
            os.makedirs(od, exist_ok=True)
            # the ingroup-vs-pool sweep is cutoff-independent: reuse it rather
            # than repeating a 7-minute skani per grid point
            src = os.path.join(W, "outgroup", c, "ingroup_vs_pool.tsv")
            dst = os.path.join(od, "ingroup_vs_pool.tsv")
            if os.path.exists(src) and not os.path.exists(dst):
                os.symlink(src, dst)
            log("[%s %s] selecting outgroups (max-ani %.1f)" % (tag, c, cut))
            run([sys.executable, os.path.join(_HERE, "80_outgroup_selector.py"),
                 "--ingroup", os.path.join(W, "armB_bench/stage10/clusters",
                                           c + ".manifest"),
                 "--pool-list", os.path.join(W, "outgroup/pool_all.txt"),
                 "--outdir", od, "--min-ani", str(args.min_ani),
                 "--max-ani", str(cut), "--min-af", str(args.min_af),
                 "-n", str(args.n_outgroups),
                 "--n-candidates", str(args.n_candidates),
                 "--threads", str(args.threads), "--reuse"])
            mf = os.path.join(od, "outgroups.manifest")
            if not (os.path.exists(mf) and os.path.getsize(mf)):
                log("  no outgroups at this cutoff")
                continue
            os.makedirs(mo, exist_ok=True)
            log("[%s %s] mapping alleles" % (tag, c))
            run([sys.executable, os.path.join(_HERE, "81_outgroup_allele_mapper.py"),
                 "--loci", os.path.join(W, "allele_recon/armB_" + c, "loci.tsv"),
                 "--alleles", os.path.join(W, "allele_recon/armB_" + c, "alleles.fna"),
                 "--seed-loci", os.path.join(W, "allele_recon/loci_armB_" + c + ".tsv"),
                 "--seed-genomes", os.path.join(W, "armB_bench/stage10/clusters",
                                                c + ".manifest"),
                 "--outgroups", mf, "--outdir", mo,
                 "--threads", str(args.threads)])

    # ---------------------------------------------------------- aggregate
    per = defaultdict(dict)
    ncall = defaultdict(dict)
    loci_all = set()
    for tag in tags:
        for c in args.clusters:
            p = os.path.join(W, "og_sens", tag, c, "outgroup_alleles.tsv")
            if not os.path.exists(p):
                continue
            by = defaultdict(list)
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["outgroup"] == ".":
                        loci_all.add(r["locus_id"])
                        continue
                    by[r["locus_id"]].append(r)
                    loci_all.add(r["locus_id"])
            for lid, rs in by.items():
                a, n, tot, agree = supported_allele(rs)
                per[lid][tag] = (a, n, tot, agree)
                ncall[lid][tag] = tot

    cols = []
    for tag in tags:
        cols += ["allele_at_%s" % tag, "support_%s" % tag, "n_callable_%s" % tag,
                 "agreement_%s" % tag]
    with open(args.out + "_per_locus.tsv", "w") as fh:
        fh.write("locus_id\t" + "\t".join(cols) +
                 "\tn_cutoffs_callable\tpolarity_stable_across_cutoffs\n")
        for lid in sorted(loci_all):
            row, alleles = [], []
            for tag in tags:
                v = per[lid].get(tag)
                if v and v[2] >= 2:
                    row += [v[0], v[1], v[2], v[3]]
                    alleles.append(v[0])
                else:
                    row += [".", 0, v[2] if v else 0, "none"]
            stable = ("yes" if len(set(alleles)) == 1 and len(alleles) >= 2
                      else "no" if len(set(alleles)) > 1
                      else "insufficient")
            fh.write("%s\t%s\t%d\t%s\n"
                     % (lid, "\t".join(map(str, row)), len(alleles), stable))

    log("=" * 78)
    log("OUTGROUP DISTANCE SENSITIVITY")
    log("")
    log("  %-10s %12s %14s %10s %10s %8s" %
        ("max-ANI", "callable>=2", "callability%", "unanimous", "majority", "split"))
    n_tot = len(loci_all)
    for tag in tags:
        ok = [lid for lid in loci_all
              if per[lid].get(tag) and per[lid][tag][2] >= 2]
        ag = Counter(per[lid][tag][3] for lid in ok)
        log("  %-10s %12d %13.1f%% %10d %10d %8d"
            % (tag, len(ok), 100.0 * len(ok) / n_tot if n_tot else 0,
               ag["unanimous"], ag["majority"], ag["split"]))
    log("")
    both = [lid for lid in loci_all
            if sum(1 for t in tags if per[lid].get(t) and per[lid][t][2] >= 2) >= 2]
    flips = [lid for lid in both
             if len({per[lid][t][0] for t in tags
                     if per[lid].get(t) and per[lid][t][2] >= 2}) > 1]
    log("  loci callable at >=2 cutoffs : %d" % len(both))
    log("  SAME allele at every cutoff  : %d  (%.1f%%)"
        % (len(both) - len(flips),
           100.0 * (len(both) - len(flips)) / len(both) if both else 0))
    log("  POLARITY FLIPS               : %d  (%.1f%%)"
        % (len(flips), 100.0 * len(flips) / len(both) if both else 0))
    log("")
    if flips:
        log("  A flip means the close (99%%) neighbours support one allele and the")
        log("  more distant ones another -- i.e. those neighbours share the event")
        log("  under test. Majority over close relatives is not ancestry.")
    log("=" * 78)
    log("wrote %s_per_locus.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

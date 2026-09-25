#!/usr/bin/env python3
"""Derive a species' ingroup ANI floor from ITS OWN distribution, by quantile.

§9 forbids hard-coding E. coli's thresholds elsewhere, and K. pneumoniae is the
case that shows why. Measured on the two species' full all-vs-all output:

    E. coli        median reported ANI 97.35     ANI 99.0 = quantile 0.927
    K. pneumoniae  median reported ANI 99.02     ANI 99.0 = quantile ~0.24

**ANI 99.0 sits on K. pneumoniae's MODE.** Applying E. coli's constant would
retain roughly three quarters of all reported pairs instead of the top ~7%, and
the neighbour count at 99.0/AF85 has a median of 2,383 -- it lumps distinct
clonal complexes into one component and leaves the threshold balanced on a peak,
where a 0.05 ANI shift moves enormous numbers of pairs.

THE TRANSFERABLE QUANTITY IS THE QUANTILE, NOT THE VALUE. E. coli's working
threshold sits at quantile 0.927 of its reported pairs; applying that same
quantile to K. pneumoniae gives **ANI 99.59**.

Both the chosen value AND its quantile are written to the params table, because
inter-species yield differences are uninterpretable without them -- a species
that yields fewer loci at a stricter effective threshold has not been shown to
carry fewer events.

CAVEAT recorded in the output: skani reports only pairs above roughly 80% ANI, so
these are quantiles of the REPORTED set, not of all pairs. The comparison is
consistent because both species were run the same way against their own full
pool, but the reported set's composition differs with pool size and diversity.

    98_species_ani_calibration.py --species kpneumoniae \\
        --shards kpneu/shards/ani_*.tsv --reference-quantile 0.9266 \\
        --out census/kpneumoniae
"""
import argparse
import bisect
import collections
import glob
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--shards", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reference-quantile", type=float, default=None,
                    help="quantile at which the REFERENCE species' working "
                         "threshold sat (E. coli 99.0 -> 0.9266). The threshold "
                         "for this species is read off its own distribution at "
                         "this quantile.")
    ap.add_argument("--reference-value", type=float, default=99.0)
    ap.add_argument("--min-af", type=float, default=85.0)
    ap.add_argument("--sample", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=3)
    args = ap.parse_args()

    random.seed(args.seed)
    files = []
    for pat in args.shards:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[")
                     else [pat])
    anis, deg_at = [], collections.Counter()
    npairs = 0
    for f in files:
        for i, line in enumerate(open(f)):
            if i == 0:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] == p[1]:
                continue
            try:
                a, af1, af2 = float(p[2]), float(p[3]), float(p[4])
            except ValueError:
                continue
            npairs += 1
            if args.sample < 1.0 and random.random() > args.sample:
                continue
            anis.append(a)
            if min(af1, af2) >= args.min_af:
                deg_at[p[1]] += 1
    if not anis:
        log("no pairs")
        return 1
    anis.sort()
    log("%s: %d reported pairs across %d shards (%d sampled)"
        % (args.species, npairs, len(files), len(anis)))

    def val_at(q):
        return anis[min(len(anis) - 1, int(len(anis) * q))]

    def quant_of(v):
        return bisect.bisect_left(anis, v) / len(anis)

    own_q = quant_of(args.reference_value)
    chosen_q = args.reference_quantile if args.reference_quantile else 0.9266
    chosen = val_at(chosen_q)

    rows = [("species", args.species),
            ("n_shards", len(files)),
            ("n_reported_pairs", npairs),
            ("median_reported_ANI", "%.3f" % val_at(0.50)),
            ("p75_ANI", "%.3f" % val_at(0.75)),
            ("p90_ANI", "%.3f" % val_at(0.90)),
            ("p95_ANI", "%.3f" % val_at(0.95)),
            ("p99_ANI", "%.3f" % val_at(0.99)),
            ("reference_value", "%.3f" % args.reference_value),
            ("reference_value_quantile_here", "%.4f" % own_q),
            ("reference_quantile", "%.4f" % chosen_q),
            ("CHOSEN_ingroup_ANI_floor", "%.3f" % chosen),
            ("pairs_retained_at_chosen", "%.2f%%" % (100 * (1 - chosen_q))),
            ("min_af", "%.1f" % args.min_af),
            ("quantiles_are_of", "skani-REPORTED pairs (>~80 ANI), not all pairs")]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_params.tsv", "w") as fh:
        fh.write("parameter\tvalue\n")
        for k, v in rows:
            fh.write("%s\t%s\n" % (k, v))

    log("=" * 88)
    log("ANI CALIBRATION -- %s" % args.species)
    log("")
    for k, v in rows:
        log("  %-32s %s" % (k, v))
    log("")
    log("  Applying the reference VALUE (%.2f) here would sit at quantile %.4f"
        % (args.reference_value, own_q))
    log("  and retain %.1f%% of reported pairs, versus %.1f%% at the transferred"
        % (100 * (1 - own_q), 100 * (1 - chosen_q)))
    log("  quantile. Record both, or inter-species yields cannot be compared.")
    if abs(own_q - chosen_q) > 0.2:
        log("")
        log("  *** THE CONSTANT DOES NOT TRANSFER. Do not reuse it. ***")
    log("=" * 88)
    log("wrote %s_params.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

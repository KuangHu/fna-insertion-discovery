#!/usr/bin/env python3
"""Frozen regression for the CLOSED Layer 3 outgroup architecture.

Layer 3's outgroup architecture was closed on 2026-08-29. The closure criteria
are deliberately NOT "coverage went up" -- an earlier criterion ("a mixed panel
must cut the 16 monophyly failures") was retired because it assumed the residual
came from outgroup choice, and the data showed it does not: intrusions fall
almost evenly on every panel member, near and far alike.

What is frozen instead is the property that actually matters, and which cost the
most to establish:

    changing taxon sampling may make a P1 call DISAPPEAR;
    it must never make a P1 call CHANGE ITS ANSWER.

Gates, with the 2026-08-29 cl0000 baseline in brackets:

  G1 topology-only taxa reach the tree      median n_taxa >= 9        [11]
     A taxon whose flanks place must enter the tree even when its allele is
     ambiguous or novel. Before the fix it was dropped entirely, which silently
     halved the far panel's tree count (39 vs 65) and made a denominator
     artefact look like a real monophyly improvement.

  G2 uncertain states are never forced      state<=topology always    [3 vs 5]
     state_informative may never exceed topology_informative.

  G3 no P1 allele flips under resampling    == 0                      [0]
     Across every panel tested (NEAR/FAR/MIXED), shared P1 loci agree on the
     ancestral allele: 19/19 and 25/25.

  G4 evidence loss withdraws, never flips   withdrawn allowed, flips 0
     P1 -> unresolved is correct behaviour and is not penalised.

  G5 persistent non-monophyly stays P0      label present             [16]
     Loci non-monophyletic against any panel are reported
     `local_focal_clade_nonmonophyly`, not rescued by re-drawing the clade.

    89_layer3_regression.py --primary polarity/all_cl0000 \\
        --perturb polarity_far2/cl0000 polarity_mixed/cl0000
"""
import argparse
import csv
import os
import statistics
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def load(d):
    p = os.path.join(d, "polarity.tsv") if os.path.isdir(d) else d
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        return {r["locus_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def num(r, k, d=0):
    try:
        return int(r.get(k) or d)
    except ValueError:
        return d


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--perturb", nargs="*", default=[])
    ap.add_argument("--min-median-taxa", type=int, default=9)
    args = ap.parse_args()

    P = load(args.primary)
    if not P:
        log("FATAL: no primary polarity table at %s" % args.primary)
        return 1
    ok = True

    taxa = [num(r, "n_taxa") for r in P.values()]
    med = statistics.median(taxa) if taxa else 0
    g1 = med >= args.min_median_taxa

    bad = [l for l, r in P.items()
           if r.get("n_state_informative_outgroups") not in (None, "")
           and num(r, "n_state_informative_outgroups")
           > num(r, "n_topology_informative_outgroups")]
    g2 = not bad

    flips, shared, withdrawn = 0, 0, 0
    detail = []
    for d in args.perturb:
        Q = load(d)
        for l, r in P.items():
            if r.get("polarity_tier") != "P1":
                continue
            q = Q.get(l)
            if q is None:
                continue
            if q.get("polarity_tier") in ("P1", "P2"):
                shared += 1
                if q.get("ancestral_allele") != r.get("ancestral_allele"):
                    flips += 1
                    detail.append((l, os.path.basename(d),
                                   r.get("ancestral_allele"),
                                   q.get("ancestral_allele")))
            else:
                withdrawn += 1
    g3 = flips == 0
    g4 = flips == 0

    notes = Counter(r.get("note") for r in P.values() if r.get("note"))
    n_nonmono = notes.get("local_focal_clade_nonmonophyly", 0) + \
        notes.get("ingroup_not_monophyletic_at_locus", 0)
    g5 = all(P[l].get("polarity_tier") == "P0" for l in P
             if P[l].get("note") in ("local_focal_clade_nonmonophyly",
                                     "ingroup_not_monophyletic_at_locus"))

    rows = [("G1 topology-only taxa reach the tree",
             "median n_taxa = %d" % med, ">= %d" % args.min_median_taxa, g1),
            ("G2 uncertain states never forced",
             "%d violations" % len(bad), "== 0", g2),
            ("G3 no P1 allele flips under resampling",
             "%d flips / %d shared" % (flips, shared), "== 0", g3),
            ("G4 evidence loss withdraws, never flips",
             "%d withdrawn, %d flipped" % (withdrawn, flips), "flips == 0", g4),
            ("G5 persistent non-monophyly stays P0",
             "%d labelled, all P0" % n_nonmono, "all P0", g5)]

    log("=" * 78)
    log("LAYER 3 REGRESSION  (outgroup architecture, frozen 2026-08-29)")
    log("")
    for name, val, thr, good in rows:
        ok &= good
        log("  %-42s %-24s %-12s %s"
            % (name, val, thr, "PASS" if good else "**FAIL**"))
    if detail:
        log("")
        log("  FLIPPED LOCI (this must be empty):")
        for l, d, a, b in detail[:10]:
            log("    %s  %s: %s -> %s" % (l, d, a, b))
    log("")
    log("  Primary panel is NEAR (98.2-99.0%% ANI). FAR and MIXED are")
    log("  perturbation validation only and must not be used to rescue")
    log("  monophyly -- the residual failures are a property of the locus.")
    log("  Layer 3 outgroup architecture: %s" % ("PASS" if ok else "FAIL"))
    log("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

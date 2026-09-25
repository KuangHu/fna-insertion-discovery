#!/usr/bin/env python3
"""Which features actually separate housekeeping repeats from true mobile elements?

A diagnostic, not a scorer. Before any S_repeat_penalty is defined, measure
whether the candidate features genuinely split the three known groups:

    TRUE_IS    ArmA+ / ISEScan+ families      (recovered known IS)
    NON_MGE    STRONG -> known_non_MGE_repeat (the false positives)
    RIBO       STRONG -> ribosomal_repeat     (rrn)

Only structural features are tested. ORF presence and transposase HMM hits are
DELIBERATELY EXCLUDED: making them part of the penalty would reintroduce the
annotation bias the pipeline exists to avoid, and would down-rank exactly the
noncanonical / RNA-guided elements that are the eventual target.

Features, all computed from the stage-20 conservation profile and family table:

  H_out_50/100/200   cross-copy k-mer Jaccard at 50/100/200 bp OUTSIDE the
                     boundary, averaged over the two edges. A true mobile
                     element inserts into unrelated sites, so homology should
                     collapse immediately outside. Sustained outside homology
                     is the signature of a segmental / housekeeping repeat.
  dH                 H_inside - H_outside, the continuous version of boundary
                     sharpness. Small dH -> the "boundary" is not a boundary.
  outside_reach_bp   how far outside homology persists at >= 0.5 Jaccard.
                     Direct measure of extended co-linearity.
  C_contrast         internal identity - mean(flank Jaccard). High C is the
                     mobile-element signature: copies alike inside, unalike
                     outside.

Profile sign convention (verified against stage 20 output): on the L edge
offset < 0 is outside; on the R edge offset > 0 is outside.

    64_repeat_penalty_diagnostic.py --families-glob 'armA_f0/*_families.tsv' \\
        --profiles-glob 'armA_f0/*_profiles.tsv' \\
        --per-is bench60_f0_per_is.tsv --annot stage50_novel_f0/element_annotation.tsv \\
        --out repdiag
"""
import argparse
import csv
import glob
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402

OUT_OFFSETS = (50, 100, 200)
INSIDE_OFFSETS = (50, 100, 150, 200)


def fnum(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def med(vals):
    vals = [v for v in vals if v == v]
    return statistics.median(vals) if vals else float("nan")


def load_profiles(paths):
    """family_id -> {'out': {bp: [jaccard...]}, 'in': {bp: [...]}}"""
    prof = defaultdict(lambda: {"out": defaultdict(list), "in": defaultdict(list)})
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                off = int(r["offset_bp"])
                j = fnum(r["median_jaccard"])
                if j != j:
                    continue
                # L: outside is negative. R: outside is positive.
                outside = off < 0 if r["edge"] == "L" else off > 0
                prof[r["family_id"]]["out" if outside else "in"][abs(off)].append(j)
    return prof


def features(fam_row, pr):
    o, i = pr["out"], pr["in"]
    h_out = {bp: med(o.get(bp, [])) for bp in OUT_OFFSETS}
    h_in = med([v for bp in INSIDE_OFFSETS for v in i.get(bp, [])])
    h_out_mean = med([v for bp in OUT_OFFSETS for v in o.get(bp, [])])
    reach = 0
    for bp in sorted(o):
        if med(o[bp]) >= 0.5:
            reach = max(reach, bp)
    ident = fnum(fam_row.get("internal_identity_median")) / 100.0
    fl = fnum(fam_row.get("left_flank_jaccard"))
    fr = fnum(fam_row.get("right_flank_jaccard"))
    contrast = ident - (fl + fr) / 2.0 if fl == fl and fr == fr else float("nan")
    return {"H_out_50": h_out[50], "H_out_100": h_out[100], "H_out_200": h_out[200],
            "H_in": h_in, "dH": h_in - h_out_mean if h_in == h_in else float("nan"),
            "outside_reach_bp": reach, "C_contrast": contrast,
            "S_mobility": fnum(fam_row.get("S_mobility")),
            "n_copies": fnum(fam_row.get("n_copies")),
            "boundary_sharpness": fnum(fam_row.get("boundary_sharpness"))}


FEATS = ["H_out_50", "H_out_100", "H_out_200", "dH", "outside_reach_bp",
         "C_contrast", "boundary_sharpness", "S_mobility", "n_copies"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families-glob", required=True)
    ap.add_argument("--profiles-glob", required=True)
    ap.add_argument("--per-is", required=True)
    ap.add_argument("--annot", required=True, help="stage 50 element_annotation.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fam = {}
    for p in sorted(glob.glob(args.families_glob)):
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                fam[r["family_id"]] = r
    prof = load_profiles(sorted(glob.glob(args.profiles_glob)))
    log("%d families, %d with profiles" % (len(fam), len(prof)))

    # group 1: Arm A families that matched an ISEScan call
    true_is = set()
    with open(args.per_is) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["found_by_armA"] == "1" and r["armA_family"] != ".":
                true_is.add(r["armA_family"])

    # groups 2 and 3 from the stage 50 grading of the novel pool
    non_mge, ribo = set(), set()
    with open(args.annot) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("discovery_verdict") != "STRONG_MOBILE_CANDIDATE":
                continue
            if r.get("class_call") == "known_non_MGE_repeat":
                non_mge.add(r["element_id"])
            elif r.get("class_call") == "ribosomal_repeat":
                ribo.add(r["element_id"])

    groups = [("TRUE_IS", true_is), ("NON_MGE", non_mge), ("RIBO", ribo)]
    rows = []
    vals = {g: defaultdict(list) for g, _ in groups}
    for gname, ids in groups:
        for fid in sorted(ids):
            if fid not in fam or fid not in prof:
                continue
            f = features(fam[fid], prof[fid])
            rows.append([gname, fid] + ["%.4f" % f[k] if f[k] == f[k] else "NA"
                                        for k in FEATS])
            for k in FEATS:
                if f[k] == f[k]:
                    vals[gname][k].append(f[k])

    with open(args.out + "_per_family.tsv", "w") as fh:
        fh.write("group\tfamily_id\t" + "\t".join(FEATS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    log("=" * 88)
    log("FEATURE SEPARATION  (median per group; n in parentheses)")
    log("")
    hdr = "  %-20s" % "feature"
    for g, _ in groups:
        hdr += "%14s" % ("%s(%d)" % (g, len(vals[g].get("S_mobility", []))))
    log(hdr + "%14s" % "NON_MGE-TRUE")
    log("  " + "-" * 84)
    for k in FEATS:
        line = "  %-20s" % k
        ms = {}
        for g, _ in groups:
            m = med(vals[g][k])
            ms[g] = m
            line += "%14s" % ("%.3f" % m if m == m else "NA")
        d = (ms["NON_MGE"] - ms["TRUE_IS"]
             if ms["NON_MGE"] == ms["NON_MGE"] and ms["TRUE_IS"] == ms["TRUE_IS"]
             else float("nan"))
        line += "%14s" % ("%+.3f" % d if d == d else "NA")
        log(line)
    log("=" * 88)
    log("wrote %s_per_family.tsv" % args.out)


if __name__ == "__main__":
    main()

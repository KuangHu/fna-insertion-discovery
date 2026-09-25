#!/usr/bin/env python3
"""Benchmark 4 -- is Layer 2 reconstruction accurate? Arm B is the independent gold.

Arm B already produces an alternative allele for every event, by a completely
different route (minigraph bubbles + per-genome genotyping). So it can grade the
allele reconstructor without any read data.

The test is only meaningful because the loci were seeded with a PAD: anchors sit
--pad bp away from Arm B's claimed junction, so Layer 2 has to relocate the
junction from sequence alone. Seeding the anchors exactly at the junction would
make the comparison circular.

FRAME. A junction is only meaningful in the frame of a genome carrying the
pre-event allele. Layer 2 reports its junctions in the frame of its own chosen
pre-event carrier; Arm B reports them in `empty_ref_genome`. An event is scored
only when those two genomes agree -- otherwise the coordinates are not
comparable and the row is counted as `frame_mismatch`, not as an error.

    dL = L2.junction_L_abs - armB.junction_L
    dR = L2.junction_R_abs - armB.junction_R

BOTH junctions must be right; one-sided agreement is reported for contrast only,
exactly as in the stage 61 read benchmark.

    71_layer2_vs_armB.py --layer2 allele_recon/armB_cl0000/loci.tsv \\
        --armb armB_bench/stage30/cl0000/armB_events.tsv --out l2bench
"""
import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def num(v, d=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        try:
            return float(v)
        except (TypeError, ValueError):
            return d


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layer2", required=True, nargs="+", help="loci.tsv file(s)")
    ap.add_argument("--armb", required=True, nargs="+", help="armB_events.tsv file(s)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tiers", default="high",
                    help="comma list of Arm B qc_tier values to score, or 'all'")
    args = ap.parse_args()

    tiers = None if args.tiers == "all" else set(args.tiers.split(","))

    gold = {}
    for p in args.armb:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gold[r["event_id"]] = r
    l2 = {}
    for p in args.layer2:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                l2[r["locus_id"]] = r
    log("Arm B events: %d ; Layer 2 loci resolved: %d" % (len(gold), len(l2)))

    rows, dl, dr, dins, cls = [], [], [], [], Counter()
    status = Counter()
    for eid, g in sorted(gold.items()):
        if tiers and g.get("qc_tier") not in tiers:
            continue
        status["armB_scored"] += 1
        r = l2.get(eid)
        if r is None:
            status["not_resolved_by_layer2"] += 1
            rows.append([eid, g.get("qc_tier", "."), "NOT_RESOLVED",
                         "", "", "", "", "", ""])
            continue
        cls[r["event_class"]] += 1
        if r["pre_event_genome"] != g["empty_ref_genome"]:
            status["frame_mismatch"] += 1
            rows.append([eid, g.get("qc_tier", "."), "FRAME_MISMATCH",
                         r["event_class"], r["pre_event_genome"],
                         g["empty_ref_genome"], "", "", ""])
            continue
        jl, jr = num(r["junction_L_abs"]), num(r["junction_R_abs"])
        gl, gr = num(g["junction_L"]), num(g["junction_R"])
        if None in (jl, jr, gl, gr) or jl < 0 or jr < 0:
            status["no_junction"] += 1
            continue
        a, b = jl - gl, jr - gr
        di = num(r["inserted_len"], 0) - num(g["insert_size"], 0)
        dl.append(a)
        dr.append(b)
        dins.append(di)
        status["scored"] += 1
        rows.append([eid, g.get("qc_tier", "."), "SCORED", r["event_class"],
                     r["pre_event_genome"], g["empty_ref_genome"], a, b, di])

    with open(args.out + "_per_event.tsv", "w") as fh:
        fh.write("event_id\tarmB_qc_tier\tstatus\tlayer2_class\t"
                 "layer2_pre_event_genome\tarmB_empty_genome\tdL\tdR\t"
                 "d_inserted_len\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    def pct(f):
        return 100.0 * sum(1 for i in range(len(dl)) if f(dl[i], dr[i])) / len(dl) \
            if dl else float("nan")

    lines = []
    n = len(dl)
    if n:
        for k in (0, 1, 2, 5, 10):
            lines.append(("both_within_%d" % k,
                          pct(lambda a, b, k=k: abs(a) <= k and abs(b) <= k)))
            lines.append(("oneside_within_%d" % k,
                          pct(lambda a, b, k=k: abs(a) <= k or abs(b) <= k)))
    with open(args.out + "_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in status.items():
            fh.write("%s\t%d\n" % (k, v))
        for k, v in lines:
            fh.write("%s_pct\t%.1f\n" % (k, v))
        if n:
            fh.write("median_abs_dL\t%.1f\n" % statistics.median(map(abs, dl)))
            fh.write("median_abs_dR\t%.1f\n" % statistics.median(map(abs, dr)))
            fh.write("median_d_inserted_len\t%.1f\n" % statistics.median(dins))
            fh.write("exact_inserted_len_pct\t%.1f\n"
                     % (100.0 * sum(1 for x in dins if x == 0) / len(dins)))

    log("=" * 74)
    log("LAYER 2 vs ARM B   (tier=%s)" % args.tiers)
    for k in ("armB_scored", "scored", "not_resolved_by_layer2",
              "frame_mismatch", "no_junction"):
        if status.get(k):
            log("  %-26s %5d" % (k, status[k]))
    log("")
    if n:
        log("  %-22s %8s %8s" % ("tolerance", "both", "one-side"))
        for k in (0, 1, 2, 5, 10):
            log("  +/- %-18d %7.1f%% %7.1f%%" % (
                k,
                pct(lambda a, b, k=k: abs(a) <= k and abs(b) <= k),
                pct(lambda a, b, k=k: abs(a) <= k or abs(b) <= k)))
        log("")
        log("  median |dL| = %.1f bp   median |dR| = %.1f bp"
            % (statistics.median(map(abs, dl)), statistics.median(map(abs, dr))))
        log("  inserted length: median delta %.1f bp, exact on %.1f%%"
            % (statistics.median(dins),
               100.0 * sum(1 for x in dins if x == 0) / len(dins)))
    log("")
    log("  Layer 2 event_class over the same events:")
    for k, v in cls.most_common():
        log("    %-28s %5d" % (k, v))
    log("=" * 74)
    log("wrote %s_{per_event,summary}.tsv" % args.out)


if __name__ == "__main__":
    main()

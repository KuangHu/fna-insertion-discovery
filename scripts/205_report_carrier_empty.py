#!/usr/bin/env python3
"""Stratified report over the carrier empty-site inference output.

Counts are reported at three DIFFERENT grains because they are different
numbers and conflating them overstates coverage:
  events            distinct db_id
  (event,carrier)   distinct (db_id, carrier_accession) pairs -- the input unit
  occurrences       distinct located intervals -- the output unit

A boundary shift is NOT automatically an error. junction_overlap_bp is the
slide span: the insert can be placed anywhere within it and still yield the
same two alleles. A shift whose magnitude falls inside that span is consistent
with the event as called; only a shift exceeding it is evidence of a real
boundary problem. Both are counted separately rather than merged.
"""
import collections
import glob
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else \
    "/global/scratch/users/kh36969/fna_ins_discovery/inferred_empty_sample"


def rows():
    cols = None
    for f in sorted(glob.glob(OUT + "/parts/occ_*.tsv")):
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if cols is None:
                import importlib.util
                sp = importlib.util.spec_from_file_location(
                    "m", os.path.dirname(os.path.abspath(__file__))
                    + "/203_infer_carrier_empty.py")
                m = importlib.util.module_from_spec(sp)
                sp.loader.exec_module(m)
                cols = m.COLS
            if len(x) != len(cols):
                continue
            yield dict(zip(cols, x))


def pct(a, b):
    return "%5.1f%%" % (100.0 * a / b) if b else "  n/a"


def tab(title, counter, total):
    print("\n%s" % title)
    for k, v in counter.most_common():
        print("   %-34s %7d  %s" % (k, v, pct(v, total)))


def main():
    R = list(rows())
    ev = set(); pair = set(); occ = set()
    loc = collections.Counter(); rec = collections.Counter()
    shift = collections.Counter(); shift_vs_slide = collections.Counter()
    ldc = collections.Counter(); spec = collections.Counter()
    bycop = collections.Counter(); byov = collections.Counter()
    trunc = collections.Counter(); ids = []
    unsup_by_copies = collections.Counter()

    for r in R:
        ev.add(r["db_id"])
        pair.add((r["db_id"], r["carrier_accession"]))
        if r.get("insert_start", ".") not in (".", ""):
            occ.add(r["occurrence_id"])
        loc[r["localization_status"]] += 1
        rec[r["reconstruction_status"]] += 1
        trunc[r["candidate_search_truncated"]] += 1
        ldc[r["len_delta_class"]] += 1
        if r["localization_status"] == "unsupported_flanks":
            spec["unsupported " + r["species"]] += 1
        try:
            us, ds = int(r["up_shift"]), int(r["dn_shift"])
        except ValueError:
            continue
        mx = max(abs(us), abs(ds))
        shift["0" if mx == 0 else ("1-3" if mx <= 3 else
              ("4-10" if mx <= 10 else ">10"))] += 1
        try:
            sl = int(r["junction_overlap_bp"])
        except ValueError:
            sl = 0
        if mx == 0:
            shift_vs_slide["no shift"] += 1
        elif mx <= sl:
            shift_vs_slide["shift within slide span (consistent)"] += 1
        else:
            shift_vs_slide["shift EXCEEDS slide span (real boundary issue)"] += 1
        try:
            nc = int(r["n_candidates_examined"])
            bycop["1 copy" if nc == 1 else ("2-5" if nc <= 5 else
                  ("6-20" if nc <= 20 else ">20"))] += 1
            if r["localization_status"] == "unsupported_flanks":
                unsup_by_copies["1 copy" if nc == 1 else ("2-5" if nc <= 5 else
                                ("6-20" if nc <= 20 else ">20"))] += 1
        except ValueError:
            pass
        byov["slide 0" if sl == 0 else ("slide 1-3" if sl <= 3 else
             ("slide 4-10" if sl <= 10 else "slide >10"))] += 1
        try:
            ids.append(float(r["up_identity"]))
            ids.append(float(r["dn_identity"]))
        except ValueError:
            pass

    n = len(R)
    print("=" * 70)
    print("GRAINS (these are different numbers; do not conflate)")
    print("   rows in table                 %7d" % n)
    print("   distinct events               %7d" % len(ev))
    print("   distinct (event,carrier)      %7d" % len(pair))
    print("   distinct occurrences located  %7d" % len(occ))
    tab("LOCALIZATION STATUS (is this the corresponding copy?)", loc, n)
    tab("RECONSTRUCTION STATUS (can a clean empty be built?)", rec, n)
    tab("CANDIDATE SEARCH TRUNCATED", trunc, n)
    tab("COPIES OF THE INSERT EXAMINED IN THE CARRIER", bycop, n)
    tab("JUNCTION BOUNDARY SHIFT (max of the two sides)", shift, n)
    tab("SHIFT JUDGED AGAINST THE SLIDE SPAN", shift_vs_slide, n)
    tab("SLIDE SPAN (junction_overlap_bp, NOT a verified TSD)", byov, n)
    tab("LENGTH-DELTA CLASS", ldc, n)
    if spec:
        tab("UNSUPPORTED FLANKS BY SPECIES", spec, n)
    if unsup_by_copies:
        tab("UNSUPPORTED FLANKS BY COPY NUMBER", unsup_by_copies,
            sum(unsup_by_copies.values()))
    if ids:
        ids.sort()
        k = len(ids)
        print("\nANCHOR IDENTITY over both sides (n=%d)" % k)
        for lab, q in (("min", 0), ("p05", int(k * .05)), ("p25", int(k * .25)),
                       ("median", k // 2), ("p75", int(k * .75)), ("max", k - 1)):
            print("   %-8s %.4f" % (lab, ids[q]))


if __name__ == "__main__":
    main()

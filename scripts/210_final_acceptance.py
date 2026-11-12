#!/usr/bin/env python3
"""Final acceptance for inferred_empty_v2: coverage, per-row gate, delta join.

PART FILES DO NOT PROVE A COMPLETE RUN. Completeness is checked as a set
relation against the v1 input -- every (db_id, carrier_accession) pair present
in v1 must appear exactly once in v2 -- not by counting files or shards.

EVERY OFFICIAL ROW IS CHECKED INDIVIDUALLY, not merely bounded in aggregate.
The aggregate "official <= localization-resolved" is necessary but not
sufficient: it would still pass if a truncated row were swapped in for a
resolved one. Each official row must satisfy ALL of:
    localization_class_revised == resolved
    reconstructed_sequence_status == resolved
    candidate_search_truncated == no
    core_recoverable == yes
    a usable breakpoint range
Any violation is listed and the script exits non-zero.

THE LENGTH-DELTA JOIN IS EVENT -> OCCURRENCE. delta is a property of the EVENT
(it comes from the allele pair in the catalogue); official is a property of an
OCCURRENCE. One unexplained event can therefore contribute many official
occurrences, so both grains are reported.

A DELTA COINCIDING WITH THE SLIDE SPAN IS NOT EXPLAINED. Sliding preserves the
length of the removed interval, so it cannot cause a length difference; those
rows are counted as unexplained here.
"""
import collections
import glob
import importlib.util
import os
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
HERE = os.path.dirname(os.path.abspath(__file__))
V1 = W + "/inferred_empty_v1"
V2 = W + "/inferred_empty_v2"


def _cols(script, attr):
    sp = importlib.util.spec_from_file_location("x", HERE + "/" + script)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return getattr(m, attr)


def read(path, cols):
    out = []
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        if first != cols:
            fh.seek(0)
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if len(x) == len(cols):
                out.append(dict(zip(cols, x)))
    return out


def main():
    c1 = _cols("203_infer_carrier_empty.py", "COLS")
    c2 = _cols("206_revise_carrier_empty.py", "NEW_COLS")
    R1 = read(V1 + "/occurrences.tsv", c1)
    R2 = read(V2 + "/occurrences.tsv", c2)

    print("=" * 74)
    print("A. INPUT COVERAGE  (set relation, not a file count)")
    k1 = collections.Counter((r["db_id"], r["carrier_accession"]) for r in R1)
    k2 = collections.Counter((r["db_id"], r["carrier_accession"]) for r in R2)
    miss = set(k1) - set(k2)
    extra = set(k2) - set(k1)
    dup = [k for k, v in k2.items() if v > 1]
    print("   v1 rows %d   v2 rows %d" % (len(R1), len(R2)))
    print("   v1 distinct (event,carrier) %d   v2 %d" % (len(k1), len(k2)))
    print("   missing from v2 %d   unexpected in v2 %d   duplicated in v2 %d"
          % (len(miss), len(extra), len(dup)))

    print("\nB. FINAL COUNTS (v2)")
    loc = collections.Counter(r["localization_class_revised"] for r in R2)
    seq = collections.Counter(r["reconstructed_sequence_status"] for r in R2)
    bp = collections.Counter(r["breakpoint_status"] for r in R2)
    off = [r for r in R2 if r["official_fasta"] == "yes"]
    print("   events                      %7d" % len(set(r["db_id"] for r in R2)))
    print("   (event,carrier) pairs       %7d" % len(k2))
    print("   localization resolved       %7d" % loc["resolved"])
    print("   sequence resolved           %7d" % seq["resolved"])
    print("   OFFICIAL FASTA              %7d" % len(off))
    for title, c in (("localization", loc), ("sequence", seq), ("breakpoint", bp)):
        print("   -- %s" % title)
        for k, v in c.most_common():
            print("      %-44s %7d" % (k, v))

    print("\nC. PER-ROW OFFICIAL GATE (each row, not just the total)")
    bad = collections.Counter()
    for r in off:
        if r["localization_class_revised"] != "resolved":
            bad["localization not resolved"] += 1
        if r["reconstructed_sequence_status"] != "resolved":
            bad["sequence not resolved"] += 1
        if r.get("candidate_search_truncated") != "no":
            bad["candidate search truncated"] += 1
        if r.get("core_recoverable") != "yes":
            bad["core not recoverable"] += 1
        if r.get("breakpoint_range_start", ".") in (".", "") or \
           r.get("breakpoint_status") in (".", "", "not_attempted"):
            bad["no usable breakpoint range"] += 1
    if bad:
        for k, v in bad.most_common():
            print("   VIOLATION %-40s %7d" % (k, v))
    else:
        print("   all %d official rows satisfy every condition" % len(off))
    print("   ceiling: official %d <= localization resolved %d  %s"
          % (len(off), loc["resolved"],
             "OK" if len(off) <= loc["resolved"] else "VIOLATED"))

    print("\nD. LENGTH-DELTA x OFFICIAL  (event property joined to occurrences)")
    dcls = {}
    p = W + "/audits/delta_nonzero_classified.tsv"
    if os.path.exists(p):
        with open(p) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                dcls[f[ix["db_id"]]] = (f[ix["class"]], int(f[ix["delta"]]))
    UNEXP = ("unexplained", "DESCRIPTIVE")
    grp = collections.Counter()
    ev_off = collections.defaultdict(set)
    for r in R2:
        cl, d = dcls.get(r["db_id"], ("delta_zero", 0))
        unexp = any(t in cl for t in UNEXP)
        band = ("delta_zero" if cl == "delta_zero" else
                ("unexplained_<=3bp" if unexp and abs(d) <= 3 else
                 ("unexplained_4-50bp" if unexp and abs(d) <= 50 else
                  ("unexplained_>50bp" if unexp else "candidate_cause"))))
        grp[(band, r["official_fasta"])] += 1
        if r["official_fasta"] == "yes":
            ev_off[band].add(r["db_id"])
    print("   %-24s %10s %10s %10s" % ("delta band", "official", "withheld", "events_off"))
    for band in ("delta_zero", "candidate_cause", "unexplained_<=3bp",
                 "unexplained_4-50bp", "unexplained_>50bp"):
        y, n = grp[(band, "yes")], grp[(band, "no")]
        if y or n:
            print("   %-24s %10d %10d %10d" % (band, y, n, len(ev_off[band])))
    big = grp[("unexplained_>50bp", "yes")]
    print("\n   official occurrences from events with an UNEXPLAINED >50bp"
          " length difference: %d" % big)
    print("   (these pass the junction and core checks at the occurrence level;"
          " the event-level length difference remains unexplained)")

    fail = bool(miss or extra or dup or bad or len(off) > loc["resolved"])
    print()
    print("ACCEPTANCE: %s" % ("FAILED" if fail else "PASSED"))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()

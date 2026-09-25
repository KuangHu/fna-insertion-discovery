#!/usr/bin/env python3
"""Frozen regression test for the LOCKED Layer 2 reconstruction.

Layer 2 was locked on 2026-08-28 against the 251 Arm B events of clusters
cl0000/cl0001/cl0002. Any later change -- outgroup polarity, target context,
more genomes, a new cluster engine -- must re-run this and clear every gate.
The point is not to reproduce the numbers exactly (allele sets shift as inputs
change) but to catch a silent regression in the parts that were validated.

Gates, with the frozen 2026-08-28 baseline in brackets:

  resolved_loci                >= 225           [234]
  catastrophic_placements      == 0             [0]
  exact_inserted_len_pct       >= 89.0          [91.7]
  tolerant_inserted_len_pct    >= 71.0          [74.0]
  all_inserted_len_pct         >= 79.0          [82.0]
  agreement_regressions        == 0             [0]  (baseline 102)
  rc_invariant_pct             >= 95.0          [96.2]
  placement_ambiguous_pct      <= 8.0           [4.3]

A "catastrophic placement" is a locus whose longest allele exceeds 20x the seed
interval AND 10 kb -- the signature of two anchors landing on different copies
of a repeat, which fabricated a 42,367 bp allele before place_locus() was
rewritten.

What this does NOT test, by design: nucleotide-resolution junction accuracy.
Arm B's own junction_span has median 35 bp, so it cannot certify +/-2 bp. That
claim requires read-level gold and is explicitly not made.

    73_layer2_regression.py --workdir /path/fna_ins_discovery \\
        --clusters cl0000 cl0001 cl0002
"""
import argparse
import csv
import os
import sys
import types
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert                                                # noqa: E402
from lib.fasta import revcomp                                    # noqa: E402
from lib.util import log                                         # noqa: E402

_AUD = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "72_layer2_disagreement_audit.py")

BASELINE = {
    "resolved_loci": (">=", 225, 234),
    "catastrophic_placements": ("==", 0, 0),
    "exact_inserted_len_pct": (">=", 89.0, 91.7),
    "tolerant_inserted_len_pct": (">=", 71.0, 74.0),
    "all_inserted_len_pct": (">=", 79.0, 82.0),
    "agreement_regressions": ("==", 0, 0),   # vs --baseline-agreements, re-set to 102
    "rc_invariant_pct": (">=", 95.0, 96.2),
    "placement_ambiguous_pct": ("<=", 8.0, 4.3),
    # --- GATE 9, added 2026-09-02 -------------------------------------------
    # The eight gates above returned BIT-IDENTICAL values after a fix that
    # changed junction_ambiguity_bp on 99 loci. That proves none of them covers
    # the field -- and after the catalogue/verification reframing, that field
    # carries the TSD signal and most of the structural polarity evidence.
    # Frozen post-fix; the real validation (76/76 agreement with the
    # independent k-sweep on exact-path loci) is not a gate and cannot catch a
    # regression on its own.
    #
    # RE-FROZEN 2026-09-24, and the reason matters more than the numbers.
    # These seven, plus agreement_regressions, had been FAILING since the day
    # they were set and nobody had re-run the harness to see it. Values as
    # frozen on 2026-09-02, and as measured today:
    #
    #     overlap_4_15_loci        85 -> 83     overlap_pos_exact     57 -> 56
    #     overlap_zero_loci        18 -> 16     overlap_pos_tolerant  58 -> 52
    #     overlap_1_3_loci         21 -> 19     target_bases_lost_0  139 -> 135
    #     overlap_gt15_loci         9 ->  6     agreement_regressions  0 ->  7
    #
    # Cause, established by sweeping --min-insert over the same three clusters:
    # at --min-insert 0 SEVEN of the values return EXACTLY (85/21/9/57/58, plus
    # resolved_loci 234 and agreement_regressions 0) and the remaining two land
    # one locus away. The baseline was therefore taken from a run in which the
    # 500 bp floor was not in force on both paths -- the floor was added the
    # same day, and the harness was never re-run after it. The eleven loci that
    # move are all sub-500 bp (1, 220, 245, 300, 300, 300, 328, 390, 392, 403,
    # 453), i.e. the floor is working and the baseline was stale.
    #
    # Re-frozen against --min-insert 500 on post-fix code. Verified not to be
    # today's change: the pre-edit and post-edit 70_ differ in 0 of 11,700
    # cells on these clusters, so the re-freeze captures behaviour that was
    # already there.
    "overlap_4_15_loci": ("==", 83, 83),      # usable TSD band -> the T line
    "overlap_zero_loci": ("==", 16, 16),      # clean point insertions
    "overlap_1_3_loci": ("==", 19, 19),       # below the k>=4 chance floor
    "overlap_gt15_loci": ("==", 6, 6),        # repeat array / segmental dup
    "overlap_pos_exact": ("==", 56, 56),      # exact path -- 0 before the fix
    "overlap_pos_tolerant": ("==", 52, 52),   # tolerant path
    "target_bases_lost_zero": ("==", 135, 135),
    # --- GATE 10, added 2026-09-24 ------------------------------------------
    # Code review of 0262ce3 found that 86_/99_/100_ each sliced the LONG
    # allele with `lcp_bp`, which is a SHORT-allele offset. The two coincide on
    # the exact path, so the bug was invisible there and wrong on 31.1% of
    # tolerant-path loci across 19 species -- 11.5% of all catalogue rows
    # carried an insert_seq whose hash disagreed with its own inserted_md5.
    #
    # Gate 9 did not catch it and could not: every field it reads was correct.
    # The corrupted value was produced downstream, by re-deriving a sequence
    # 70_ had already computed. This gate asserts the re-derivation agrees, so
    # it fails if either 70_ stops emitting insert_start_in_long_bp or a
    # consumer goes back to slicing with lcp_bp.
    # Deliberately NOT a count. How many loci decompose is exactly the quantity
    # that has drifted here (133 when gate 9 was frozen, 124 today), so a count
    # would couple this gate to an unrelated question. What must hold is that
    # re-deriving the insert never needs the md5-anchored repair and never
    # fails: on the pre-fix code these read 13 and 0, on the post-fix code 0
    # and 0, so the pair discriminates the bug it was written for.
    "insert_frame_needing_repair": ("==", 0, 0),
    "insert_frame_unresolved": ("==", 0, 0),
}


def _helpers():
    src = open(_AUD).read().split("def main(")[0]
    m = types.ModuleType("aud")
    m.__dict__["__file__"] = _AUD
    exec(compile(src, "aud", "exec"), m.__dict__)
    return m


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--clusters", nargs="+", default=["cl0000", "cl0001", "cl0002"])
    ap.add_argument("--recon-prefix", default="allele_recon/armB_")
    ap.add_argument("--seed-prefix", default="allele_recon/loci_armB_")
    ap.add_argument("--armb-dir", default="armB_bench/stage30")
    ap.add_argument("--baseline-agreements", type=int, default=102,
                    help="loci that agreed with Arm B in the frozen run. Was "
                         "109, taken before the 500 bp floor applied on both "
                         "paths; 7 of the 11 loci the floor removes had agreed")
    args = ap.parse_args()

    W = args.workdir
    m = _helpers()
    loci, alleles, gold, seed = {}, defaultdict(dict), {}, {}
    for c in args.clusters:
        d = os.path.join(W, args.recon_prefix + c)
        with open(os.path.join(d, "loci.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                loci[r["locus_id"]] = r
        for k, v in m.read_fasta(os.path.join(d, "alleles.fna")).items():
            p = k.split("|")
            if len(p) >= 2:
                alleles[p[0]][p[1]] = v
        with open(os.path.join(W, args.armb_dir, c, "armB_events.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gold[r["event_id"]] = r
        with open(os.path.join(W, args.seed_prefix + c + ".tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                seed[r["locus_id"]] = int(r["end"]) - int(r["start"])

    got = {"resolved_loci": len(loci)}

    # GATE 10: the insert a consumer re-derives must be the insert 70_ hashed.
    n_fr = {"coord": 0, "lcp": 0, "repaired": 0, "unresolved": 0}
    for lid, r in loci.items():
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        ls = alleles.get(lid, {}).get(r["longest_allele"], "")
        if not ls:
            n_fr["unresolved"] += 1
            continue
        _, st = lib_insert.extract_insert(r, ls)
        n_fr[st] = n_fr.get(st, 0) + 1
    log("  insert frame: %s" % dict(n_fr))
    got["insert_frame_needing_repair"] = n_fr["repaired"]
    got["insert_frame_unresolved"] = n_fr["unresolved"]

    cat = 0
    for lid, r in loci.items():
        sd, L = seed.get(lid), int(r["longest_len"])
        if sd and L > 20 * sd and L > 10000:
            cat += 1
    got["catastrophic_placements"] = cat

    amb = sum(1 for r in loci.values() if r["placement_status"] == "ambiguous")
    got["placement_ambiguous_pct"] = 100.0 * amb / len(loci) if loci else 0.0

    by, agree = defaultdict(list), 0
    for lid, r in loci.items():
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        g = gold.get(lid)
        if not g:
            continue
        d = int(r["inserted_len"]) - int(g["insert_size"])
        by[r["decomposition_method"]].append(d)
        if d == 0:
            agree += 1

    def pct(v):
        return 100.0 * sum(1 for x in v if x == 0) / len(v) if v else 0.0
    got["exact_inserted_len_pct"] = pct(by.get("exact", []))
    got["tolerant_inserted_len_pct"] = pct(by.get("tolerant_alignment", []))
    allv = [x for v in by.values() for x in v]
    got["all_inserted_len_pct"] = pct(allv)
    got["agreement_regressions"] = max(0, args.baseline_agreements - agree)

    same = n = 0
    for lid, r in loci.items():
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        av = alleles.get(lid, {})
        s_, l_ = av.get(r["shortest_allele"]), av.get(r["longest_allele"])
        if not s_ or not l_ or len(l_) <= len(s_):
            continue
        f = m.decomp(s_, l_)
        rv = m.decomp(revcomp(s_), revcomp(l_))
        if not f or not rv:
            continue
        n += 1
        same += (f[0] == rv[0])
    got["rc_invariant_pct"] = 100.0 * same / n if n else 0.0

    # --- GATE 9: the junction-overlap distribution --------------------------
    band = Counter()
    pos_by_method = Counter()
    for lid, r in loci.items():
        if (r["event_class"] not in ("insertion_target_retained", "replacement")
                or int(r.get("inserted_len") or 0) <= 0):
            continue
        a = int(r.get("junction_ambiguity_bp") or 0)
        band["0" if a == 0 else "1_3" if a <= 3 else "4_15" if a <= 15
             else "gt15"] += 1
        if a > 0:
            pos_by_method[r.get("decomposition_method", "?")] += 1
    got["overlap_4_15_loci"] = band["4_15"]
    got["overlap_zero_loci"] = band["0"]
    got["overlap_1_3_loci"] = band["1_3"]
    got["overlap_gt15_loci"] = band["gt15"]
    got["overlap_pos_exact"] = pos_by_method["exact"]
    got["overlap_pos_tolerant"] = pos_by_method["tolerant_alignment"]
    got["target_bases_lost_zero"] = sum(
        1 for r in loci.values() if int(r.get("target_bases_lost") or 0) == 0)

    ok = True
    log("=" * 76)
    log("LAYER 2 REGRESSION  (frozen 2026-08-28 baseline)")
    log("")
    log("  %-30s %10s %10s %8s %s" % ("gate", "value", "threshold", "frozen", ""))
    log("  " + "-" * 72)
    for k, (op, thr, frozen) in BASELINE.items():
        v = got[k]
        good = (v >= thr if op == ">=" else v <= thr if op == "<=" else v == thr)
        ok &= good
        log("  %-30s %10.1f %4s %5.1f %8.1f %s"
            % (k, v, op, thr, frozen, "PASS" if good else "**FAIL**"))
    log("")
    log("  Layer 2 assembly-tier reconstruction: %s" % ("PASS" if ok else "FAIL"))
    log("  Nucleotide-resolution junction: NOT TESTED, NOT CLAIMED (needs read gold)")
    log("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

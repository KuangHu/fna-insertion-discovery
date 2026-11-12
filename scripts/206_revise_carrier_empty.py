#!/usr/bin/env python3
"""Revision pass over inferred_empty_v1: carrier-local breakpoint equivalence.

THIS NEVER WRITES INTO v1. It reads v1's occurrences.tsv and emits a separate
versioned directory, plus a change log naming every row whose classification
moved and why.

WHAT v1 GOT WRONG. v1 marked EVERY nonzero anchor shift boundary_ambiguous.
The slide-equivalence distinction existed only in the report script, never in
the output logic -- confirmed arithmetically on the sample: boundary_ambiguous
874 == within-slide 567 + beyond 307, and no-shift 4986 == resolved 4645 +
built_but_locus_unresolved 341.

WHY THE CATALOGUE'S SLIDE SPAN CANNOT BE REUSED HERE. junction_overlap_bp comes
from _slide() computed on the allele pair of ANOTHER assembly. A SNP or indel
in this carrier changes which placements are actually equivalent HERE. The
equivalence range is therefore recomputed on the carrier sequence itself, and
the catalogue value is carried only as a reference column.

THE RANGE IS SIGNED. Moving the interval left by one base preserves the product
iff carrier[s-1] == carrier[e-1]; moving right by one iff carrier[s] ==
carrier[e]. These extend independently, so the legal set is [-left, +right],
not a symmetric band. Comparing abs(shift) against a total width would accept a
shift in the illegal direction.

EQUIVALENCE IS VERIFIED ON SEQUENCE, NOT ON NUMBERS. For a fixed outer window
[s-F, e+F), the product at shift k is
    carrier[s-F : s+k] + carrier[e+k : e+F]
and two placements are equivalent only when those strings are byte-identical.
The arithmetic range is used to bound the search; the string comparison decides.

A COHERENT BOUNDARY MOVE SHIFTS BOTH SIDES TOGETHER. v1 scored the two anchors
independently, so up_shift and dn_shift can disagree. A genuine slide moves the
whole interval, giving up_shift == dn_shift. Disagreement is not a slide and is
never called equivalent.

TWO SEPARATE STATUSES. reconstructed_sequence_status answers "is the empty
sequence determined?"; breakpoint_status answers "is the cut point unique, a
range, or unresolved?". A canonical coordinate chosen inside an equivalent
range is a representative, NOT the unique biological breakpoint, so the range
is always reported alongside it.
"""
import argparse
import collections
import importlib.util
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("m203", HERE + "/203_infer_carrier_empty.py")
M = importlib.util.module_from_spec(_s)
_s.loader.exec_module(M)

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"

NEW_COLS = M.COLS + [
    "carrier_equiv_left", "carrier_equiv_right", "carrier_equiv_range",
    "shift_coherent", "chosen_shift", "shift_in_equiv_range",
    "sequence_equivalence_verified",
    "reconstructed_sequence_status", "breakpoint_status",
    "breakpoint_range_start", "breakpoint_range_end",
    "localization_class_revised", "official_fasta", "revision_note",
]


def equiv_range(cs, s, e, cap=200):
    """Signed [-left, +right] over which the interval may slide on THIS carrier."""
    left = 0
    while left < cap and s - 1 - left >= 0 and e - 1 - left >= 0 \
            and cs[s - 1 - left] == cs[e - 1 - left]:
        left += 1
    right = 0
    while right < cap and e + right < len(cs) and s + right < len(cs) \
            and cs[s + right] == cs[e + right]:
        right += 1
    return -left, right


def product(cs, s, e, k, F):
    lo = max(0, s - F)
    hi = min(len(cs), e + F)
    if s + k < lo or e + k > hi:
        return None
    return cs[lo:s + k] + cs[e + k:hi]


def relabel(up_id, dn_id, cand, trust):
    hi, lo = max(up_id, dn_id), min(up_id, dn_id)
    if hi >= trust and lo < cand:
        return "one_sided_anchor_support"
    if hi < cand:
        return "low_support_both_sides"
    return "insufficient_bilateral_support"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=W + "/inferred_empty_v1")
    ap.add_argument("--out", default=W + "/inferred_empty_v2")
    ap.add_argument("--manifest", default=SRC + "/genome_map.tsv")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--trust-id", type=float, default=0.90)
    ap.add_argument("--cand-min", type=float, default=0.80)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out + "/parts", exist_ok=True)

    man = M.load_map(a.manifest)
    rows = []
    with open(a.src + "/occurrences.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if len(x) == len(hdr):
                rows.append(dict(zip(hdr, x)))
    by_car = collections.defaultdict(list)
    for r in rows:
        by_car[r["carrier_accession"]].append(r)
    gs = sorted(by_car)
    mine = [g for i, g in enumerate(gs) if i % a.nshard == a.shard]
    sys.stderr.write("shard %d: %d carriers, %d rows\n"
                     % (a.shard, len(mine), sum(len(by_car[g]) for g in mine)))

    ftsv = open("%s/parts/occ_%03d.tsv" % (a.out, a.shard), "w")
    ffna = open("%s/parts/inferred_%03d.fna" % (a.out, a.shard), "w")
    fchg = open("%s/parts/changes_%03d.tsv" % (a.out, a.shard), "w")
    stat = collections.Counter()
    F = a.flank

    for n, g in enumerate(mine, 1):
        ctgs = None
        p = man.get(g)
        if p:
            ctgs = M.load_genome(p)
        for r in by_car[g]:
            out = dict(r)
            note = []
            loc = r["localization_status"]
            # descriptive relabel of the non-resolved localization classes
            if loc in ("unsupported_flanks", "one_sided_homology",
                       "no_homologous_locus", "low_confidence"):
                try:
                    ui, di = float(r["up_identity"]), float(r["dn_identity"])
                    new = (relabel(ui, di, a.cand_min, a.trust_id)
                           if loc != "low_confidence" else "low_confidence")
                except ValueError:
                    new = loc
                out["localization_class_revised"] = new
                if new != loc:
                    note.append("localization %s -> %s" % (loc, new))
            else:
                out["localization_class_revised"] = loc

            if r.get("insert_start", ".") in (".", "") or ctgs is None:
                out.update(reconstructed_sequence_status="not_attempted",
                           breakpoint_status="not_attempted",
                           official_fasta="no",
                           revision_note="; ".join(note) or ".")
                stat["not_attempted"] += 1
                ftsv.write("\t".join(str(out.get(c, ".")) for c in NEW_COLS) + "\n")
                continue

            ctg = r["carrier_contig"]
            cs = ctgs.get(ctg)
            s, e = int(r["insert_start"]), int(r["insert_end"])
            if cs is None:
                out.update(reconstructed_sequence_status="carrier_contig_missing",
                           breakpoint_status="not_attempted", official_fasta="no",
                           revision_note="; ".join(note) or ".")
                stat["contig_missing"] += 1
                ftsv.write("\t".join(str(out.get(c, ".")) for c in NEW_COLS) + "\n")
                continue

            lo_r, hi_r = equiv_range(cs, s, e)
            us, ds = int(r["up_shift"]), int(r["dn_shift"])
            coherent = (us == ds)
            k = us if coherent else None
            out["carrier_equiv_left"] = lo_r
            out["carrier_equiv_right"] = hi_r
            out["carrier_equiv_range"] = "[%d,%d]" % (lo_r, hi_r)
            out["shift_coherent"] = "yes" if coherent else "no"
            out["chosen_shift"] = k if coherent else "."

            in_range = coherent and (lo_r <= k <= hi_r)
            out["shift_in_equiv_range"] = "yes" if in_range else "no"

            verified = "."
            if coherent and k == 0:
                verified = "trivial"
            elif in_range:
                p0, pk = product(cs, s, e, 0, F), product(cs, s, e, k, F)
                verified = "yes" if (p0 is not None and p0 == pk) else "no"
            out["sequence_equivalence_verified"] = verified

            if verified in ("trivial", "yes"):
                seq_status = "resolved"
                bp = "unique" if (lo_r == 0 and hi_r == 0) else "equivalent_range"
            elif not coherent:
                seq_status = "unresolved"
                bp = "discordant_sides"
                note.append("up_shift %+d != dn_shift %+d" % (us, ds))
            else:
                seq_status = "unresolved"
                bp = "beyond_equivalence_range"
                note.append("shift %+d outside carrier range [%d,%d]" % (k, lo_r, hi_r))
            out["reconstructed_sequence_status"] = seq_status
            out["breakpoint_status"] = bp
            out["breakpoint_range_start"] = s + lo_r
            out["breakpoint_range_end"] = s + hi_r

            core_ok = r.get("core_recoverable") == "yes"
            trunc_ok = r.get("candidate_search_truncated") == "no"
            official = (out["localization_class_revised"] == "resolved"
                        and seq_status == "resolved" and core_ok and trunc_ok)
            out["official_fasta"] = "yes" if official else "no"
            stat["official" if official else "withheld"] += 1
            stat["bp_" + bp] += 1
            stat["seq_" + seq_status] += 1

            was_official = (r["localization_status"] == "resolved"
                            and r["reconstruction_status"] == "resolved")
            if official and not was_official:
                note.append("ADDED to official (breakpoint %s)" % bp)
                stat["added_vs_v1"] += 1
            if was_official and not official:
                note.append("REMOVED from official")
                stat["removed_vs_v1"] += 1
            out["revision_note"] = "; ".join(note) or "."
            if note:
                fchg.write("%s\t%s\t%s\n" % (r["occurrence_id"],
                                             r["localization_status"],
                                             out["revision_note"]))
            if official:
                seq = product(cs, s, e, k if coherent else 0, F)
                if seq:
                    ffna.write(">%s loc=%s bp=%s range=[%d,%d] shift=%s\n"
                               % (r["occurrence_id"],
                                  out["localization_class_revised"], bp,
                                  s + lo_r, s + hi_r, out["chosen_shift"]))
                    for i in range(0, len(seq), 80):
                        ffna.write(seq[i:i + 80] + "\n")
            ftsv.write("\t".join(str(out.get(c, ".")) for c in NEW_COLS) + "\n")
        if n % 200 == 0:
            sys.stderr.write("  %d/%d carriers\n" % (n, len(mine)))
    for fh in (ftsv, ffna, fchg):
        fh.close()
    print("shard %d: %s" % (a.shard, dict(stat)))


if __name__ == "__main__":
    main()

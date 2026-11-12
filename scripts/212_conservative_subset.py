#!/usr/bin/env python3
"""Flag the full v2 set, then carve a conservative analysis subset from it.

TWO DIFFERENT CLAIMS, KEPT APART.

  PASSES THE INFERENCE RULES  -- what official_fasta already means: the
  occurrence was localized, its sequence was determined, its candidate search
  was not truncated, and it carries a valid breakpoint interval.

  FIT FOR PRIMARY ANALYSIS    -- additionally, no UNRESOLVED anomaly attaches
  to the evidence THIS occurrence used. An anomaly is only disqualifying when
  it touches the reference locus this row was built from, or the event-level
  length arithmetic this row inherits.

Nothing is deleted. The full v2 set keeps every row and gains the flag
columns; the conservative subset is a filtered view with a stated rule.

EXCLUSIONS ARE A UNION OVER occurrence_id, NOT A SUM. The Audit A set and the
length-difference set overlap, so adding their sizes double-counts.
"""
import argparse
import collections
import glob
import os

W = "/global/scratch/users/kh36969/fna_ins_discovery"
V2 = W + "/inferred_empty_v2"

# Audit A verdicts that contradict emptiness at the locus this row used.
A_CONTRADICT = {"element_present_at_locus",
                "element_present_right_end_inexact"}
# Audit A verdicts that are unresolved rather than refuted: a left-junction
# match whose span does not reconstruct the element. Repeat context makes
# these ambiguous; they are neither cleared nor condemned.
A_UNRESOLVED = {"partial_element_at_locus",
                "left_junction_only_repeat_ambiguous",
                "left_junction_only_divergent_span",
                "locus_truncated_at_contig_end"}


def delta_bands():
    """event -> band, replicating the acceptance script's join exactly."""
    out = {}
    p = W + "/audits/delta_nonzero_classified.tsv"
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cl, d = f[ix["class"]], int(f[ix["delta"]])
            unexp = ("unexplained" in cl) or ("DESCRIPTIVE" in cl)
            if not unexp:
                out[f[ix["db_id"]]] = "candidate_cause"
            elif abs(d) <= 3:
                out[f[ix["db_id"]]] = "unexplained_<=3bp"
            elif abs(d) <= 50:
                out[f[ix["db_id"]]] = "unexplained_4-50bp"
            else:
                out[f[ix["db_id"]]] = "unexplained_>50bp"
    return out


def audit_a(d):
    """occurrence_id -> verdict, merged across shards."""
    out = {}
    files = sorted(glob.glob(d + "/auditA_occ_*.tsv"))
    if not files:
        raise SystemExit("no Audit A shard output in " + d)
    for p in files:
        with open(p) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out[f[ix["occurrence_id"]]] = f[ix["verdict"]]
    print("   Audit A shards merged: %d files, %d flagged occurrences"
          % (len(files), len(out)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-dir", default=W + "/auditA_occ/full")
    ap.add_argument("--out", default=V2)
    ap.add_argument("--exclude-4-50", action="store_true",
                    help="also drop the unexplained 4-50bp band")
    a = ap.parse_args()

    bands = delta_bands()
    av = audit_a(a.audit_dir)

    src = V2 + "/occurrences.tsv"
    fl = a.out + "/occurrences_flagged.tsv"
    keep = set()
    n = collections.Counter()
    xt = collections.Counter()
    ev_keep, ev_off = set(), set()
    bp_keep = collections.Counter()
    with open(src) as fh, open(fl, "w") as o:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        o.write("\t".join(cols + ["delta_band", "auditA_verdict",
                                  "exclusion_reason", "analysis_set"]) + "\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            occ, db = f[ix["occurrence_id"]], f[ix["db_id"]]
            off = f[ix["official_fasta"]] == "yes"
            band = bands.get(db, "delta_zero")
            ver = av.get(occ, "not_examined" if not off else "locus_empty")
            reasons = []
            if ver in A_CONTRADICT:
                reasons.append("auditA_element_at_reference_locus")
            elif ver in A_UNRESOLVED:
                reasons.append("auditA_unresolved_at_reference_locus")
            if band == "unexplained_>50bp":
                reasons.append("unexplained_length_difference_gt50bp")
            elif band == "unexplained_4-50bp" and a.exclude_4_50:
                reasons.append("unexplained_length_difference_4_50bp")
            if not off:
                aset = "withheld_by_inference_rules"
            elif reasons:
                aset = "flagged_for_review"
            else:
                aset = "conservative"
            if off:
                n[aset] += 1
                ev_off.add(db)
                xt[(band, "flagged" if reasons else "conservative")] += 1
                for r in reasons:
                    n["reason:" + r] += 1
                ina = any(r.startswith("auditA") for r in reasons)
                ind = any(r.startswith("unexplained") for r in reasons)
                if ina and ind:
                    n["overlap:auditA_AND_length"] += 1
                elif ina:
                    n["overlap:auditA_only"] += 1
                elif ind:
                    n["overlap:length_only"] += 1
                if not reasons:
                    keep.add(occ)
                    ev_keep.add(db)
                    bp_keep[f[ix["breakpoint_status"]]] += 1
            o.write("\t".join(f + [band, ver, ";".join(reasons) or ".",
                                   aset]) + "\n")

    # ---- conservative FASTA, filtered from the full one by header token ----
    src_fa, dst_fa = V2 + "/inferred_empty.fna", a.out + "/conservative.fna"
    w, nfa = False, 0
    with open(src_fa) as fh, open(dst_fa, "w") as o:
        for line in fh:
            if line[0] == ">":
                w = line[1:].split()[0] in keep
                if w:
                    nfa += 1
            if w:
                o.write(line)

    print("\nFULL v2 (unchanged, now flagged): %s" % fl)
    print("CONSERVATIVE SUBSET:              %s" % dst_fa)
    print("\n   official occurrences            %7d  over %d events"
          % (n["conservative"] + n["flagged_for_review"], len(ev_off)))
    print("   conservative                    %7d  over %d events"
          % (n["conservative"], len(ev_keep)))
    print("   flagged for review              %7d" % n["flagged_for_review"])
    print("\n   exclusion reasons (OVERLAPPING -- these sum to more than the"
          " %d excluded rows)" % n["flagged_for_review"])
    for k, v in sorted(n.items()):
        if k.startswith("reason:"):
            print("      %-44s %7d" % (k[7:], v))
    print("\n   the two exclusion sets, as a partition of the excluded rows")
    for k in ("overlap:auditA_only", "overlap:length_only",
              "overlap:auditA_AND_length"):
        print("      %-44s %7d" % (k[8:], n[k]))
    print("      %-44s %7d" % ("UNION (= flagged_for_review)",
                               n["flagged_for_review"]))
    print("\n   breakpoint_status within the conservative subset")
    tot = max(1, n["conservative"])
    for k, v in sorted(bp_keep.items()):
        print("      %-28s %7d  %5.2f%%" % (k, v, 100.0 * v / tot))
    print("\n   delta band x analysis set (official rows)")
    print("      %-24s %12s %12s" % ("band", "conservative", "flagged"))
    for b in ("delta_zero", "candidate_cause", "unexplained_<=3bp",
              "unexplained_4-50bp", "unexplained_>50bp"):
        c, g = xt[(b, "conservative")], xt[(b, "flagged")]
        if c or g:
            print("      %-24s %12d %12d" % (b, c, g))
    assert nfa == len(keep), "FASTA %d != kept %d" % (nfa, len(keep))
    print("\n   conservative FASTA records %d == kept occurrence ids %d  OK"
          % (nfa, len(keep)))


if __name__ == "__main__":
    main()

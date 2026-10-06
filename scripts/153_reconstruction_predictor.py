#!/usr/bin/env python3
"""When is an excision-reconstructed empty site CORRECT? Find a usable rule.

152_ measured 81.0% flank-120 accuracy overall, with a wide spread by stratum
(0.0% where target bases were lost, 42.2% at overlap 0, 85.2% at overlap 4-15).
A corpus cannot ship at 81% because a wrong flank is indistinguishable from a
right one downstream. The question is whether some condition identifies the
correct ones precisely enough to be worth labelling.

SPLIT-SAMPLE, because this project has twice had a post-hoc rule fail by
misspecification. Panels are partitioned by hash into DERIVE and TEST halves.
Every rate quoted for a candidate rule is its TEST-half rate; the derive half
only proposes. A rule that looks good on derive and collapses on test is a
rule that was fitted to noise, and that outcome is reported rather than hidden.

THE ASSUMED RELATIONSHIP, stated before looking (per the pre-registration
lesson -- a threshold alone is not a hypothesis):

  Excision reproduces the empty site when the element's BOUNDARY is
  unambiguous and no target sequence was destroyed. So correctness should
  increase with decomposition confidence and class stability, and should fall
  when the decomposition had equivalent alternatives, when bases were lost, or
  when the junction's microhomology makes the boundary slide.

  Features are therefore tested in that order of prior plausibility, and a
  feature that predicts well WITHOUT fitting that story is treated as
  suspicious rather than as a discovery.

Reported per candidate rule: precision (share of kept loci whose flank is
right) and coverage (share of all loci kept). A rule is only useful if both
are high -- 99% precision over 2% of loci adds nothing.
"""
import collections
import csv
import glob
import hashlib
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI

HALF = 60
dirs = []
for pat in sys.argv[1].split(","):
    dirs.extend(sorted(glob.glob(pat)))

rows = []
for d in dirs:
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)):
        continue
    panel = LI.panel_label(d)
    half = int(hashlib.md5(panel.encode()).hexdigest()[:8], 16) % 2
    aseq = LI.read_alleles(ap)
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        av = aseq.get(r["locus_id"], {})
        lng, sht = av.get(r["longest_allele"], ""), av.get(r["shortest_allele"], "")
        if not lng or not sht:
            continue
        i0, _ = LI.locate_insert(r, lng)
        if i0 < 0:
            continue
        ilen = int(r["inserted_len"])
        ins = lng[i0:i0 + ilen]
        j = lng.find(ins)
        if j < 0:
            continue
        multi = lng.find(ins, j + 1) >= 0
        rec = lng[:j] + lng[j + ilen:]
        if j < HALF or len(sht) - j < HALF or len(rec) - j < HALF:
            continue
        ok = sht[j - HALF:j + HALF] == rec[j - HALF:j + HALF]

        def num(k, d=0.0):
            try:
                return float(r.get(k) or d)
            except ValueError:
                return d
        rows.append({
            "half": half, "ok": ok,
            "lost": num("target_bases_lost"),
            "ov": num("junction_ambiguity_bp"),
            "method": r.get("decomposition_method", ""),
            "dconf": num("decomposition_confidence"),
            "nequiv": num("n_equivalent_decompositions", 1),
            "stable": r.get("class_stable", ""),
            "place": r.get("placement_status", ""),
            "pmargin": num("placement_score_margin"),
            "hident": num("homologous_identity"),
            "ilen": num("inserted_len"),
            "nall": num("n_alleles"),
            "multi": multi})

print("loci scored %d   derive %d / test %d"
      % (len(rows), sum(1 for x in rows if x["half"] == 0),
         sum(1 for x in rows if x["half"] == 1)))
base = sum(1 for x in rows if x["ok"]) / max(1, len(rows))
print("baseline flank120 accuracy: %.1f%%" % (100 * base))
print()


def ev(rule, label):
    for h, name in ((0, "derive"), (1, "TEST")):
        sub = [x for x in rows if x["half"] == h]
        kept = [x for x in sub if rule(x)]
        if not kept:
            print("  %-46s %-7s kept 0" % (label if h == 0 else "", name))
            continue
        prec = 100.0 * sum(1 for x in kept if x["ok"]) / len(kept)
        cov = 100.0 * len(kept) / max(1, len(sub))
        print("  %-46s %-7s precision %5.1f%%  coverage %5.1f%%  n=%d"
              % (label if h == 0 else "", name, prec, cov, len(kept)))


print("single conditions")
ev(lambda x: True, "no filter")
ev(lambda x: x["lost"] == 0, "lost == 0")
ev(lambda x: 4 <= x["ov"] <= 15, "overlap 4-15")
ev(lambda x: not x["multi"], "insert occurs once in the carrier")
ev(lambda x: x["method"] == "exact", "exact decomposition path")
ev(lambda x: x["nequiv"] <= 1, "no equivalent decompositions")
ev(lambda x: x["stable"] == "True", "class_stable")
ev(lambda x: x["place"] != "ambiguous", "placement not ambiguous")
ev(lambda x: x["dconf"] >= 0.99, "decomposition_confidence >= 0.99")
print()
print("combinations, in the order the stated hypothesis predicts")
ev(lambda x: x["lost"] == 0 and not x["multi"], "lost==0 AND single-hit")
ev(lambda x: x["lost"] == 0 and not x["multi"] and x["nequiv"] <= 1,
   "  + no equivalent decompositions")
ev(lambda x: x["lost"] == 0 and not x["multi"] and x["nequiv"] <= 1
   and x["stable"] == "True", "  + class_stable")
ev(lambda x: x["lost"] == 0 and not x["multi"] and x["nequiv"] <= 1
   and x["stable"] == "True" and x["place"] != "ambiguous", "  + placement ok")
ev(lambda x: x["lost"] == 0 and not x["multi"] and x["nequiv"] <= 1
   and x["stable"] == "True" and x["place"] != "ambiguous" and 4 <= x["ov"] <= 15,
   "  + overlap 4-15")
print()
print("A rule is useful only if precision AND coverage are high on TEST.")

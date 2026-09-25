#!/usr/bin/env python3
"""Benchmark 2 -- the FNA-first replacement decision. TIERED, not pooled.

The 70 read-audited events are NOT one gold standard. They stratify:

    T1  nucleotide-precise           n=24   <- the only real precision gold
    T2  near-precise                 n=5
    T3  located but not nt-resolved  n=41   <- detection gold only

Pooling them would let 41 loosely-resolved events dominate a precision number
they cannot support. So each tier answers a different question:

  Primary precision    T1 (24)      can FNA replace reads for exact flanks?
  Secondary precision  T1+T2 (29)   is it usable at +/-5 bp?
  Detection recall     all 70       does FNA find known events at all?
  Architecture         asymmetric   does FNA recover p_L != p_R?

BOTH JUNCTIONS, NEVER A SINGLE POS. An insertion has a left and a right
junction, and comparing against one Sniffles POS hides asymmetry:

    dL = pL_fna - pL_read        dR = pR_fna - pR_read

The headline metric requires BOTH to be right:

    max(|dL|, |dR|) <= k

One-sided agreement is reported separately, for contrast only -- it is not a
pass. Events like dL=+59/dR=0 or dL=-56/dR=0 are the whole reason: a
single-POS comparison would score them as a near miss or a hit depending on
which end it happened to use.

PRE-REGISTERED DECISION RULE (fixed before looking at results, --decide to
override). On the T1 tier:

    detection recall            >= 90%
    both junctions within +/-5  >= 80%
    inserted-seq identity       >= 95%
    inserted-seq coverage       >= 90%

All four met -> FNA-first is sufficient for large-scale discovery, and reads
drop to gold validation of a few representative families and unusual junction
architectures.

GOLD TSV columns (override with --gold-cols):
    event_id, sample, chrom, pos_l, pos_r, tier, architecture,
    insert_len, insert_seq
  pos_r may be absent -- then pos_r is taken as pos_l and the architecture test
  is skipped for that event, which is recorded rather than silently ignored.

TEST TSV: armB_events.tsv or universal_events.tsv (needs junction_L/junction_R).
"""
import argparse
import csv
import os
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402

DEFAULT_GOLD = "event_id,sample,chrom,pos_l,pos_r,tier,architecture,insert_len,insert_seq"
DEFAULT_TEST = ("event_id,empty_ref_genome,empty_ref_contig,junction_L,junction_R,"
                "architecture,insert_len_direct,")


def load(path, spec, keys):
    cols = dict(zip(keys, spec.split(",")))
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            rec = {}
            for k, v in cols.items():
                rec[k] = r.get(v, "") if v else ""
            rec["_raw"] = r
            out.append(rec)
    return out


def as_int(v, default=None):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def read_fasta(path):
    seqs, name = {}, None
    if not path or not os.path.exists(path):
        return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


# ------------------------------------------------------- sequence comparison
def kmers(s, k=15):
    s = s.upper()
    return {s[i:i + k] for i in range(len(s) - k + 1) if "N" not in s[i:i + k]}


def seq_compare(a, b, minimap2=None):
    """Return (identity_pct, coverage_pct_of_gold). minimap2 when available,
    k-mer containment otherwise -- the fallback is labelled in the output."""
    if not a or not b:
        return float("nan"), float("nan"), "none"
    if minimap2:
        try:
            with tempfile.TemporaryDirectory() as td:
                fa, fb = os.path.join(td, "a.fa"), os.path.join(td, "b.fa")
                open(fa, "w").write(">gold\n%s\n" % a)
                open(fb, "w").write(">test\n%s\n" % b)
                p = subprocess.run([minimap2, "-c", "-x", "asm20",
                                    "--secondary=no", fa, fb],
                                   capture_output=True, text=True, check=True)
                best = None
                for line in p.stdout.splitlines():
                    f = line.split("\t")
                    if len(f) < 12:
                        continue
                    nm, al = int(f[9]), int(f[10])
                    cov = (int(f[8]) - int(f[7])) / len(a) * 100.0
                    ident = 100.0 * nm / al if al else 0.0
                    if best is None or cov > best[1]:
                        best = (ident, cov)
                if best:
                    return best[0], min(100.0, best[1]), "minimap2"
        except Exception:
            pass
    ka, kb = kmers(a), kmers(b)
    if not ka or not kb:
        return float("nan"), float("nan"), "none"
    inter = len(ka & kb)
    return 100.0 * inter / len(kb), 100.0 * inter / len(ka), "kmer"


# --------------------------------------------------------------- matching
def build_index(test):
    idx = defaultdict(list)
    for t in test:
        idx[t.get("chrom", "")].append(t)
        idx[""].append(t)
    return idx


def match(g, idx, window, gseq, tseq, min_fallback_ident):
    """Pair a gold event with a test call.

    1. coordinate match -- same contig, both junctions inside `window`.
    2. cross-frame fallback -- SAME INSERTED SEQUENCE. Matching on insert
       length alone is unsafe and was measured to be: every IS1 copy is 768 bp,
       so a length-only fallback invents a match for every unmatched event and
       reports 100% recall. The fallback therefore requires reciprocal k-mer
       containment >= min_fallback_ident, and is skipped entirely when no
       sequence is available for either side.
    """
    gl = as_int(g["pos_l"])
    gr = as_int(g["pos_r"], gl)
    best, bscore, mode = None, None, "none"
    for t in idx.get(g.get("chrom", ""), ()):
        tl, tr = as_int(t["junction_L"]), as_int(t["junction_R"])
        if tl is None or gl is None:
            continue
        d = max(abs(tl - gl), abs((tr if tr is not None else tl) - gr))
        if d <= window and (bscore is None or d < bscore):
            best, bscore, mode = t, d, "coord"
    if best is not None:
        return best, mode

    a = g["insert_seq"] or gseq.get(g["event_id"], "")
    if not a:
        return None, "none_no_seq"
    ka = kmers(a)
    if not ka:
        return None, "none_no_seq"
    bestc = 0.0
    for t in idx.get("", ()):
        b = t["insert_seq"] or tseq.get(t["event_id"], "")
        if not b:
            continue
        kb = kmers(b)
        if not kb:
            continue
        inter = len(ka & kb)
        c = min(inter / len(ka), inter / len(kb)) * 100.0
        if c >= min_fallback_ident and c > bestc:
            best, bestc, mode = t, c, "seq_crossframe"
    return best, (mode if best is not None else "none")


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--gold-cols", default=DEFAULT_GOLD)
    ap.add_argument("--test-cols", default=DEFAULT_TEST)
    ap.add_argument("--gold-inserts", default=None, help="FASTA keyed by gold event_id")
    ap.add_argument("--test-inserts", default=None, help="FASTA keyed by test event_id")
    ap.add_argument("--window", type=int, default=2000,
                    help="coordinate search radius for pairing (default 2000)")
    ap.add_argument("--fallback-min-identity", type=float, default=95.0,
                    help="min reciprocal k-mer containment %% for the cross-frame "
                         "fallback; length-only matching is deliberately not offered")
    ap.add_argument("--primary-tier", default="T1")
    ap.add_argument("--secondary-tiers", default="T1,T2")
    ap.add_argument("--minimap2", default="minimap2")
    # pre-registered thresholds
    ap.add_argument("--decide-recall", type=float, default=90.0)
    ap.add_argument("--decide-both-bp", type=int, default=5)
    ap.add_argument("--decide-both-pct", type=float, default=80.0)
    ap.add_argument("--decide-identity", type=float, default=95.0)
    ap.add_argument("--decide-coverage", type=float, default=90.0)
    args = ap.parse_args()

    gkeys = ["event_id", "sample", "chrom", "pos_l", "pos_r", "tier",
             "architecture", "insert_len", "insert_seq"]
    tkeys = ["event_id", "sample", "chrom", "junction_L", "junction_R",
             "architecture", "insert_len", "insert_seq"]
    gold = load(args.gold, args.gold_cols, gkeys)
    test = load(args.test, args.test_cols, tkeys)
    gseq = read_fasta(args.gold_inserts)
    tseq = read_fasta(args.test_inserts)
    mm2 = args.minimap2 if subprocess.run(["which", args.minimap2],
                                          capture_output=True).returncode == 0 else None
    log("gold=%d  test=%d  minimap2=%s" % (len(gold), len(test), bool(mm2)))

    idx = build_index(test)
    rows = []
    for g in gold:
        gl, gr = as_int(g["pos_l"]), as_int(g["pos_r"])
        one_sided_gold = gr is None
        if one_sided_gold:
            gr = gl
        t, mode = match(g, idx, args.window, gseq, tseq,
                        args.fallback_min_identity)
        if t is None:
            rows.append(dict(g, matched=0, mode="none", dL="", dR="", dmax="",
                             ident="", cov="", cmp_mode="", test_event="",
                             test_span="", gold_span=(gr - gl) if gl is not None else "",
                             one_sided_gold=int(one_sided_gold)))
            continue
        tl, tr = as_int(t["junction_L"]), as_int(t["junction_R"])
        dL = tl - gl if (tl is not None and gl is not None) else None
        dR = tr - gr if (tr is not None and gr is not None) else None
        dmax = max(abs(dL), abs(dR)) if (dL is not None and dR is not None) else None
        a = g["insert_seq"] or gseq.get(g["event_id"], "")
        b = t["insert_seq"] or tseq.get(t["event_id"], "")
        ident, cov, cmode = seq_compare(a, b, mm2)
        rows.append(dict(g, matched=1, mode=mode,
                         dL="" if dL is None else dL, dR="" if dR is None else dR,
                         dmax="" if dmax is None else dmax,
                         ident="" if ident != ident else round(ident, 2),
                         cov="" if cov != cov else round(cov, 2),
                         cmp_mode=cmode, test_event=t["event_id"],
                         test_span="" if (tl is None or tr is None) else tr - tl,
                         gold_span=gr - gl, one_sided_gold=int(one_sided_gold)))

    # ------------------------------------------------------------ per-event
    with open(args.out + "_per_event.tsv", "w") as fh:
        cols = ["event_id", "sample", "tier", "architecture", "matched", "mode",
                "dL", "dR", "dmax", "gold_span", "test_span", "one_sided_gold",
                "ident", "cov", "cmp_mode", "test_event"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ------------------------------------------------------------- summary
    def tier_stats(sel, label):
        n = len(sel)
        m = [r for r in sel if r["matched"]]
        both = [r for r in m if r["dmax"] != ""]
        st = {"tier": label, "n_gold": n, "n_matched": len(m),
              "recall_pct": round(100.0 * len(m) / n, 1) if n else float("nan")}
        for k in (1, 2, 5, 10):
            st["both_within_%d_pct" % k] = (
                round(100.0 * sum(1 for r in both if abs(r["dmax"]) <= k) / n, 1)
                if n else float("nan"))
            # one-sided, for contrast only -- not a pass criterion
            st["oneside_within_%d_pct" % k] = (
                round(100.0 * sum(1 for r in m
                                  if r["dL"] != "" and r["dR"] != ""
                                  and min(abs(r["dL"]), abs(r["dR"])) <= k) / n, 1)
                if n else float("nan"))
        for key, name in (("ident", "identity"), ("cov", "coverage")):
            vals = [float(r[key]) for r in m if r[key] != ""]
            st["insert_%s_median" % name] = round(statistics.median(vals), 2) if vals else float("nan")
        dl = [abs(r["dL"]) for r in m if r["dL"] != ""]
        dr = [abs(r["dR"]) for r in m if r["dR"] != ""]
        st["median_absdL"] = round(statistics.median(dl), 1) if dl else float("nan")
        st["median_absdR"] = round(statistics.median(dr), 1) if dr else float("nan")
        ln = [abs(as_int(r["insert_len"], 0) - 0) for r in m if r["insert_len"]]
        st["n_len_reported"] = len(ln)
        st["n_coord_matched"] = sum(1 for r in m if r["mode"] == "coord")
        st["n_crossframe_matched"] = sum(1 for r in m if r["mode"] == "seq_crossframe")
        return st

    by_tier = defaultdict(list)
    for r in rows:
        by_tier[r["tier"] or "UNSPECIFIED"].append(r)
    prim = args.primary_tier
    sec = args.secondary_tiers.split(",")
    stats = [tier_stats(by_tier.get(prim, []), "PRIMARY(%s)" % prim),
             tier_stats([r for r in rows if r["tier"] in sec],
                        "SECONDARY(%s)" % "+".join(sec)),
             tier_stats(rows, "DETECTION(all)")]
    for t in sorted(by_tier):
        stats.append(tier_stats(by_tier[t], "tier_%s" % t))

    with open(args.out + "_summary.tsv", "w") as fh:
        keys = list(stats[0])
        fh.write("\t".join(keys) + "\n")
        for st in stats:
            fh.write("\t".join(str(st.get(k, "")) for k in keys) + "\n")

    # -------------------------------------------------------- architecture
    asym_gold = [r for r in rows
                 if (not r["one_sided_gold"]) and abs(as_int(r["gold_span"], 0)) > 2]
    arch = {"n_asymmetric_gold": len(asym_gold),
            "n_matched": sum(1 for r in asym_gold if r["matched"]),
            "n_test_also_asymmetric": sum(
                1 for r in asym_gold if r["matched"] and r["test_span"] != ""
                and abs(int(r["test_span"])) > 2),
            "n_span_within_5bp": sum(
                1 for r in asym_gold if r["matched"] and r["test_span"] != ""
                and abs(int(r["test_span"]) - as_int(r["gold_span"], 0)) <= 5)}
    with open(args.out + "_architecture.tsv", "w") as fh:
        for k, v in arch.items():
            fh.write("%s\t%s\n" % (k, v))

    # ----------------------------------------------------------- report
    def p(msg=""):
        sys.stderr.write(msg + "\n")
    p("=" * 74)
    p("TIERED READ-GOLD BENCHMARK")
    p("=" * 74)
    hdr = ("%-22s %6s %8s %8s %9s %9s %9s %9s" %
           ("tier", "n", "recall%", "+/-1bp", "+/-2bp", "+/-5bp", "ident%", "cov%"))
    p(hdr)
    p("-" * 74)
    for st in stats:
        p("%-22s %6d %8s %8s %9s %9s %9s %9s" % (
            st["tier"], st["n_gold"], st["recall_pct"],
            st["both_within_1_pct"], st["both_within_2_pct"],
            st["both_within_5_pct"], st["insert_identity_median"],
            st["insert_coverage_median"]))
    p("-" * 74)
    pr = stats[0]
    p("+/- columns require BOTH junctions: max(|dL|,|dR|) <= k")
    p("match provenance on %s: %d by coordinate, %d by cross-frame sequence"
      % (prim, pr["n_coord_matched"], pr["n_crossframe_matched"]))
    p("  (cross-frame matches cannot support a junction delta -- if that count is")
    p("   non-zero, lift the gold into the test coordinate frame before deciding)")
    p("one-sided-only within +/-5bp on %s: %s%%  (NOT a pass criterion,"
      % (prim, pr["oneside_within_5_pct"]))
    p("  shown because a single-POS comparison would report this number)")
    p()
    p("ARCHITECTURE (asymmetric gold events, |gold_span| > 2 bp)")
    for k, v in arch.items():
        p("  %-28s %s" % (k, v))
    p()
    p("=" * 74)
    p("PRE-REGISTERED DECISION on %s" % prim)
    p("=" * 74)
    checks = [
        ("detection recall", pr["recall_pct"], args.decide_recall),
        ("both junctions within +/-%d bp" % args.decide_both_bp,
         pr["both_within_%d_pct" % args.decide_both_bp], args.decide_both_pct),
        ("inserted-seq identity", pr["insert_identity_median"], args.decide_identity),
        ("inserted-seq coverage", pr["insert_coverage_median"], args.decide_coverage),
    ]
    allpass = True
    for name, got, need in checks:
        ok = (got == got) and got >= need
        allpass &= bool(ok)
        p("  %-34s %8s  >= %-6s  %s"
          % (name, got, need, "PASS" if ok else "FAIL"))
    p()
    p("  VERDICT: %s" % ("FNA-FIRST SUFFICIENT for large-scale discovery; reads "
                         "become gold validation of representative families "
                         "and unusual junction architectures."
                         if allpass else
                         "NOT YET -- keep reads in the main path for the failing "
                         "criteria above."))
    p("=" * 74)
    log("wrote %s_{per_event,summary,architecture}.tsv" % args.out)
    return 0 if allpass else 1


if __name__ == "__main__":
    sys.exit(main())

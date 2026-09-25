#!/usr/bin/env python3
"""Layer 1 + Layer 2 validity on a species with NO frozen baseline. Self-contained.

`73_layer2_regression.py` is a REGRESSION test: its gates are frozen against
E. coli values (234 loci, 91.7% exact, 96.2% revcomp invariance) and verify that
Layer 2 still does what it did on 2026-08-28. Those numbers are species-specific,
so they say nothing about whether a new species' output is correct. This module
is the check that can run anywhere, needs no truth set, and no external tool.

CHECK 1 -- RECONSTRUCTION IDENTITY. The strongest self-contained test in the
    pipeline: `pre[:offset] + insert + pre[offset:]` must regenerate the OBSERVED
    derived allele. It ties the decomposition and the offset together, so if it
    fails everything downstream is void. E. coli reference: median 1.0000,
    min 0.9915.

CHECK 2 -- ARM B <-> LAYER 2 CONCORDANCE, and the ONLY correctness check
    available for Layer 1. Arm B measures insert size from the graph bubble;
    Layer 2 measures it from sequence. Two independent routes to one quantity.
    Layer 1's three assertions catch empty input; they do not verify that a
    bubble is real. **Sequence comparison is orientation-agnostic**: Arm B stores
    in backbone orientation and Layer 2 in the carrier's, and a forward-only
    comparison returns ~0.50 -- the score of two unrelated sequences. That bug
    once nearly produced a false 19.2% error rate.

CHECK 3 -- CATASTROPHIC PLACEMENT. Anchors landing on different copies of a
    repeat once fabricated a 42,367 bp allele. The guard is the interval cap, so
    count placements AT the cap and inspect the inserted_len range. A
    repeat-richer species is exactly the condition that produced the original
    failure.

CHECK 4 -- WHAT DID NOT DECOMPOSE. Classify the non-decomposable loci by Layer
    2's own event_class. `length_polymorphism` is the known limitation -- alleles
    align end-to-end at ~85%, so short-plus-insert is the wrong model. But
    `substitution_only` or `monomorphic` in quantity means Arm B called bubbles
    where there is no insertion, which is a LAYER 1 PRECISION problem.

    101_layer12_validity.py --recon <dirs> --armb <dirs> --out validity
"""
import argparse
import collections
import csv
import glob
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_alleles(path):
    out, cur = collections.defaultdict(dict), None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "").upper() for a, v in d.items()}
            for k, d in out.items()}


def read_fasta(path):
    out, cur = {}, None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}



def expand(paths):
    """Expand glob patterns that survived the shell.

    `slurm/analysis.sh` sets `set -f` so patterns reach the script intact
    (without it the submitting shell expands them and argparse sees dozens of
    positional args). The cost is that any script taking a directory list MUST
    glob for itself -- and when a chained job is submitted BEFORE its input
    directories exist, the pattern is all it ever gets. 100_event_dedup.py once
    read 0 panels this way and the whole A. baumannii chain stalled.
    """
    out = []
    for p in paths:
        out.extend(sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p])
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--armb", nargs="+", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--interval-cap", type=int, default=7000)
    ap.add_argument("--len-tol", type=int, default=10)
    args = ap.parse_args()
    args.recon = expand(args.recon)
    args.armb = expand(args.armb)

    from Bio.Align import PairwiseAligner
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-6, extend_gap_score=-0.5)

    def ident(a, b):
        if a == b:
            return 1.0
        if not a or not b:
            return float("nan")
        try:
            x = al.align(a, b)[0]
        except Exception:
            return float("nan")
        A, B = x.aligned
        m = n = 0
        for (s0, s1), (t0, t1) in zip(A, B):
            m += sum(1 for u, v in zip(a[s0:s1], b[t0:t1]) if u == v)
            n += s1 - s0
        return m / n if n else float("nan")

    loci, aseq = {}, {}
    for d in args.recon:
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            loci[r["locus_id"]] = r
        aseq.update(read_alleles(os.path.join(d, "alleles.fna")))
    log("%d Layer-2 loci" % len(loci))

    # ---------------- CHECK 4 first: it is free ---------------------------
    dec, nodec = [], collections.Counter()
    for lid, r in loci.items():
        if (r["event_class"] in ("insertion_target_retained", "replacement")
                and int(r["inserted_len"] or 0) > 0):
            dec.append(lid)
        else:
            nodec[r["event_class"]] += 1

    # ---------------- CHECK 1: reconstruction identity --------------------
    ids, exact_hit, fail = [], 0, []
    for lid in dec:
        r = loci[lid]
        av = aseq.get(lid, {})
        s_, l_ = av.get(r["shortest_allele"], ""), av.get(r["longest_allele"], "")
        if not s_ or not l_:
            continue
        off, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
        ins = l_[off:off + ilen]
        rebuilt = s_[:off] + ins + s_[off:]
        if rebuilt == l_:
            ids.append(1.0)
            exact_hit += 1
            continue
        i = ident(rebuilt, l_)
        ids.append(i)
        if i == i and i < 0.99:
            fail.append((lid, i, r.get("decomposition_method")))

    # ---------------- CHECK 3: catastrophic placement ---------------------
    at_cap = [lid for lid in dec
              if int(loci[lid]["inserted_len"] or 0) >= args.interval_cap - 50]
    lens = [int(loci[lid]["inserted_len"]) for lid in dec]
    pstat = collections.Counter(loci[l].get("placement_status", ".") for l in dec)
    pconf = collections.Counter(loci[l].get("placement_confidence", ".") for l in dec)

    # ---------------- CHECK 2: Arm B <-> Layer 2 --------------------------
    ev, ins_fa = {}, {}
    for d in args.armb:
        p = os.path.join(d, "armB_events.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                ev[r["event_id"]] = r
        ins_fa.update(read_fasta(os.path.join(d, "armB_inserts.fna")))
    matched = [lid for lid in dec if lid in ev]
    lag = []
    seqid, seq_rc = [], 0
    for lid in matched:
        a = int(ev[lid].get("insert_size") or ev[lid].get("insert_len_direct") or 0)
        b = int(loci[lid]["inserted_len"])
        if a > 0:
            lag.append(b - a)
        av = aseq.get(lid, {})
        l_ = av.get(loci[lid]["longest_allele"], "")
        off, ilen = int(loci[lid]["lcp_bp"]), int(loci[lid]["inserted_len"])
        mine = l_[off:off + ilen]
        theirs = ins_fa.get(lid, "")
        if mine and theirs:
            f = ident(mine, theirs)
            rv = ident(mine, rc(theirs))
            # ORIENTATION-AGNOSTIC: forward-only scores ~0.50 on the same seq
            if rv == rv and (f != f or rv > f):
                seq_rc += 1
                seqid.append(rv)
            else:
                seqid.append(f)

    def qs(a, p):
        a = sorted(a)
        return a[min(len(a) - 1, int(len(a) * p))] if a else float("nan")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_validity.tsv", "w") as fh:
        fh.write("check\tmetric\tvalue\n")
        fh.write("recon_identity\tn\t%d\n" % len(ids))
        if ids:
            fh.write("recon_identity\tmedian\t%.4f\n" % statistics.median(ids))
            fh.write("recon_identity\tmin\t%.4f\n" % min(ids))
            fh.write("recon_identity\texact_bytewise\t%d\n" % exact_hit)
            fh.write("recon_identity\tbelow_0.99\t%d\n" % len(fail))
        fh.write("armb_concordance\tn_matched\t%d\n" % len(matched))
        if lag:
            fh.write("armb_concordance\tlen_exact\t%d\n"
                     % sum(1 for x in lag if x == 0))
            fh.write("armb_concordance\tlen_within_tol\t%d\n"
                     % sum(1 for x in lag if abs(x) <= args.len_tol))
        if seqid:
            fh.write("armb_concordance\tseq_identity_median\t%.4f\n"
                     % statistics.median(seqid))
            fh.write("armb_concordance\trevcomp_orientation\t%d\n" % seq_rc)
        fh.write("placement\tat_interval_cap\t%d\n" % len(at_cap))
        fh.write("placement\tinserted_len_max\t%d\n" % (max(lens) if lens else 0))
        for k, v in nodec.items():
            fh.write("not_decomposed\t%s\t%d\n" % (k, v))

    log("=" * 82)
    log("LAYER 1 + LAYER 2 VALIDITY  (no frozen baseline, no truth set)")
    log("")
    log("  CHECK 1 -- RECONSTRUCTION IDENTITY   pre[:off] + insert + pre[off:]")
    if ids:
        log("    n                       %6d" % len(ids))
        log("    byte-identical          %6d   %5.1f%%"
            % (exact_hit, 100.0 * exact_hit / len(ids)))
        log("    median identity         %.4f" % statistics.median(ids))
        log("    q1 / min                %.4f / %.4f" % (qs(ids, .25), min(ids)))
        log("    below 0.99              %6d   %5.1f%%"
            % (len(fail), 100.0 * len(fail) / len(ids)))
        if fail:
            log("    worst: " + "  ".join("%s=%.3f(%s)" % f for f in
                                          sorted(fail, key=lambda x: x[1])[:4]))
    log("")
    log("  CHECK 2 -- ARM B <-> LAYER 2   (two independent routes to insert size)")
    log("    loci matched to an Arm B event   %6d of %d" % (len(matched), len(dec)))
    if lag:
        log("    insert length EXACT              %6d   %5.1f%%"
            % (sum(1 for x in lag if x == 0), 100.0 * sum(1 for x in lag if x == 0) / len(lag)))
        log("    within +/-%d bp                   %6d   %5.1f%%"
            % (args.len_tol, sum(1 for x in lag if abs(x) <= args.len_tol),
               100.0 * sum(1 for x in lag if abs(x) <= args.len_tol) / len(lag)))
        log("    Layer2 - ArmB: median %d  q1 %d  q3 %d"
            % (statistics.median(lag), qs(lag, .25), qs(lag, .75)))
    if seqid:
        log("    insert SEQUENCE identity (orientation-agnostic): median %.4f"
            % statistics.median(seqid))
        log("    scored on the reverse complement %6d   %5.1f%%  <- expected, "
            "not an error" % (seq_rc, 100.0 * seq_rc / len(seqid)))
    log("")
    log("  CHECK 3 -- CATASTROPHIC PLACEMENT")
    log("    at the %d bp interval cap        %6d" % (args.interval_cap, len(at_cap)))
    if lens:
        log("    inserted_len  min %d  median %d  max %d"
            % (min(lens), statistics.median(lens), max(lens)))
    log("    placement_status: " + "  ".join("%s=%d" % kv for kv in pstat.most_common(4)))
    log("    placement_confidence: " + "  ".join("%s=%d" % kv for kv in pconf.most_common(4)))
    log("")
    log("  CHECK 4 -- WHAT DID NOT DECOMPOSE   (a LAYER 1 precision question)")
    tot_nd = sum(nodec.values())
    for k, v in nodec.most_common():
        tag = ("known limitation" if k == "length_polymorphism"
               else "*** Arm B called a bubble with no insertion ***"
               if k in ("substitution_only", "monomorphic") else "")
        log("    %-26s %6d   %5.1f%%  %s"
            % (k, v, 100.0 * v / tot_nd if tot_nd else 0, tag))
    log("=" * 82)
    log("wrote %s_validity.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

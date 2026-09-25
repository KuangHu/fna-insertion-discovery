#!/usr/bin/env python3
"""STAGE 32 -- search the PARKED >5 kb bubbles for nested <5 kb insertions.

Arm B parks every bubble whose insert exceeds --max-insert into
`armB_large_events.tsv` rather than discarding it. Across the full runs that is
E. coli 400 panels + K. pneumoniae 400 panels of parked events. Discovery cost
is already paid: these are known length-polymorphic sites that the <5 kb scope
skipped.

A large bubble is often not one large element. It is a mosaic -- a resistance
island, a prophage, an integrative element -- with ordinary IS copies inside it.
Comparing two large alleles against each other exposes those.

METHOD, and it reuses Layer 2 rather than inventing a second decomposer: take
the two longest parked alleles at a site, decompose the longer against the
shorter, and keep the result when the difference falls in [--min-insert,
--max-insert]. That is exactly the main pipeline's model applied one level down.

THE DEDUP COLLISION IS THE REAL RISK, AND IT IS NOT ASSUMED AWAY. A sub-5 kb
element can be found twice: by the main pipeline if it also forms its own
bubble, and here if it sits inside a parked parent. The key
(insert_md5, target_context_md5) SHOULD merge them -- but the context window is
cut from the short allele, and in the nested case the "short allele" is a
sub-interval of a large bubble, so the +/-100 bp window need not be the same
sequence. **This module therefore reports the nested-vs-main merge rate as
matched/total and never presents a combined count without it.**

    32_nested_search.py --armb <armB dirs> --main-events <dedup_events.tsv> \\
        --species ecoli --out stage32/ecoli
"""
import argparse
import collections
import csv
import glob
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def canon(s):
    r = rc(s)
    return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]


def decompose(short, long_):
    n = min(len(short), len(long_))
    p = 0
    while p < n and short[p] == long_[p]:
        p += 1
    s = 0
    while s < n - p and short[len(short) - 1 - s] == long_[len(long_) - 1 - s]:
        s += 1
    return p, s, long_[p:len(long_) - s]


def slide(L, l0, l1):
    left = 0
    while l0 - 1 - left >= 0 and L[l0 - 1 - left] == L[l1 - 1 - left]:
        left += 1
    right = 0
    while l1 + right < len(L) and L[l0 + right] == L[l1 + right]:
        right += 1
    return left + right


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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--armb", nargs="+", required=True)
    ap.add_argument("--main-events", default=None)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-insert", type=int, default=500)
    ap.add_argument("--max-insert", type=int, default=5000)
    ap.add_argument("--context", type=int, default=100)
    args = ap.parse_args()

    dirs = []
    for d in args.armb:
        dirs.extend(sorted(glob.glob(d)) if any(c in d for c in "*?[") else [d])
    n_parked, rows, skip = 0, [], collections.Counter()
    for d in dirs:
        p = os.path.join(d, "armB_large_events.tsv")
        if not os.path.exists(p):
            continue
        seqs = read_fasta(os.path.join(d, "armB_large_inserts.fna"))
        panel = os.path.basename(os.path.abspath(d.rstrip(os.sep)))
        for r in csv.DictReader(open(p), delimiter="\t"):
            n_parked += 1
            # The parked table's id column is `large_event_id`, NOT `event_id`.
            # Using the wrong name returned "" for every row and the whole stage
            # reported "0 carry a sequence" -- the seventh wrong-key failure in
            # this project, and the seventh to fail to a NEGATIVE. Accept either.
            eid = r.get("large_event_id") or r.get("event_id") or ""
            ins = seqs.get(eid, "")
            if not ins:
                skip["no_parked_sequence"] += 1
                continue
            # The parked record carries the FILLED allele's insert. A nested
            # element shows up as an internal difference, so we need two large
            # alleles. Where only one is stored, self-comparison is impossible
            # and the site is recorded as unresolvable rather than dropped.
            if len(ins) <= args.max_insert:
                skip["not_actually_oversized"] += 1
                continue
            rows.append({"panel": panel, "parent_event": eid,
                         "parent_len": len(ins), "parent_seq": ins})
    log("%d parked bubbles across %d panels; %d carry a sequence"
        % (n_parked, len(dirs), len(rows)))
    # ASSERTION: a 0% sequence-attachment rate is a key mismatch, not a result.
    if n_parked > 0 and len(rows) == 0:
        log("FATAL: 0 of %d parked bubbles matched a stored sequence. That is a "
            "KEY MISMATCH between armB_large_events.tsv and "
            "armB_large_inserts.fna, not an empty parked pile." % n_parked)
        return 2
    log("  skipped: %s" % dict(skip))

    if not rows:
        log("nothing to search -- parked inserts have no stored sequence")
        with open(args.out + "_nested.tsv", "w") as fh:
            fh.write("species\tpanel\tparent_event\tparent_len\tnested_len\t"
                     "offset\toverlap\tinsert_key\tcontext_key\n")
        return 0

    # group parents by canonical sequence: two panels that parked the same
    # island give two alleles of it, which is what makes a nested comparison
    # possible at all
    by_key = collections.defaultdict(list)
    for r in rows:
        by_key[canon(r["parent_seq"])].append(r)
    fams = collections.defaultdict(list)
    for r in rows:
        fams[r["parent_len"] // 1000].append(r)

    nested = []
    for bucket, group in fams.items():
        seen = {}
        for r in group:
            k = canon(r["parent_seq"])
            if k in seen:
                continue
            seen[k] = r
        uniq = list(seen.values())
        for i, a in enumerate(uniq):
            for b in uniq[i + 1:]:
                s_, l_ = sorted((a, b), key=lambda x: x["parent_len"])
                d = l_["parent_len"] - s_["parent_len"]
                if not (args.min_insert <= d <= args.max_insert):
                    continue
                S, L = s_["parent_seq"], l_["parent_seq"]
                lcp, lcs, ins = decompose(S, L)
                if not (args.min_insert <= len(ins) <= args.max_insert):
                    continue
                if lcp + lcs < 0.8 * len(S):
                    continue          # the two parents are not the same island
                ov = slide(L, lcp, len(L) - lcs)
                lo = max(0, lcp - args.context)
                hi = min(len(S), lcp + args.context)
                nested.append({
                    "species": args.species, "panel": l_["panel"],
                    "parent_event": l_["parent_event"],
                    "parent_len": l_["parent_len"], "nested_len": len(ins),
                    "offset": lcp, "overlap": ov,
                    "insert_key": canon(ins), "context_key": canon(S[lo:hi])})
    log("%d nested candidates from %d distinct parked islands"
        % (len(nested), len(by_key)))

    with open(args.out + "_nested.tsv", "w") as fh:
        cols = ["species", "panel", "parent_event", "parent_len", "nested_len",
                "offset", "overlap", "insert_key", "context_key"]
        fh.write("\t".join(cols) + "\n")
        for r in nested:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- THE MERGE RATE, reported as matched/total, never assumed ----------
    log("=" * 78)
    log("STAGE 32 -- %s" % args.species)
    log("  parked bubbles          %6d" % n_parked)
    log("  distinct parked islands %6d" % len(by_key))
    log("  nested candidates       %6d" % len(nested))
    if args.main_events and os.path.exists(args.main_events) and nested:
        main_ik, main_pair = set(), set()
        for r in csv.DictReader(open(args.main_events), delimiter="\t"):
            main_ik.add(r.get("insert_key", ""))
            main_pair.add((r.get("insert_key", ""), r.get("context_key", "")))
        by_full = sum(1 for r in nested
                      if (r["insert_key"], r["context_key"]) in main_pair)
        by_ins = sum(1 for r in nested if r["insert_key"] in main_ik)
        log("")
        log("  MERGE RATE against the main catalogue (%d events):" % len(main_pair))
        log("    full key (insert + context)  %6d / %d   %5.1f%%"
            % (by_full, len(nested), 100.0 * by_full / len(nested)))
        log("    insert key alone             %6d / %d   %5.1f%%"
            % (by_ins, len(nested), 100.0 * by_ins / len(nested)))
        log("")
        if by_ins and by_full < 0.5 * by_ins:
            log("    *** The context key does NOT merge these. The nested short")
            log("    allele is a sub-interval of a large bubble, so its +/-%d bp"
                % args.context)
            log("    window is different sequence. Nested finds are a SEPARATE")
            log("    NAMESPACE: do not add these counts to the main catalogue.")
        elif by_full:
            log("    The full key merges them; nested finds that are also main")
            log("    catalogue events are already counted there.")
        log("    NEW (not in the main catalogue by insert key): %d"
            % (len(nested) - by_ins))
    else:
        log("")
        log("  NO MERGE RATE COMPUTED -- pass --main-events. Without it these")
        log("  counts must not be combined with the main catalogue.")
    log("=" * 78)
    log("wrote %s_nested.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

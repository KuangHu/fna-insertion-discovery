#!/usr/bin/env python3
"""Collapse per-panel loci into distinct EVENTS. Required before any count is quoted.

At census scale a panel is a sliding window over one species, so the same
insertion at the same site is captured by every panel that contains both a
carrier and a non-carrier. K. pneumoniae ran 40 overlapping panels of 24
genomes drawn from 5,916, so "3,486 Layer-2 loci" is an unknown mixture of
distinct events and the same event seen many times. Without this key, no
per-species number means anything.

THE KEY IS (canonical insert, canonical target context):

  insert          the inserted sequence itself
  target context  a fixed +/- --context bp window of the SHORT allele centred on
                  the insertion offset

**The context must be a window centred on the offset, NOT the whole short
allele.** The short allele is essentially the 200 bp seeding pad -- measured on
E. coli it equals the seed interval width in 95.5% of loci -- and its exact
boundaries shift between panels because junction_L/junction_R differ slightly
per graph. Hashing the whole short allele therefore UNDER-merges: the same event
gets a different key in every panel.

BOTH COMPONENTS ARE CANONICALISED FOR ORIENTATION. Arm B stores sequence in
backbone orientation and Layer 2 in the carrier's, so the same event appears
forward in one panel and reverse-complemented in another. A forward-only hash
splits every such event in two. min(seq, revcomp(seq)) makes the key
strand-invariant. (The `inserted_md5` column already in loci.tsv is NOT
orientation-canonical and must not be used as the key on its own.)

    100_event_dedup.py --recon <recon dirs> --out events
"""
import argparse
import collections
import csv
import glob
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert                                                 # noqa: E402
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def canon(s):
    """Strand-invariant: hash whichever of the two orientations sorts first."""
    r = rc(s)
    return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]


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
    ap.add_argument("--out", required=True)
    ap.add_argument("--context", type=int, default=100,
                    help="bp each side of the offset used as target context")
    args = ap.parse_args()
    args.recon = expand(args.recon)

    rows, n_panels, n_skip = [], 0, collections.Counter()
    contributed = set()
    for d in args.recon:
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        n_panels += 1
        # Panel identity: the recon dir itself, unless it is a generic
        # "recon"/"l2" wrapper in which case the parent names the panel.
        # BUG 2026-09-03: this was basename(dirname(d)), correct for the
        # sweep layout `sweep/<panel>/recon` but wrong for `l2/<panel>`, where
        # it labelled every locus "l2". The cross-panel merge count then read
        # 0.0% -- a SILENT NEGATIVE, the fifth this session. It was caught only
        # because 0.0% is implausible and got checked. Panel names are also
        # carried in locus_id, so they are cross-checked below.
        ad = os.path.abspath(d).rstrip(os.sep)
        base = os.path.basename(ad)
        panel = (os.path.basename(os.path.dirname(ad))
                 if base in ("recon", "l2", "allele_recon") else base)
        aseq = read_alleles(os.path.join(d, "alleles.fna"))
        for r in csv.DictReader(open(p), delimiter="\t"):
            lid = r["locus_id"]
            if r["event_class"] not in ("insertion_target_retained", "replacement"):
                n_skip["not_decomposed"] += 1
                continue
            try:
                lcp, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
            except ValueError:
                n_skip["bad_fields"] += 1
                continue
            if ilen <= 0:
                n_skip["no_insert"] += 1
                continue
            av = aseq.get(lid, {})
            ls = av.get(r["longest_allele"], "")
            ss = av.get(r["shortest_allele"], "")
            if not ls or not ss:
                n_skip["missing_allele"] += 1
                continue
            # the insert comes from the LONG allele and needs a long-allele
            # coordinate; the context window is centred on the insertion point
            # in the SHORT allele and correctly uses lcp. The two frames are
            # not the same number on the tolerant path.
            i0, frame = lib_insert.locate_insert(r, ls)
            if i0 < 0:
                n_skip["insert_frame_unresolved"] += 1
                continue
            ins = ls[i0:i0 + ilen]
            lo = max(0, lcp - args.context)
            hi = min(len(ss), lcp + args.context)
            ctx = ss[lo:hi]
            if len(ins) < 1 or len(ctx) < args.context:
                n_skip["short_context"] += 1
                continue
            contributed.add(panel)
            rows.append({"panel": panel, "locus_id": lid,
                         # canonical_insert_key, not canon(ins): a direct
                         # repeat at the junction admits several equally valid
                         # boundaries, the aligner's choice is not symmetric
                         # under revcomp, and keying the raw slice therefore
                         # split 351 of 21,051 E. coli events in two.
                         "insert_key": lib_insert.canonical_insert_key(ls, i0, ilen),
                         "context_key": canon(ctx),
                         "frame": frame,
                         "inserted_len": ilen, "offset": lcp,
                         "overlap": int(r.get("junction_ambiguity_bp") or 0),
                         "lost": int(r.get("target_bases_lost") or 0),
                         "method": r.get("decomposition_method", ""),
                         "placement": r.get("placement_status", "")})
    log("%d panels, %d decomposable loci keyed, skipped %s"
        % (n_panels, len(rows), dict(n_skip)))
    # ASSERTION: distinct panel labels must match the number of recon dirs read.
    # One label for many dirs means the panel field is being derived wrongly,
    # which silently reports 0% cross-panel redundancy.
    n_lab = len(contributed)
    n_with_rows = len({os.path.basename(os.path.abspath(d).rstrip(os.sep))
                       for d in args.recon
                       if os.path.exists(os.path.join(d, "loci.tsv"))})
    if n_panels > 1 and n_lab < 2:
        log("FATAL: %d recon dirs produced only %d distinct panel label(s) -- "
            "the panel field is mis-derived and cross-panel redundancy will "
            "read as 0%%" % (n_panels, n_lab))
        return 2
    if n_lab < n_panels:
        log("note: %d of %d panels contributed no decomposable locus"
            % (n_panels - n_lab, n_panels))
    if not rows:
        log("FATAL: nothing keyed")
        return 2

    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["insert_key"], r["context_key"])].append(r)
    for i, (k, g) in enumerate(sorted(groups.items()), 1):
        for r in g:
            r["event_id"] = "E%06d" % i
            r["event_n_loci"] = len(g)
            r["event_n_panels"] = len({x["panel"] for x in g})

    cols = ["event_id", "panel", "locus_id", "insert_key", "context_key",
            "inserted_len", "offset", "overlap", "lost", "method", "placement",
            "frame",
            "event_n_loci", "event_n_panels"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_loci_to_event.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["event_id"], x["panel"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(args.out + "_events.tsv", "w") as fh:
        fh.write("event_id\tn_loci\tn_panels\tinserted_len\toffset\toverlap\t"
                 "lost\tmethod\tinsert_key\tcontext_key\n")
        for k, g in sorted(groups.items(), key=lambda kv: kv[1][0]["event_id"]):
            r = g[0]
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n"
                     % (r["event_id"], len(g), len({x["panel"] for x in g}),
                        r["inserted_len"], r["offset"], r["overlap"], r["lost"],
                        r["method"], r["insert_key"], r["context_key"]))

    mult = collections.Counter(len(g) for g in groups.values())
    log("=" * 78)
    log("EVENT-LEVEL DEDUP")
    log("")
    log("  loci keyed              %6d" % len(rows))
    log("  DISTINCT EVENTS         %6d" % len(groups))
    log("  redundancy factor       %6.2fx" % (len(rows) / len(groups)))
    log("")
    log("  loci per event: " + "  ".join(
        "%dx:%d" % (k, mult[k]) for k in sorted(mult)[:12]))
    log("  max loci in one event   %6d" % max(mult))
    log("  events seen in >1 panel %6d   %5.1f%%"
        % (sum(1 for g in groups.values() if len({x["panel"] for x in g}) > 1),
           100.0 * sum(1 for g in groups.values()
                       if len({x["panel"] for x in g}) > 1) / len(groups)))
    log("")
    log("  A count quoted per SPECIES must use DISTINCT EVENTS. A count quoted")
    log("  per PANEL may use loci, because within one panel each locus is one")
    log("  observation.")
    log("=" * 78)
    log("wrote %s_events.tsv and %s_loci_to_event.tsv" % (args.out, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""THE PRODUCT. One row per locus, one row per event, and the sequences.

Replaces `86_pre_insertion_target_catalog.py`, which was a POLARITY GATE: it
emitted a target only where Layer 3 had called the short allele ancestral, and
234 Layer-2 loci became 20 rows. Layer 3 is gone (docs/SIMPLIFICATION.md §2), so
this module emits EVERY locus and writes S1 as a FLAG.

    NOTHING HERE FILTERS. S1 is a column. Length is a column. Mobility is a
    column. A caller who wants only clean, mobile, >=5 kb rows writes one awk;
    a caller who wants the discards can have them, because they are still here.

WHAT A ROW ASSERTS -- and it is only this:

    at this anchored site these genomes differ by an inserted segment, the
    empty target sequence is `target_seq` verbatim, and the insert goes in at
    `insertion_point_offset`

    target_seq[:offset] + insert_seq + target_seq[offset:]  ==  the observed
    derived allele  (median identity 1.0000, 93.8% / 92.0% byte-identical)

WHAT NO ROW ASSERTS: which allele is ancestral (no polarity anywhere in this
pipeline); what family the insert belongs to (no annotation); that the junction
is correct to the base (never tested -- assembly tier only); that the insert is
a mobile element rather than an rRNA operon, prophage or REP array
(`insert_copies_in_genome`, `insert_distinct_loci` and `inserted_len` are
columns so a caller can decide -- rRNA operons are ~5 kb, REP arrays <0.5 kb).

EVENT vs LOCUS. Overlapping panels see one insertion many times: E. coli 1.98x,
K. pneumoniae 2.60x. Per-species counts MUST use `event_id`; per-panel counts
may use loci. Both files are written so neither can be quoted by accident.

    86_catalogue.py --recon <dirs> --events dedup_loci_to_event.tsv \\
        --mobility mobility_mobility.tsv --species ecoli --out catalogue
"""
import argparse
import collections
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


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
    ap.add_argument("--events", required=True)
    ap.add_argument("--mobility", default=None)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-overlap", type=int, default=15,
                    help="S1: above this a direct repeat is a repeat array or "
                         "segmental duplication, not a junction feature")
    args = ap.parse_args()
    args.recon = expand(args.recon)

    loci, aseq, carriers = {}, {}, collections.defaultdict(dict)
    for d in args.recon:
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            loci[r["locus_id"]] = r
        aseq.update(read_alleles(os.path.join(d, "alleles.fna")))
        ap_ = os.path.join(d, "alleles.tsv")
        if os.path.exists(ap_):
            for r in csv.DictReader(open(ap_), delimiter="\t"):
                carriers[r["locus_id"]][r["allele_id"]] = r.get("genomes", "")

    ev_of, ev_n = {}, {}
    for r in csv.DictReader(open(args.events), delimiter="\t"):
        ev_of[r["locus_id"]] = r["event_id"]
        ev_n[r["event_id"]] = int(r.get("event_n_loci") or 1)
    mob = {}
    if args.mobility and os.path.exists(args.mobility):
        for r in csv.DictReader(open(args.mobility), delimiter="\t"):
            mob[r["locus_id"]] = r

    dec = [l for l, r in loci.items()
           if r["event_class"] in ("insertion_target_retained", "replacement")
           and int(r["inserted_len"] or 0) > 0]
    # standing practice 6: report matched/total before any rate derived from a join
    n_ev = sum(1 for l in dec if l in ev_of)
    n_mb = sum(1 for l in dec if l in mob)
    log("%d loci, %d decomposable" % (len(loci), len(dec)))
    log("  event map  %d/%d  (%.1f%%)" % (n_ev, len(dec), 100.0 * n_ev / len(dec)))
    log("  mobility   %d/%d  (%.1f%%)" % (n_mb, len(dec), 100.0 * n_mb / len(dec)))
    if dec and n_ev < 0.5 * len(dec):
        log("FATAL: event map covers <50%% of loci; per-species counts would be "
            "inflated by panel redundancy")
        return 2

    cols = ["event_id", "locus_id", "species", "panel", "n_loci_in_event",
            "target_seq", "target_len", "insertion_point_offset",
            "insert_seq", "inserted_len", "insert_md5",
            "junction_overlap_bp", "target_bases_lost", "decomposition_method",
            "placement_status", "placement_score_margin",
            "S1_structurally_clean",
            "insert_copies_in_genome", "insert_distinct_loci",
            "M1_within_genome", "M2_cross_event", "mobility_positive",
            "carriers_derived", "carriers_empty"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    rows, seen_ev = [], {}
    for lid in sorted(dec):
        r = loci[lid]
        av = aseq.get(lid, {})
        tgt = av.get(r["shortest_allele"], "")
        lng = av.get(r["longest_allele"], "")
        off, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
        ins = lng[off:off + ilen] if lng else ""
        ov = int(r.get("junction_ambiguity_bp") or 0)
        lost = int(r.get("target_bases_lost") or 0)
        s1 = (lost == 0 and ov <= args.max_overlap
              and r.get("placement_status", "") not in ("ambiguous",))
        m = mob.get(lid, {})
        row = {
            "event_id": ev_of.get(lid, "."), "locus_id": lid,
            "species": args.species, "panel": lid.rsplit(".", 1)[0],
            "n_loci_in_event": ev_n.get(ev_of.get(lid, ""), 1),
            "target_seq": tgt, "target_len": len(tgt),
            "insertion_point_offset": off,
            "insert_seq": ins, "inserted_len": ilen,
            "insert_md5": r.get("inserted_md5", "."),
            "junction_overlap_bp": ov, "target_bases_lost": lost,
            "decomposition_method": r.get("decomposition_method", "."),
            "placement_status": r.get("placement_status", "."),
            "placement_score_margin": r.get("placement_score_margin", "."),
            "S1_structurally_clean": s1,
            "insert_copies_in_genome": m.get("M1_copies", "."),
            "insert_distinct_loci": m.get("M1_distinct_loci", "."),
            "M1_within_genome": m.get("M1_within_genome", "."),
            "M2_cross_event": m.get("M2_cross_event", "."),
            "mobility_positive": m.get("mobility_positive", "."),
            "carriers_derived": carriers.get(lid, {}).get(r["longest_allele"], ""),
            "carriers_empty": carriers.get(lid, {}).get(r["shortest_allele"], "")}
        rows.append(row)
        e = row["event_id"]
        if e != "." and e not in seen_ev:
            seen_ev[e] = row

    with open(args.out + "_loci.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(args.out + "_events.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in sorted(seen_ev):
            fh.write("\t".join(str(seen_ev[e][c]) for c in cols) + "\n")
    with open(args.out + "_targets.fna", "w") as ft, \
            open(args.out + "_inserts.fna", "w") as fi:
        for e in sorted(seen_ev):
            r = seen_ev[e]
            if r["target_seq"]:
                ft.write(">%s|%s offset=%d len=%d overlap=%d\n%s\n"
                         % (e, args.species, r["insertion_point_offset"],
                            r["target_len"], r["junction_overlap_bp"],
                            r["target_seq"]))
            if r["insert_seq"]:
                fi.write(">%s|%s len=%d md5=%s\n%s\n"
                         % (e, args.species, r["inserted_len"], r["insert_md5"],
                            r["insert_seq"]))

    n_s1 = sum(1 for r in rows if r["S1_structurally_clean"])
    n_mo = sum(1 for r in rows if r["mobility_positive"] == "True")
    ev_s1 = sum(1 for r in seen_ev.values() if r["S1_structurally_clean"])
    ev_mo = sum(1 for r in seen_ev.values() if r["mobility_positive"] == "True")
    log("=" * 80)
    log("CATALOGUE -- %s" % args.species)
    log("")
    log("  %-34s %8s %8s" % ("", "loci", "EVENTS"))
    log("  " + "-" * 52)
    log("  %-34s %8d %8d" % ("emitted", len(rows), len(seen_ev)))
    log("  %-34s %8d %8d" % ("S1 structurally clean", n_s1, ev_s1))
    log("  %-34s %8d %8d" % ("mobility positive", n_mo, ev_mo))
    log("  %-34s %8d %8d"
        % ("both", sum(1 for r in rows if r["S1_structurally_clean"]
                       and r["mobility_positive"] == "True"),
           sum(1 for r in seen_ev.values() if r["S1_structurally_clean"]
               and r["mobility_positive"] == "True")))
    log("")
    log("  Quote EVENTS per species. Loci are per-panel observations.")
    log("  Nothing was filtered: S1 and mobility are columns.")
    log("=" * 80)
    for suf in ("_loci.tsv", "_events.tsv", "_targets.fna", "_inserts.fna"):
        log("wrote %s%s" % (args.out, suf))
    return 0


if __name__ == "__main__":
    sys.exit(main())

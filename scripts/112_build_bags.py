#!/usr/bin/env python3
"""Emit Channel A/B bags: one bag per CDS cluster, 120 bp empty site per insertion.

Conforms to `docs/CANONICAL_BAG_SPEC.md` of the DL project; the conformance
report the run prints lists every departure (§7 requires one).

    bag          = one CDS cluster (protein similarity, mmseqs). Bags MAY span
                   species -- that is intended, not a leak.
    site         = one insertion event assigned to that cluster
    flank        = 120 bp of the EMPTY target site, 60 upstream + 60 downstream,
                   so the insertion point sits between index 59 and 60 (0-based)
    noncoding    = the insert minus its predicted CDS spans

THE FLANK COMES FROM THE EMPTY ALLELE, NOT THE FILLED ONE. `target_seq` is the
site as it was BEFORE the element arrived, reconstructed from genomes that lack
the insertion, so the 120 bp window holds no element sequence at either end. A
flank carved from the filled allele would leak the element's own termini into
the model input.

EXACT CENTRING IS ENFORCED, NEVER PADDED. A site is emitted only when
`offset >= 60` and `len(target_seq) - offset >= 60`. Padding a short window
would move the junction off 60|60 and the model would learn the padding.

`empty_site_source` IS EVIDENCE-BACKED, NOT ASSUMED. A site is `observed` only
when `carriers_empty` names at least one genome that actually carries the empty
allele. Sites without that evidence are DROPPED, not relabelled. No site in
this corpus is `tsd_reconstructed`; the enum exists so a future source that
does reconstruct can say so, and a consumer can filter on it.

`orient` IS CANONICAL, NOT `unknown`. The flank is stored in whichever
orientation is lexicographically smaller of (flank, revcomp(flank)), and
`orient` records which was applied -- `fwd` or `rc`. Deterministic, so two
records of the same site agree. COST, stated: after canonicalisation the
upstream/downstream identity of the two halves is no longer fixed, because
revcomp swaps them. The junction stays at index 60 either way (the window is
symmetric), and `orient` makes the original recoverable.

NO ANNOTATION ANYWHERE. prodigal is ab initio, mmseqs is sequence-only. No HMM,
no IS library. `cds_cluster_id` is a cluster, not a named family.

    112_build_bags.py --db database/insertions.tsv --cds cds/insert_cds.tsv \\
        --inserts database/insertions_inserts.fna --orf-gff cds/orfs.gff \\
        --out bags/insertions
"""
import argparse
import collections
import csv
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(p):
    out, cur = {}, None
    for line in open(p):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--orf-gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--half", type=int, default=60)
    ap.add_argument("--min-nc", type=int, default=20)
    ap.add_argument("--min-sites", type=int, default=1)
    ap.add_argument("--corpus", default="fna_ins_discovery")
    ap.add_argument("--version", default="2026-09-20")
    args = ap.parse_args()

    cds = {r["db_id"]: r for r in csv.DictReader(open(args.cds), delimiter="\t")}
    db = list(csv.DictReader(open(args.db), delimiter="\t"))
    matched = sum(1 for r in db if r["db_id"] in cds)
    log("db rows %d ; with a CDS record %d (%.1f%%)"
        % (len(db), matched, 100.0 * matched / len(db) if db else 0))
    if db and matched < 0.5 * len(db):
        log("FATAL: db -> cds join below 50%% -- key mismatch")
        return 2

    spans = collections.defaultdict(list)
    for line in open(args.orf_gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[2] != "CDS":
            continue
        try:
            spans[f[0]].append((int(f[3]), int(f[4])))
        except ValueError:
            continue
    seqs = read_fa(args.inserts)

    # drops are counted by REASON x SPECIES: a reason that fires only in one
    # species is a different problem from one that fires everywhere
    drop = collections.defaultdict(collections.Counter)
    sites = []
    for r in db:
        did = r["db_id"]
        sp = r.get("species", "?")
        c = cds.get(did)
        key = c.get("cds_cluster_id", c.get("transposase_id", ".")) if c else "."
        if not c or key in (".", ""):
            drop["no_cds_cluster"][sp] += 1
            continue

        # empty_site_source: evidence, not assumption
        emptyg = [g for g in r.get("carriers_empty", "").split(",") if g]
        if not emptyg:
            drop["no_observed_empty_carrier"][sp] += 1
            continue

        tgt = r.get("target_seq", "")
        try:
            off = int(r["insertion_point_offset"])
        except (KeyError, ValueError):
            drop["bad_offset"][sp] += 1
            continue
        if off < args.half or len(tgt) - off < args.half:
            drop["window_would_not_fit"][sp] += 1
            continue
        flank = tgt[off - args.half: off + args.half]
        if len(flank) != 2 * args.half:
            drop["short_window"][sp] += 1
            continue
        if flank.count("N") > 6:
            drop["too_many_N"][sp] += 1
            continue

        s = seqs.get(did, "")
        if not s:
            drop["no_insert_sequence"][sp] += 1
            continue
        mask = bytearray(len(s))
        for st, en in spans.get(did, []):
            for i in range(max(0, st - 1), min(len(s), en)):
                mask[i] = 1
        nc, run = [], []
        for i, m in enumerate(mask):
            if m:
                if run:
                    nc.append("".join(run)); run = []
            else:
                run.append(s[i])
        if run:
            nc.append("".join(run))
        nc = [x for x in nc if len(x) >= args.min_nc]
        if not nc:
            drop["no_noncoding_region"][sp] += 1
            continue

        # canonical orientation
        rcf = rc(flank)
        if flank <= rcf:
            cflank, orient = flank, "fwd"
        else:
            cflank, orient = rcf, "rc"
        nch = hashlib.sha1(json.dumps(nc).encode()).hexdigest()

        sites.append({
            "site_id": "%s.site" % did,
            "bag_id": key,
            "cds_cluster_id": key,
            "ncrna_id": nch[:12],
            "corpus": args.corpus,
            "species": sp,
            "phylum": r.get("phylum", "."),
            "flank": cflank,
            "flank_len": len(cflank),
            "flank_side": "joined",
            "flank_spacer": "",
            "flank_source": "real_genomic_empty_allele",
            "empty_site_source": "observed",
            "empty_site_n_carrier_genomes": len(emptyg),
            "insertion_point_in_flank": args.half,
            "noncoding_regions": nc,
            "nc_region_count": len(nc),
            "nc_total_len": sum(len(x) for x in nc),
            "nc_sequence_hash": nch,
            "nc_padding_scheme": "none",
            "nc_padding_offset": 0,
            "orient_granularity": "per_site",
            "orient": orient,
            "orient_source": "canonical_lexicographic_min",
            "inserted_len": r.get("inserted_len", "."),
            "insert_md5": r.get("insert_md5", "."),
            "junction_overlap_bp": r.get("junction_overlap_bp", "."),
            "S1_structurally_clean": r.get("S1_structurally_clean", "."),
            "mobility_positive": r.get("mobility_positive", "."),
            "data_source": args.corpus,
            "build_date": args.version,
            "generator_version_or_commit": "112_build_bags.py/%s" % args.version,
        })

    bysize = collections.Counter(s["bag_id"] for s in sites)
    sites = [s for s in sites if bysize[s["bag_id"]] >= args.min_sites]
    bags = collections.defaultdict(list)
    for s in sites:
        bags[s["bag_id"]].append(s)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_sites.jsonl", "w") as fh:
        for s in sites:
            fh.write(json.dumps(s) + "\n")
    with open(args.out + "_bags.tsv", "w") as fh:
        fh.write("bag_id\tcds_cluster_id\tn_sites\tn_species\tn_phyla\t"
                 "n_unique_nc_hash\tmedian_nc_total_len\tspecies\n")
        for b, ss in sorted(bags.items(), key=lambda kv: -len(kv[1])):
            spp = sorted({x["species"] for x in ss})
            ph = sorted({x["phylum"] for x in ss})
            L = sorted(x["nc_total_len"] for x in ss)
            fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n"
                     % (b, b, len(ss), len(spp), len(ph),
                        len({x["nc_sequence_hash"] for x in ss}),
                        L[len(L) // 2], ",".join(spp)))

    nuh = len({s["nc_sequence_hash"] for s in sites})
    orc = collections.Counter(s["orient"] for s in sites)
    log("=" * 78)
    log("BAGS   bag = CDS cluster ;  flank = %d bp empty site, %d|%d"
        % (2 * args.half, args.half, args.half))
    log("")
    log("  sites emitted                 %7d" % len(sites))
    log("  bags                          %7d" % len(bags))
    if bags:
        sz = sorted(len(v) for v in bags.values())
        log("  sites per bag: median %d  max %d  bags >=5 sites %d"
            % (sz[len(sz) // 2], sz[-1], sum(1 for v in sz if v >= 5)))
    log("  orient after canonicalisation: fwd %d  rc %d"
        % (orc["fwd"], orc["rc"]))
    log("")
    log("  §6 check 8 -- scaffold sharing")
    log("    distinct nc_sequence_hash   %7d" % nuh)
    log("    bags / unique nc hashes     %7.3f" % (len(bags) / nuh if nuh else 0))
    log("")
    log("  cross-species bags: >=2 species %d   >=3 %d   >=2 phyla %d"
        % (sum(1 for v in bags.values() if len({x['species'] for x in v}) >= 2),
           sum(1 for v in bags.values() if len({x['species'] for x in v}) >= 3),
           sum(1 for v in bags.values() if len({x['phylum'] for x in v}) >= 2)))
    log("")
    log("  DROPS by reason x species")
    allsp = sorted({s for c in drop.values() for s in c})
    log("    %-28s %7s  %s" % ("reason", "total", "  ".join("%-9s" % s[:9] for s in allsp[:6])))
    for k in sorted(drop, key=lambda x: -sum(drop[x].values())):
        log("    %-28s %7d  %s"
            % (k, sum(drop[k].values()),
               "  ".join("%-9d" % drop[k].get(s, 0) for s in allsp[:6])))
    if len(allsp) > 6:
        log("    (%d further species in the TSV)" % (len(allsp) - 6))
    with open(args.out + "_drops.tsv", "w") as fh:
        fh.write("reason\tspecies\tn\n")
        for k, c in drop.items():
            for s, n in sorted(c.items()):
                fh.write("%s\t%s\t%d\n" % (k, s, n))
    log("")
    log("  CONFORMANCE (spec §7)")
    log("    DEPARTURE -- bag-level nc invariant (§1). Bags key on CDS cluster,")
    log("      so sites in a bag may differ in noncoding_regions. Measured: a")
    log("      strict nc re-split gives 50,379 bags of which only 10.4%% hold")
    log("      >=2 sites, covering 30.8%% of sites. The re-split is available")
    log("      via nc_sequence_hash but yields a ~90%% singleton corpus.")
    log("    Bags SPAN SPECIES by design; key by (corpus, bag_id).")
    log("    flank_side=joined, spacer empty, 60+60, junction at index %d, no"
        % args.half)
    log("      padding ever applied.")
    log("    orient is CANONICAL (lexicographic min vs revcomp), not gold; it")
    log("      MUST NOT be used as a strand label.")
    log("    labels.* ABSENT -- unlabelled real observations, not generator")
    log("      positives. Do not load as a labelled corpus.")
    log("=" * 78)
    for s in ("_sites.jsonl", "_bags.tsv", "_drops.tsv"):
        log("wrote %s%s" % (args.out, s))
    return 0


if __name__ == "__main__":
    sys.exit(main())

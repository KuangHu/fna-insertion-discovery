#!/usr/bin/env python3
"""BLOCKING CHECK: is "the consumer can re-split on nc_sequence_hash" a real escape?

`112_build_bags.py` keys bags on CDS cluster, which departs from
CANONICAL_BAG_SPEC §1 (all sites in a bag must share IDENTICAL
`noncoding_regions`). The departure was going to be justified by saying a
consumer can recover the strict form by re-splitting on `nc_sequence_hash`.

That escape only exists if the strict split leaves usable bags. So measure it
BEFORE any bag is written:

    after re-splitting every CDS cluster on exact nc_sequence_hash,
      - what FRACTION OF BAGS have >= 2 sites?
      - what FRACTION OF SITES do those bags hold?

If bags of size >= 2 are a single-digit percentage, the escape does not exist
and the bag key has to be reconsidered rather than documented around.

Nothing is emitted here. One number, then a decision.

    113_nc_resplit_check.py --db database/insertions.tsv \\
        --cds cds/insert_cds.tsv --inserts database/insertions_inserts.fna \\
        --orf-gff cds/orfs.gff
"""
import argparse
import collections
import csv
import hashlib
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


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
    ap.add_argument("--min-nc", type=int, default=20)
    args = ap.parse_args()

    cds = {r["db_id"]: r for r in csv.DictReader(open(args.cds), delimiter="\t")}
    db = list(csv.DictReader(open(args.db), delimiter="\t"))
    matched = sum(1 for r in db if r["db_id"] in cds)
    log("db rows %d ; with a CDS record %d (%.1f%%)"
        % (len(db), matched, 100.0 * matched / len(db) if db else 0))
    if db and matched < 0.5 * len(db):
        log("FATAL: db -> cds join below 50%%; key mismatch")
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

    # cluster -> list of nc hashes
    byclu = collections.defaultdict(list)
    n_used = 0
    for r in db:
        did = r["db_id"]
        c = cds.get(did)
        if not c:
            continue
        clu = c.get("transposase_id", c.get("cds_cluster_id", "."))
        if clu in (".", ""):
            continue
        s = seqs.get(did, "")
        if not s:
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
            continue
        byclu[clu].append(hashlib.sha1(json.dumps(nc).encode()).hexdigest())
        n_used += 1

    log("%d sites carry a CDS cluster and >=1 non-coding region" % n_used)
    if n_used == 0:
        log("FATAL: nothing to measure")
        return 2

    # the bag as proposed: one per CDS cluster
    loose_bags = len(byclu)
    loose_ge2 = sum(1 for v in byclu.values() if len(v) >= 2)
    loose_sites_ge2 = sum(len(v) for v in byclu.values() if len(v) >= 2)

    # the strict spec form: re-split each cluster on exact nc hash
    strict = collections.Counter()
    for clu, hs in byclu.items():
        for h, n in collections.Counter(hs).items():
            strict[(clu, h)] = n
    s_bags = len(strict)
    s_ge2 = sum(1 for v in strict.values() if v >= 2)
    s_sites_ge2 = sum(v for v in strict.values() if v >= 2)

    log("=" * 74)
    log("NC RE-SPLIT CHECK -- does the strict spec form leave usable bags?")
    log("")
    log("  %-34s %10s %10s" % ("", "CDS-cluster", "strict nc"))
    log("  " + "-" * 56)
    log("  %-34s %10d %10d" % ("bags", loose_bags, s_bags))
    log("  %-34s %10d %10d" % ("bags with >=2 sites", loose_ge2, s_ge2))
    log("  %-34s %9.1f%% %9.1f%%"
        % ("  as a share of bags",
           100.0 * loose_ge2 / loose_bags if loose_bags else 0,
           100.0 * s_ge2 / s_bags if s_bags else 0))
    log("  %-34s %10d %10d" % ("sites in those bags", loose_sites_ge2, s_sites_ge2))
    log("  %-34s %9.1f%% %9.1f%%"
        % ("  as a share of sites",
           100.0 * loose_sites_ge2 / n_used if n_used else 0,
           100.0 * s_sites_ge2 / n_used if n_used else 0))
    log("")
    sz = sorted(strict.values())
    log("  strict bag sizes: median %d  max %d" % (sz[len(sz) // 2], sz[-1]))
    log("")
    pb = 100.0 * s_ge2 / s_bags if s_bags else 0
    ps = 100.0 * s_sites_ge2 / n_used if n_used else 0
    if pb < 10 or ps < 10:
        log("  *** THE RE-SPLIT ESCAPE DOES NOT EXIST. Only %.1f%% of strict bags"
            % pb)
        log("  hold >=2 sites, covering %.1f%% of sites. Saying 'a consumer can" % ps)
        log("  re-split on nc_sequence_hash' would hand over a corpus of")
        log("  singletons. The bag key must be decided, not documented around.")
    else:
        log("  The strict form retains %.1f%% of bags at >=2 sites, covering" % pb)
        log("  %.1f%% of sites -- the re-split escape is real." % ps)
    log("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())

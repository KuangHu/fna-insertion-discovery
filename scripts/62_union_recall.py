#!/usr/bin/env python3
"""Benchmark 3 -- Recall(A u B). The number that actually matters.

Arm A is a HIGH-COPY engine: measured P(detect | ncopy) is 10% at 2 copies,
37% at 3, and 89-95% from 5 up. Its low-copy blind spot is a design property,
not a bug, and the fix is not to loosen Arm A -- it is Arm B, which keys on
between-genome empty/filled alleles and does not care how many copies exist.

So the question is never Recall(A). It is:

    Recall(A u B)   on known complete IS,
    while the ArmA+/ISEScan- pool stays clean.

MATCHING. Arm A is matched to ISEScan positionally (stage 60 already did that,
and its per_is table is the input here). Arm B CANNOT be matched the same way:
bubble coordinates live in backbone space, and the backbone is often not the
genome carrying the IS. Matching Arm B positionally in the wrong frame is the
same mistake that made junction spans come out element-sized. So Arm B is
matched by SEQUENCE: the ISEScan IS instance is aligned against the Arm B
inserted sequences from the cluster that contains its genome. That asks the
right question -- did Arm B independently discover this element as an
empty/filled polymorphism?

SCOPE. Only genomes that stage 10 actually clustered can be scored: Arm B is
structurally unable to see an unclustered genome, so including one would
charge Arm B for a genome it was never given. Restricting the denominator is
reported, not silent.

    62_union_recall.py --per-is bench60_per_is.tsv \\
        --clusters armB_bench/stage10/clusters.tsv \\
        --inserts-glob 'armB_bench/stage30/*/armB_inserts.fna' \\
        --genome-dir isescan_run/fna --isescan-glob 'isescan_run/tsv/*.tsv' \\
        --out bench62
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def read_clusters(path):
    """assembly path -> cluster_id, keyed by sample name."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = os.path.basename(r["assembly"])
            for suf in (".gz", ".fna", ".fa", ".fasta"):
                if s.endswith(suf):
                    s = s[: -len(suf)]
            out[s] = r["cluster_id"]
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-is", required=True, help="bench60_per_is.tsv")
    ap.add_argument("--clusters", required=True, help="stage 10 clusters.tsv")
    ap.add_argument("--inserts-glob", required=True,
                    help="armB_inserts.fna per cluster; cluster id taken from "
                         "the parent directory name")
    ap.add_argument("--genome-dir", required=True,
                    help="directory of <sample>.fna used by the ISEScan run")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="fraction of the IS instance covered by the Arm B insert")
    ap.add_argument("--exclude-family", default="ISNCY",
                    help="comma-separated ISEScan labels that are not real "
                         "families; ISNCY is a catch-all whose members do not "
                         "align to each other and must not be scored")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()

    excl = {x for x in args.exclude_family.split(",") if x}
    clu = read_clusters(args.clusters)

    # ---- Arm B inserts, grouped by the cluster they came from --------------
    ins_by_cluster = {}
    for p in glob.glob(args.inserts_glob):
        cid = os.path.basename(os.path.dirname(p))
        ins_by_cluster[cid] = p
    log("Arm B insert sets: %s" % ", ".join(sorted(ins_by_cluster)) or "none")

    # ---- the ISEScan calls to score ----------------------------------------
    rows = []
    with open(args.per_is) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["is_type"] != "c" or r["class"] != "multi":
                continue
            if r["is_family"] in excl:
                continue
            rows.append(r)
    scored = [r for r in rows if clu.get(r["sample"])]
    log("%d multi-copy complete IS calls; %d in clustered genomes (%d genomes "
        "clustered of %d seen)"
        % (len(rows), len(scored), len({r['sample'] for r in scored}),
           len({r['sample'] for r in rows})))

    # ---- extract each IS instance, align to its cluster's Arm B inserts ----
    tmp = tempfile.mkdtemp(prefix="union_")
    per_cluster = defaultdict(list)
    for i, r in enumerate(scored):
        per_cluster[clu[r["sample"]]].append((i, r))

    armb_hit = set()
    for cid, items in sorted(per_cluster.items()):
        ins = ins_by_cluster.get(cid)
        if not ins or not os.path.getsize(ins):
            log("  %s: no Arm B inserts -- %d calls unscorable by Arm B"
                % (cid, len(items)))
            continue
        q = os.path.join(tmp, "%s_is.fna" % cid)
        with open(q, "w") as out:
            for i, r in items:
                fa = os.path.join(args.genome_dir, r["sample"] + ".fna")
                reg = "%s:%s-%s" % (r["contig"], r["is_start"], r["is_end"])
                try:
                    seq = subprocess.run(
                        [args.samtools, "faidx", fa, reg],
                        capture_output=True, text=True, check=True).stdout
                except subprocess.CalledProcessError:
                    continue
                body = "".join(seq.split("\n")[1:])
                if body:
                    out.write(">%d\n%s\n" % (i, body))
        paf = subprocess.run(
            [args.minimap2, "-x", "asm10", "-c", "-N", "200", "-p", "0",
             "--secondary=yes", ins, q],
            capture_output=True, text=True).stdout
        for line in paf.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            qlen, nmatch, blen = int(f[1]), int(f[9]), int(f[10])
            ident = 100.0 * nmatch / max(1, blen)
            cov = nmatch / max(1, qlen)
            if ident >= args.min_ident and cov >= args.min_cov:
                armb_hit.add(int(f[0]))
        log("  %s: %d IS calls, %d matched an Arm B insert"
            % (cid, len(items), sum(1 for i, _ in items if i in armb_hit)))

    # ---- tabulate -----------------------------------------------------------
    def unit(r):
        return (r["sample"], r["is_family"])

    per_copy = {"A": 0, "B": 0, "U": 0, "n": 0}
    fam_A, fam_B, fam_all = set(), set(), set()
    detail = []
    for i, r in enumerate(scored):
        a = r["found_by_armA"] == "1"
        b = i in armb_hit
        per_copy["n"] += 1
        per_copy["A"] += a
        per_copy["B"] += b
        per_copy["U"] += (a or b)
        u = unit(r)
        fam_all.add(u)
        if a:
            fam_A.add(u)
        if b:
            fam_B.add(u)
        detail.append([r["sample"], r["contig"], r["is_start"], r["is_end"],
                       r["is_family"], r["is_copy_number"], int(a), int(b),
                       int(a or b)])

    with open(args.out + "_per_is.tsv", "w") as fh:
        fh.write("sample\tcontig\tis_start\tis_end\tis_family\tis_copy_number\t"
                 "found_armA\tfound_armB\tfound_union\n")
        for d in detail:
            fh.write("\t".join(map(str, d)) + "\n")

    def pct(a, b):
        return 100.0 * a / b if b else float("nan")

    fam_U = fam_A | fam_B
    with open(args.out + "_summary.tsv", "w") as fh:
        fh.write("level\tarm\tdetected\ttotal\trecall_pct\n")
        for lvl, d in (("per_copy", per_copy),):
            for k in "ABU":
                fh.write("%s\t%s\t%d\t%d\t%.1f\n"
                         % (lvl, k, d[k], d["n"], pct(d[k], d["n"])))
        for k, s in (("A", fam_A), ("B", fam_B), ("U", fam_U)):
            fh.write("family\t%s\t%d\t%d\t%.1f\n"
                     % (k, len(s), len(fam_all), pct(len(s), len(fam_all))))

    log("=" * 72)
    log("RECALL on multi-copy COMPLETE IS, clustered genomes only (%s excluded)"
        % ",".join(sorted(excl)))
    log("")
    log("  FAMILY-LEVEL (the honest number -- read this one)")
    log("    Arm A %5.1f%%   Arm B %5.1f%%   A u B %5.1f%%   (n=%d units)"
        % (pct(len(fam_A), len(fam_all)), pct(len(fam_B), len(fam_all)),
           pct(len(fam_U), len(fam_all)), len(fam_all)))
    log("")
    log("  per-copy (Arm B column is OPTIMISTIC -- see note)")
    log("    Arm A %5.1f%%   Arm B %5.1f%%   A u B %5.1f%%   (n=%d calls)"
        % (pct(per_copy["A"], per_copy["n"]), pct(per_copy["B"], per_copy["n"]),
           pct(per_copy["U"], per_copy["n"]), per_copy["n"]))
    log("")
    log("  NOTE: Arm A is matched POSITIONALLY, Arm B by SEQUENCE, so the two")
    log("  per-copy columns do not mean the same thing. One Arm B insert")
    log("  matches every instance of that element in the cluster, so Arm B is")
    log("  credited with copies it never localised: here %d calls matched from"
        % per_copy["B"])
    log("  only %d distinct (genome, family) units. Compare arms at FAMILY"
        % len(fam_B))
    log("  level; per-copy Arm B is an upper bound, not a detection rate.")
    log("")
    rescued = fam_B - fam_A
    log("  Arm B RESCUES %d family units Arm A missed:" % len(rescued))
    byfam = defaultdict(int)
    for _, f in rescued:
        byfam[f] += 1
    for f, n in sorted(byfam.items(), key=lambda kv: -kv[1]):
        log("      %-14s %d" % (f, n))
    still = fam_all - fam_U
    log("  still missed by BOTH: %d" % len(still))
    byfam = defaultdict(int)
    for _, f in still:
        byfam[f] += 1
    for f, n in sorted(byfam.items(), key=lambda kv: -kv[1])[:8]:
        log("      %-14s %d" % (f, n))
    log("=" * 72)
    log("wrote %s_{per_is,summary}.tsv" % args.out)


if __name__ == "__main__":
    main()

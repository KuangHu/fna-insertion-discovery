#!/usr/bin/env python3
"""Stage 10 -- group assemblies into ANI clusters with skani.

Empty-vs-filled comparison is only meaningful between genomes whose loci are
actually homologous. Throwing E. coli and Pseudomonas into one graph produces
"insertions" that are really different genomic architecture.

Two guards, both required:
  --min-ani  sequence identity   (default 99.0)
  --min-af   ALIGNED FRACTION    (default 85.0)   <-- the one that matters

Aligned fraction is the real protection: two genomes can be 99% identical over
the 30% of themselves that aligns and share no architecture at all.

Outputs:
  clusters.tsv    assembly -> cluster_id, plus per-genome QC
  clusters/<id>.manifest   one FNA path per line, backbone first
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import UnionFind, log, run                        # noqa: E402


def assembly_qc(path):
    """n_contigs, total_len, N50, n_frac -- read once, cheaply."""
    lens, ns, tot = [], 0, 0
    cur = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                s = line.strip()
                cur += len(s)
                ns += s.upper().count("N")
    if cur:
        lens.append(cur)
    tot = sum(lens)
    lens.sort(reverse=True)
    acc, n50 = 0, 0
    for L in lens:
        acc += L
        if acc >= tot / 2:
            n50 = L
            break
    return len(lens), tot, n50, (ns / tot if tot else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fna-list", required=True,
                    help="file with one assembly FASTA path per line")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-ani", type=float, default=99.0)
    ap.add_argument("--min-af", type=float, default=85.0,
                    help="min aligned fraction %% (the real homology guard)")
    ap.add_argument("--min-cluster-size", type=int, default=4,
                    help="a cluster needs enough genomes to see both alleles")
    ap.add_argument("--max-cluster-size", type=int, default=200,
                    help="graph engines degrade past a few hundred genomes")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--min-n50", type=int, default=20000,
                    help="drop assemblies below this N50 (junctions unresolvable)")
    ap.add_argument("--max-contigs", type=int, default=1000)
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "clusters"), exist_ok=True)
    paths = [l.strip() for l in open(args.fna_list) if l.strip()]
    log("%d assemblies in" % len(paths))

    # ---- QC gate first: a shattered assembly cannot support a junction call
    qc, kept = {}, []
    for p in paths:
        nc, tot, n50, nf = assembly_qc(p)
        qc[p] = (nc, tot, n50, nf)
        if n50 >= args.min_n50 and nc <= args.max_contigs:
            kept.append(p)
    log("%d pass QC (N50>=%d, contigs<=%d); %d dropped"
        % (len(kept), args.min_n50, args.max_contigs, len(paths) - len(kept)))

    keep_list = os.path.join(args.outdir, "qc_pass.txt")
    with open(keep_list, "w") as fh:
        fh.write("\n".join(kept) + "\n")

    # ---- skani all-vs-all
    dist = os.path.join(args.outdir, "skani_dist.tsv")
    if not os.path.exists(dist):
        log("skani triangle ...")
        cmd = [args.skani, "dist", "--ql", keep_list, "--rl", keep_list,
               "-t", str(args.threads), "-o", dist,
               "--min-af", str(min(args.min_af, 50.0))]
        run(cmd)

    # ---- single-linkage on (ANI, AF)
    uf = UnionFind()
    for p in kept:
        uf.add(p)
    n_edge = 0
    with open(dist) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            r, q = f[col["Ref_file"]], f[col["Query_file"]]
            if r == q:
                continue
            ani = float(f[col["ANI"]])
            af = min(float(f[col["Align_fraction_ref"]]),
                     float(f[col["Align_fraction_query"]]))
            if ani >= args.min_ani and af >= args.min_af:
                uf.union(r, q)
                n_edge += 1
    log("%d edges at ANI>=%.2f AF>=%.1f" % (n_edge, args.min_ani, args.min_af))

    groups = [g for g in uf.groups() if len(g) >= args.min_cluster_size]
    groups.sort(key=len, reverse=True)
    log("%d clusters with >=%d members" % (len(groups), args.min_cluster_size))

    out = open(os.path.join(args.outdir, "clusters.tsv"), "w")
    out.write("cluster_id\tassembly\tis_backbone\tn_contigs\ttotal_len\tn50\tn_frac\n")
    for i, g in enumerate(groups):
        cid = "cl%04d" % i
        # backbone = most contiguous assembly; graph is built on it first
        g = sorted(g, key=lambda p: (-qc[p][2], qc[p][0]))
        if len(g) > args.max_cluster_size:
            log("cluster %s: %d members, keeping the %d most contiguous"
                % (cid, len(g), args.max_cluster_size))
            g = g[:args.max_cluster_size]
        with open(os.path.join(args.outdir, "clusters", cid + ".manifest"), "w") as mf:
            mf.write("\n".join(g) + "\n")
        for j, p in enumerate(g):
            nc, tot, n50, nf = qc[p]
            out.write("%s\t%s\t%d\t%d\t%d\t%d\t%.5f\n"
                      % (cid, p, int(j == 0), nc, tot, n50, nf))
    out.close()
    log("wrote clusters.tsv + %d manifests" % len(groups))


if __name__ == "__main__":
    main()

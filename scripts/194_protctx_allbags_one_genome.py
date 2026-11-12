#!/usr/bin/env python3
"""All 6,139 bags vs ONE genome in a single tblastn, then cut the contexts.

WHY ONE SEARCH FOR ALL BAGS. tblastn translates the subject in six frames once
and reuses it for every query, so cost is far below one-search-per-bag.
Measured on a 4.7Mb genome at 8 threads: 1 query 5.4s, 100 queries 10.7s,
6139 queries 514s -- i.e. ~15,000 CPU-hours for the corpus against ~74,000 for
6,139 separate passes.

WHY ONE REPRESENTATIVE PER BAG. Using every member's dominant ORF would be
105,048 queries, ~320,000 CPU-hours. The query id IS the bag, so no join is
needed afterwards.

OUTPUT IS ONE FILE PER WORKER, NOT PER BAG. 6,139 bags x 16,309 genomes would
create up to 100M files on a Lustre filesystem that would not survive it.
Each record is one line (bag, header, sequence) and the demultiplex into
per-bag files happens once, later, in a streaming pass.
"""
import argparse
import collections
import os
import subprocess
import sys

# IDENTITY IS THE LENGTH-WEIGHTED MEAN OVER THE MERGED HSPs, NOT THE MAXIMUM.
# max(pident) reads a locus whose HSPs are 40aa at 95% and 200aa at 35% as
# "95% identical". The weighted mean calls it 45%, which is what it is. Both
# are written to the hits table; the GATE uses the weighted mean.
GAP_MAX = 600
SPAN_SLACK = 1.6
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, k, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k:
                    out[k] = "".join(buf)
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k:
        out[k] = "".join(buf)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--out-recs", required=True, help="append records here")
    ap.add_argument("--out-hits", required=True, help="append hit rows here")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--min-ident", type=float, default=80.0)
    ap.add_argument("--min-qcov", type=float, default=80.0)
    ap.add_argument("--evalue", default="1e-10")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--tblastn", default="tblastn")
    a = ap.parse_args()

    gname = os.path.basename(a.genome)
    for ext in (".fna", ".fa", ".fasta"):
        if gname.endswith(ext):
            gname = gname[:-len(ext)]
            break

    p = subprocess.run(
        [a.tblastn, "-query", a.queries, "-subject", a.genome, "-outfmt",
         "6 qseqid sseqid pident length evalue bitscore sstart send qstart "
         "qend qlen", "-evalue", a.evalue, "-seg", "yes",
         "-num_threads", a.threads],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    if p.returncode != 0:
        print("FAIL %s: %s" % (gname, p.stderr.strip()[:300]), file=sys.stderr)
        sys.exit(1)

    byloc = collections.defaultdict(list)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        bag = f[0].split("|")[0]
        ss, se = int(f[6]), int(f[7])
        strand = "plus" if ss <= se else "minus"
        lo, hi = (ss, se) if ss <= se else (se, ss)
        byloc[(f[0], f[1], strand)].append(
            dict(bag=bag, lo=lo, hi=hi, pid=float(f[2]), bits=float(f[5]),
                 qs=int(f[8]), qe=int(f[9]), qlen=int(f[10])))
    if not byloc:
        return

    loci = []
    for (qid, sid, strand), hs in byloc.items():
        hs.sort(key=lambda x: x["lo"])
        max_span = int(hs[0]["qlen"] * 3 * SPAN_SLACK) + 60
        cur = [hs[0]]
        for h in hs[1:]:
            if (h["lo"] - cur[-1]["hi"] <= GAP_MAX
                    and h["hi"] - cur[0]["lo"] + 1 <= max_span):
                cur.append(h)
            else:
                loci.append((qid, sid, strand, cur))
                cur = [h]
        loci.append((qid, sid, strand, cur))

    # dedup WITHIN a bag only: two different bags hitting one locus is a real
    # observation about that locus, not a duplicate to collapse
    best = {}
    for qid, sid, strand, hs in loci:
        lo = min(h["lo"] for h in hs)
        hi = max(h["hi"] for h in hs)
        qlen = hs[0]["qlen"]
        acov = sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)
        rec = dict(q=qid, bag=hs[0]["bag"], s=sid, strand=strand, lo=lo, hi=hi,
                   bits=sum(h["bits"] for h in hs),
                   ident=sum(h["pid"] * (abs(h["qe"] - h["qs"]) + 1)
                             for h in hs)
                   / max(1, sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)),
                   ident_best=max(h["pid"] for h in hs),
                   qcov=100.0 * min(acov, qlen) / qlen, nhsp=len(hs))
        key = (rec["bag"], sid, strand, lo // 500, hi // 500)
        if key not in best or rec["bits"] > best[key]["bits"]:
            best[key] = rec

    seqs = read_fasta(a.genome)
    fr = open(a.out_recs, "a")
    fh = open(a.out_hits, "a")
    for r in best.values():
        ctg = seqs.get(r["s"])
        if ctg is None:
            continue
        clen = len(ctg)
        left = min(a.flank, r["lo"] - 1)
        right = min(a.flank, clen - r["hi"])
        full = "yes" if (left == a.flank and right == a.flank) else "no"
        ok = r["ident"] >= a.min_ident and r["qcov"] >= a.min_qcov
        fh.write("\t".join(str(v) for v in [
            r["bag"], r["q"], gname, r["s"], clen, r["lo"], r["hi"],
            r["strand"], "%.1f" % r["ident"], "%.1f" % r["ident_best"],
            "%.1f" % r["qcov"], r["nhsp"],
            "%.0f" % r["bits"], left, right, full,
            "yes" if ok else "no"]) + "\n")
        if not ok:
            continue
        sub = ctg[r["lo"] - 1 - left: r["hi"] + right]
        if r["strand"] == "minus":
            sub = revcomp(sub)
            up, dn = right, left
        else:
            up, dn = left, right
        hdr = ("%s|%s|%s|%d-%d|%s|ident=%.1f|qcov=%.1f|up=%d|dn=%d|full=%s"
               " query=%s" % (r["bag"], gname, r["s"], r["lo"], r["hi"],
                              r["strand"], r["ident"], r["qcov"], up, dn,
                              full, r["q"]))
        # one line per record: bag, header, sequence -- demultiplexed later
        fr.write("%s\t%s\t%s\n" % (r["bag"], hdr, sub))
    fr.close()
    fh.close()


if __name__ == "__main__":
    main()

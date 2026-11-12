#!/usr/bin/env python3
"""tblastn every bag against ONE genome and cut [flank][CDS][flank] in place.

WHY PER GENOME. The first design concatenated all 16,309 assemblies into one
62GB BLAST database. That is a large fixed cost which only pays off when the
database is searched many times; here it is searched three times with 19 short
proteins. Both the concatenation and makeblastdb are single threaded, so a
64-core node ran at load 1.04 for 90 minutes. Measured: tblastn -subject
against one 4.7MB genome takes ~2s per bag, so the whole corpus is ~27
CPU-hours -- about 25 minutes across 64 cores, with no intermediate at all.

The HSP-merging and flank rules are unchanged from 190_, which the 60-genome
smoke test validated at 100% centred placement on all three bags.
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
GAP_MAX = 600        # max subject gap between consecutive HSPs of one gene
SPAN_SLACK = 1.6     # a locus may reach SPAN_SLACK x the query coding length
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
    ap.add_argument("--queries", nargs="+", required=True,
                    help="bag faa files; the bag name is the basename")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--min-ident", type=float, default=80.0)
    ap.add_argument("--min-qcov", type=float, default=80.0)
    ap.add_argument("--evalue", default="1e-10")
    ap.add_argument("--tblastn", default="tblastn")
    a = ap.parse_args()

    gname = os.path.basename(a.genome)
    for ext in (".fna", ".fa", ".fasta"):
        if gname.endswith(ext):
            gname = gname[:-len(ext)]
            break

    seqs = None    # loaded lazily: most genomes have no hit for some bags
    for q in a.queries:
        bag = os.path.basename(q)
        if bag.endswith(".faa"):
            bag = bag[:-4]
        p = subprocess.run(
            [a.tblastn, "-query", q, "-subject", a.genome, "-outfmt",
             "6 qseqid sseqid pident length evalue bitscore sstart send "
             "qstart qend qlen", "-evalue", a.evalue, "-seg", "yes"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
        if p.returncode != 0:
            print("FAIL %s %s: %s" % (gname, bag, p.stderr.strip()[:200]),
                  file=sys.stderr)
            sys.exit(1)
        rows = [l.split("\t") for l in p.stdout.splitlines() if l.strip()]
        if not rows:
            continue

        byloc = collections.defaultdict(list)
        for f in rows:
            if len(f) < 11:
                continue
            ss, se = int(f[6]), int(f[7])
            strand = "plus" if ss <= se else "minus"
            lo, hi = (ss, se) if ss <= se else (se, ss)
            byloc[(f[0], f[1], strand)].append(
                dict(lo=lo, hi=hi, pid=float(f[2]), bits=float(f[5]),
                     qs=int(f[8]), qe=int(f[9]), qlen=int(f[10])))

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

        # one record per genomic copy: several member proteins of a bag hit the
        # same copy, and emitting each would duplicate the locus per member
        best = {}
        for qid, sid, strand, hs in loci:
            lo = min(h["lo"] for h in hs)
            hi = max(h["hi"] for h in hs)
            bits = sum(h["bits"] for h in hs)
            qlen = hs[0]["qlen"]
            acov = sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)
            rec = dict(q=qid, s=sid, strand=strand, lo=lo, hi=hi, bits=bits,
                       ident=sum(h["pid"] * (abs(h["qe"] - h["qs"]) + 1)
                             for h in hs)
                   / max(1, sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)),
                   ident_best=max(h["pid"] for h in hs),
                       qcov=100.0 * min(acov, qlen) / qlen, nhsp=len(hs))
            key = (sid, strand, lo // 500, hi // 500)
            if key not in best or bits > best[key]["bits"]:
                best[key] = rec
        if not best:
            continue

        if seqs is None:
            seqs = read_fasta(a.genome)

        os.makedirs(a.outdir, exist_ok=True)
        fh_fa = open("%s/%s.%s.fna" % (a.outdir, bag, gname), "w")
        fh_tsv = open("%s/%s.%s.tsv" % (a.outdir, bag, gname), "w")
        for r in sorted(best.values(), key=lambda x: -x["bits"]):
            ctg = seqs.get(r["s"])
            if ctg is None:
                continue
            clen = len(ctg)
            left = min(a.flank, r["lo"] - 1)
            right = min(a.flank, clen - r["hi"])
            full = "yes" if (left == a.flank and right == a.flank) else "no"
            ok = r["ident"] >= a.min_ident and r["qcov"] >= a.min_qcov
            fh_tsv.write("\t".join(str(v) for v in [
                bag, r["q"], gname, r["s"], clen, r["lo"], r["hi"],
                r["strand"], "%.1f" % r["ident"], "%.1f" % r["ident_best"],
                "%.1f" % r["qcov"], r["nhsp"], "%.0f" % r["bits"], left, right,
                full,
                "yes" if ok else "no"]) + "\n")
            if not ok:
                continue
            sub = ctg[r["lo"] - 1 - left: r["hi"] + right]
            if r["strand"] == "minus":
                sub = revcomp(sub)
                up, dn = right, left    # after revcomp the far side leads
            else:
                up, dn = left, right
            fh_fa.write(">%s|%s|%s|%d-%d|%s|ident=%.1f|qcov=%.1f|up=%d|dn=%d"
                        "|full=%s query=%s\n%s\n"
                        % (bag, gname, r["s"], r["lo"], r["hi"], r["strand"],
                           r["ident"], r["qcov"], up, dn, full, r["q"], sub))
        fh_fa.close()
        fh_tsv.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Find each query protein across a genome and emit [flank][CDS][flank] JSON.

Built for the NP2 candidate set, whose FASTA headers carry metadata tblastn
discards: qseqid is only the first token, so bags=, roles=, kinds= and Pfam=
are parsed from the header here and carried into every record. Without that a
hit could not be attributed to its bag.

JSONL, NOT ONE JSON ARRAY. A single array has to be read whole before any of
it is usable; at this record count that means holding the entire file in
memory. One object per line streams, greps and splits.

The HSP merge bounds and flank accounting are the ones validated by the
60-genome smoke test: HSPs of one gene merge, separate copies do not, and a
flank that runs off the contig end is reported rather than silently short.
"""
import argparse
import collections
import json
import os
import subprocess
import sys

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


def query_meta(path):
    """qseqid -> the metadata in its header, which tblastn drops."""
    meta = {}
    for line in open(path):
        if not line.startswith(">"):
            continue
        parts = line[1:].split()
        qid = parts[0]
        d = {"query": qid}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                d[k] = v
        meta[qid] = d
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--min-ident", type=float, default=80.0)
    ap.add_argument("--min-qcov", type=float, default=80.0)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--tblastn", default="tblastn")
    a = ap.parse_args()

    gname = os.path.basename(a.genome)
    for ext in (".fna", ".fa", ".fasta"):
        if gname.endswith(ext):
            gname = gname[:-len(ext)]
            break
    qmeta = query_meta(a.queries)

    p = subprocess.run(
        [a.tblastn, "-query", a.queries, "-subject", a.genome, "-outfmt",
         "6 qseqid sseqid pident length evalue bitscore sstart send qstart "
         "qend qlen", "-evalue", a.evalue, "-seg", "yes", "-num_threads", "1"],
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
        ss, se = int(f[6]), int(f[7])
        strand = "plus" if ss <= se else "minus"
        lo, hi = (ss, se) if ss <= se else (se, ss)
        byloc[(f[0], f[1], strand)].append(
            dict(lo=lo, hi=hi, pid=float(f[2]), bits=float(f[5]),
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

    best = {}
    for qid, sid, strand, hs in loci:
        lo = min(h["lo"] for h in hs)
        hi = max(h["hi"] for h in hs)
        qlen = hs[0]["qlen"]
        acov = sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)
        wsum = sum(abs(h["qe"] - h["qs"]) + 1 for h in hs)
        rec = dict(q=qid, s=sid, strand=strand, lo=lo, hi=hi,
                   bits=sum(h["bits"] for h in hs),
                   ident=sum(h["pid"] * (abs(h["qe"] - h["qs"]) + 1)
                             for h in hs) / max(1, wsum),
                   ident_best=max(h["pid"] for h in hs),
                   qcov=100.0 * min(acov, qlen) / qlen, nhsp=len(hs))
        key = (qid, sid, strand, lo // 500, hi // 500)
        if key not in best or rec["bits"] > best[key]["bits"]:
            best[key] = rec

    seqs = read_fasta(a.genome)
    with open(a.out, "a") as out:
        for r in best.values():
            ctg = seqs.get(r["s"])
            if ctg is None:
                continue
            clen = len(ctg)
            left = min(a.flank, r["lo"] - 1)
            right = min(a.flank, clen - r["hi"])
            if r["ident"] < a.min_ident or r["qcov"] < a.min_qcov:
                continue
            sub = ctg[r["lo"] - 1 - left: r["hi"] + right]
            if r["strand"] == "minus":
                sub = revcomp(sub)
                up, dn = right, left
            else:
                up, dn = left, right
            m = qmeta.get(r["q"], {})
            out.write(json.dumps({
                "query": r["q"],
                "bag": m.get("bags", "."),
                "role": m.get("roles", "."),
                "kind": m.get("kinds", "."),
                "query_aa": int(m.get("aa", 0) or 0),
                "pfam": m.get("Pfam", "none"),
                "genome": gname,
                "species": a.species,
                "contig": r["s"],
                "contig_len": clen,
                "cds_start": r["lo"],
                "cds_end": r["hi"],
                "cds_len": r["hi"] - r["lo"] + 1,
                "strand": r["strand"],
                "pct_ident": round(r["ident"], 1),
                "pct_ident_best": round(r["ident_best"], 1),
                "q_cov": round(r["qcov"], 1),
                "n_hsp": r["nhsp"],
                "bitscore": round(r["bits"]),
                "flank_up": up,
                "flank_dn": dn,
                "flank_full": (up == a.flank and dn == a.flank),
                "seq_len": len(sub),
                "seq": sub,
            }, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()

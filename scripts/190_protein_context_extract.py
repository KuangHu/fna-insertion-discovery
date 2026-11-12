#!/usr/bin/env python3
"""From tblastn hits of a bag's insert protein, cut [2000bp][CDS][2000bp].

INPUT  tblastn -outfmt 6 against a BLAST db built from the whole genome corpus.
OUTPUT one FASTA per bag, each record a protein-coding hit with its genomic
       context, oriented so the protein reads FORWARD in the emitted sequence.

WHY HSPs ARE MERGED. tblastn splits one gene into several HSPs whenever a
frameshift, an intron-like gap or a low-complexity stretch breaks the
alignment. Treating each HSP as a separate locus multiplies one real copy into
three or four, which is the same inflation that made the cross-species counts
unusable. HSPs from one (query, subject, strand) whose subject ranges are
within MERGE_GAP are one locus, and the locus span is min(start) to max(end).

WHY THE FLANK IS RECORDED, NOT ASSUMED. A hit near a contig end cannot have
2000bp on both sides. The record is still emitted -- dropping it would bias
the set toward closed genomes -- but left_bp/right_bp say what was actually
obtained, so a downstream filter can require the full window.

THRESHOLDS ARE NOT BAKED IN. Every hit that passed BLAST's own e-value is
written to <bag>_hits.tsv. The FASTA applies --min-ident/--min-qcov, and both
are recorded in the header, so a different cut needs no re-run of the search.
"""
import argparse
import collections
import os
import subprocess
import sys

# HSPs of ONE gene merge; copies of a repeated element must not. The bound is
# the gene itself: a locus may not exceed the query's own coding length by more
# than SPAN_SLACK, and consecutive HSPs may not be further apart than GAP_MAX.
#
# MEASURED. A flat 3000bp gap merged ~10 tandem copies of CDS06776 into single
# loci -- median span 6369bp for a 648bp protein -- and 42 of 49 records then
# failed to translate back to the query. The two bags whose copies are isolated
# (spans 558bp and 690bp) were unaffected, so the error was invisible in
# aggregate and only showed up per-record.
GAP_MAX = 600        # max subject gap between consecutive HSPs of one gene
SPAN_SLACK = 1.6     # a locus may reach SPAN_SLACK x the query coding length


def extract_batch(db, jobs, blastdbcmd, tmp):
    """Pull every region in ONE blastdbcmd call via -entry_batch.

    One subprocess per locus costs ~50ms of fork and db-open overhead, which a
    high-copy transposase turns into hours: 100,000 loci is 1.4 hours of pure
    process churn. -entry_batch opens the database once. The batch file format
    is `seqid range strand`, one per line, and results come back IN INPUT
    ORDER, which is what lets them be zipped back onto the job list.
    """
    if not jobs:
        return {}
    with open(tmp, "w") as fh:
        for (sid, s, e, strand) in jobs:
            fh.write("%s %d-%d %s\n" % (sid, s, e, strand))
    cmd = [blastdbcmd, "-db", db, "-entry_batch", tmp, "-outfmt", "%s"]
    p = subprocess.run(cmd, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    if p.returncode != 0:
        print("WARN: blastdbcmd batch failed: %s" % p.stderr.strip()[:300],
              file=sys.stderr)
        return {}
    seqs = [l.strip() for l in p.stdout.splitlines() if l.strip()]
    if len(seqs) != len(jobs):
        # a refused entry yields no line, so positional zip would shift every
        # later sequence onto the wrong locus -- silently corrupting the file
        print("WARN: asked %d regions, got %d; falling back to per-entry"
              % (len(jobs), len(seqs)), file=sys.stderr)
        out = {}
        for j in jobs:
            q = subprocess.run([blastdbcmd, "-db", db, "-entry", j[0],
                                "-range", "%d-%d" % (j[1], j[2]),
                                "-strand", j[3], "-outfmt", "%s"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
            if q.returncode == 0:
                out[j] = "".join(q.stdout.split())
        return out
    return dict(zip(jobs, seqs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="tblastn outfmt-6 table")
    ap.add_argument("--db", required=True, help="BLAST nucleotide db prefix")
    ap.add_argument("--map", required=True, help="seqid<TAB>genome<TAB>contig")
    ap.add_argument("--bag", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--min-ident", type=float, default=30.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--blastdbcmd", default="blastdbcmd")
    a = ap.parse_args()

    meta = {}
    for line in open(a.map):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            meta[f[0]] = (f[1], f[2], int(f[3]))   # genome, contig, length

    # ---- group HSPs into loci
    cols = "qseqid sseqid pident length evalue bitscore sstart send qstart qend qlen".split()
    byloc = collections.defaultdict(list)
    n_hsp = 0
    for line in open(a.hits):
        f = line.rstrip("\n").split("\t")
        if len(f) < len(cols):
            continue
        r = dict(zip(cols, f))
        n_hsp += 1
        ss, se = int(r["sstart"]), int(r["send"])
        strand = "plus" if ss <= se else "minus"
        lo, hi = (ss, se) if ss <= se else (se, ss)
        r["_lo"], r["_hi"], r["_strand"] = lo, hi, strand
        byloc[(r["qseqid"], r["sseqid"], strand)].append(r)

    loci = []
    for (q, s, strand), hsps in byloc.items():
        hsps.sort(key=lambda x: x["_lo"])
        max_span = int(int(hsps[0]["qlen"]) * 3 * SPAN_SLACK) + 60
        cur = [hsps[0]]
        for h in hsps[1:]:
            gap = h["_lo"] - cur[-1]["_hi"]
            span = h["_hi"] - cur[0]["_lo"] + 1
            if gap <= GAP_MAX and span <= max_span:
                cur.append(h)
            else:
                loci.append((q, s, strand, cur))
                cur = [h]
        loci.append((q, s, strand, cur))

    # ---- one record per (subject, strand, span): keep the best query only.
    # Several member proteins of one bag hit the same genomic copy; emitting
    # each would duplicate the locus once per member.
    best = {}
    for q, s, strand, hsps in loci:
        lo = min(h["_lo"] for h in hsps)
        hi = max(h["_hi"] for h in hsps)
        bits = sum(float(h["bitscore"]) for h in hsps)
        ident = max(float(h["pident"]) for h in hsps)
        qlen = int(hsps[0]["qlen"])
        acov = sum(abs(int(h["qend"]) - int(h["qstart"])) + 1 for h in hsps)
        qcov = 100.0 * min(acov, qlen) / qlen
        key = (s, strand, lo // 500, hi // 500)
        rec = dict(q=q, s=s, strand=strand, lo=lo, hi=hi, bits=bits,
                   ident=ident, qcov=qcov, nhsp=len(hsps))
        if key not in best or bits > best[key]["bits"]:
            best[key] = rec

    # plan every extraction first so they go out as one batch
    plan, meta_by_rec = [], {}
    for rec in sorted(best.values(), key=lambda r: -r["bits"]):
        g, c, clen = meta.get(rec["s"], ("?", rec["s"], 0))
        left = min(a.flank, rec["lo"] - 1)
        right = min(a.flank, clen - rec["hi"]) if clen else a.flank
        ok = rec["ident"] >= a.min_ident and rec["qcov"] >= a.min_qcov
        meta_by_rec[id(rec)] = (g, c, clen, left, right, ok)
        if ok:
            plan.append((rec["s"], rec["lo"] - left, rec["hi"] + right,
                         rec["strand"]))
    got = extract_batch(a.db, plan, a.blastdbcmd, a.out + ".batch.txt")
    print("  extraction: requested %d regions, got %d" % (len(plan), len(got)),
          file=sys.stderr)

    ht = open(a.out + "_hits.tsv", "w")
    ht.write("\t".join(["bag", "query", "genome", "contig", "contig_len",
                        "cds_start", "cds_end", "strand", "pct_ident",
                        "q_cov", "n_hsp", "bitscore", "left_bp", "right_bp",
                        "flank_complete", "passed"]) + "\n")
    fa = open(a.out + ".fna", "w")
    n_pass = n_emit = n_trunc = n_fail = 0
    for rec in sorted(best.values(), key=lambda r: -r["bits"]):
        g, c, clen, left, right, ok = meta_by_rec[id(rec)]
        full = "yes" if (left == a.flank and right == a.flank) else "no"
        ht.write("\t".join(str(v) for v in [
            a.bag, rec["q"], g, c, clen, rec["lo"], rec["hi"], rec["strand"],
            "%.1f" % rec["ident"], "%.1f" % rec["qcov"], rec["nhsp"],
            "%.0f" % rec["bits"], left, right, full,
            "yes" if ok else "no"]) + "\n")
        if not ok:
            continue
        n_pass += 1
        seq = got.get((rec["s"], rec["lo"] - left, rec["hi"] + right,
                       rec["strand"]))
        if not seq:
            n_fail += 1
            if n_fail <= 5:
                print("WARN: extraction failed for %s %d-%d" %
                      (rec["s"], rec["lo"], rec["hi"]), file=sys.stderr)
            elif n_fail == 6:
                print("WARN: ... further extraction failures suppressed",
                      file=sys.stderr)
            continue
        n_emit += 1
        if full == "no":
            n_trunc += 1
        # on the minus strand blastdbcmd already revcomps, so the CDS reads
        # forward and the UPSTREAM flank is the one printed first either way
        up, dn = (left, right) if rec["strand"] == "plus" else (right, left)
        fa.write(">%s|%s|%s|%d-%d|%s|ident=%.1f|qcov=%.1f|up=%d|dn=%d|full=%s"
                 " query=%s\n%s\n" % (
                     a.bag, g, c, rec["lo"], rec["hi"], rec["strand"],
                     rec["ident"], rec["qcov"], up, dn, full, rec["q"], seq))
    ht.close()
    fa.close()
    print("%s: %d HSPs -> %d loci -> %d passed (ident>=%.0f qcov>=%.0f)"
          " -> %d written, %d with a truncated flank, %d extraction failures"
          % (a.bag, n_hsp, len(best), n_pass, a.min_ident, a.min_qcov,
             n_emit, n_trunc, n_fail))
    if n_pass and not n_emit:
        sys.exit("FATAL: %d loci passed but NOTHING was extracted. The BLAST "
                 "db almost certainly lacks -parse_seqids, so blastdbcmd "
                 "cannot resolve any id. An empty FASTA with exit 0 is worse "
                 "than a crash -- refusing to report success." % n_pass)
    if n_pass and n_fail > n_pass * 0.05:
        sys.exit("FATAL: %d of %d extractions failed (>5%%); output would be "
                 "silently incomplete." % (n_fail, n_pass))


if __name__ == "__main__":
    main()

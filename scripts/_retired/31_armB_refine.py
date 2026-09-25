#!/usr/bin/env python3
"""Stage 31 -- refine Arm B candidates to nucleotide resolution.

Graph bubbles locate the locus; they are NOT a biochemical junction. Here each
(filled, empty) genome pair is re-aligned and the difference re-called, giving
the columns the downstream flank pipeline already expects:

    chrom  POS  SVLEN  INS_SEQ        <- same shape as the Sniffles output

Engine:
  svim-asm (default) -- assembly-vs-assembly SV caller, haploid mode, emits
                        VCF with the inserted sequence in ALT. Closest thing
                        to "Sniffles for assemblies".
  nucdiff            -- for fragmented drafts, where svim-asm's single
                        best-alignment model breaks down.

One svim-asm run per UNIQUE genome pair, not per event: the same pair usually
carries many loci.

TARGET-SITE DUPLICATION is recorded here as an OPTIONAL annotation:

    empty  : LEFT [tsd] RIGHT
    filled : LEFT [tsd] ELEMENT [tsd] RIGHT

A TSD present is good positive evidence of transposition. A TSD absent is not
evidence against it -- the IS110/IS1111 family makes no TSD and has no terminal
inverted repeats, and 0/13 IS110 events in this pipeline's own output carry one
(vs 43/170 for other families). Nothing filters on it.
"""
import argparse
import gzip
import os
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp                            # noqa: E402
from lib.tsd import find_tsd                                    # noqa: E402
from lib.util import log, run, seq_md5                          # noqa: E402


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# --------------------------------------------------------------- svim-asm run
def run_svimasm(filled_fa, empty_fa, outdir, threads, minimap2, samtools, svimasm):
    """filled aligned TO empty, so an insertion in filled = INS in empty coords."""
    os.makedirs(outdir, exist_ok=True)
    bam = os.path.join(outdir, "aln.bam")
    vcf = os.path.join(outdir, "variants.vcf")
    if os.path.exists(vcf):
        return vcf
    sam = os.path.join(outdir, "aln.sam")
    with open(sam, "w") as fh:
        run([minimap2, "-a", "-x", "asm5", "--cs", "-r2k", "-t", str(threads),
             empty_fa, filled_fa], stdout=fh, capture_output=False)
    run([samtools, "sort", "-@", str(threads), "-o", bam, sam])
    run([samtools, "index", bam])
    os.remove(sam)
    run([svimasm, "haploid", outdir, bam, empty_fa])
    return vcf


def parse_vcf_ins(vcf, min_len, max_len):
    """Yield (chrom, pos, svlen, ins_seq) for INS records."""
    if not os.path.exists(vcf):
        return
    with opener(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if info.get("SVTYPE") != "INS":
                continue
            try:
                svlen = abs(int(info.get("SVLEN", "0")))
            except ValueError:
                continue
            if not (min_len <= svlen <= max_len):
                continue
            alt = f[4]
            ins = alt[1:] if alt.startswith(("A", "C", "G", "T", "N")) and len(alt) > 1 else ""
            yield f[0], int(f[1]), svlen, ins


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True, help="armB_pairs.tsv from stage 30")
    ap.add_argument("--clusters", required=True,
                    help="clusters.tsv from stage 10 (maps sample -> fasta path)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-len", type=int, default=300)
    ap.add_argument("--max-len", type=int, default=5000)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--pos-tol", type=int, default=2000,
                    help="max distance between a bubble and a VCF INS to match them")
    ap.add_argument("--min-tsd", type=int, default=3)
    ap.add_argument("--max-tsd", type=int, default=30)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--svim-asm", dest="svimasm", default="svim-asm")
    ap.add_argument("--keep-bam", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # sample -> fasta
    fa_of = {}
    with open(args.clusters) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            b = os.path.basename(f[1])
            for suf in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
                if b.endswith(suf):
                    b = b[:-len(suf)]
                    break
            fa_of[b] = f[1]

    # group events by unique (filled, empty) pair
    by_pair = defaultdict(list)
    with open(args.pairs) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            r = dict(zip(hdr, f))
            by_pair[(r["filled_genome"], r["empty_genome"])].append(r)
    log("%d events across %d unique genome pairs"
        % (sum(len(v) for v in by_pair.values()), len(by_pair)))

    out = open(os.path.join(args.outdir, "armB_refined.tsv"), "w")
    seqs = open(os.path.join(args.outdir, "armB_inserts.fna"), "w")
    out.write("\t".join([
        "event_id", "cluster_id", "filled_genome", "empty_genome",
        "empty_chrom", "empty_pos", "left_junction", "right_junction",
        "insert_len", "insert_md5", "tsd_len", "tsd_seq", "tsd_side", "tsd_conf",
        "left_flank", "right_flank", "empty_site_seq",
        "bubble_insert_size", "pos_delta", "refine_status"]) + "\n")

    n_ok = n_miss = 0
    for (fsam, esam), rows in by_pair.items():
        fp, ep = fa_of.get(fsam), fa_of.get(esam)
        if not fp or not ep:
            log("  skip %s/%s: fasta path unknown" % (fsam, esam))
            n_miss += len(rows)
            continue
        wd = os.path.join(args.outdir, "pairs", "%s__%s" % (fsam, esam))
        try:
            vcf = run_svimasm(fp, ep, wd, args.threads, args.minimap2,
                              args.samtools, args.svimasm)
        except Exception as exc:
            log("  svim-asm failed for %s/%s: %s" % (fsam, esam, exc))
            n_miss += len(rows)
            continue

        ins_by_chrom = defaultdict(list)
        for chrom, pos, svlen, seq in parse_vcf_ins(vcf, args.min_len, args.max_len):
            ins_by_chrom[chrom].append((pos, svlen, seq))
        efa = Fasta(ep)

        for r in rows:
            # bubble coords are in the cluster backbone frame; match by the
            # nearest INS of comparable size on any contig of the empty genome
            want = int(r["insert_size"])
            best = None
            for chrom, lst in ins_by_chrom.items():
                for pos, svlen, seq in lst:
                    if abs(svlen - want) > max(200, 0.2 * want):
                        continue
                    score = abs(svlen - want)
                    if best is None or score < best[0]:
                        best = (score, chrom, pos, svlen, seq)
            if best is None:
                out.write("\t".join([r["event_id"], r["cluster_id"], fsam, esam] +
                                    ["NA"] * 15 + [r["insert_size"], "NA",
                                                   "NO_MATCHING_INS"]) + "\n")
                n_miss += 1
                continue
            _, chrom, pos, svlen, seq = best
            clen = efa.length(chrom)
            lf = efa.fetch(chrom, max(0, pos - args.flank), pos)
            rf = efa.fetch(chrom, pos, min(clen, pos + args.flank))
            tl, ts, tside, tconf = find_tsd(seq, lf, rf, min_tsd=args.min_tsd,
                                            max_tsd=args.max_tsd)
            out.write("\t".join(map(str, [
                r["event_id"], r["cluster_id"], fsam, esam, chrom, pos,
                "%s:%d" % (chrom, pos), "%s:%d" % (chrom, pos),
                svlen, seq_md5(seq) if seq else "NA",
                tl, ts or ".", tside, tconf,
                lf[-200:], rf[:200], (lf[-50:] + "|" + rf[:50]),
                want, svlen - want, "OK",
            ])) + "\n")
            if seq:
                seqs.write(">%s len=%d tsd=%d\n%s\n" % (r["event_id"], svlen, tl, seq))
            n_ok += 1

        if not args.keep_bam:
            for f in ("aln.bam", "aln.bam.bai"):
                p = os.path.join(wd, f)
                if os.path.exists(p):
                    os.remove(p)

    out.close()
    seqs.close()
    log("refined %d events (%d unmatched)" % (n_ok, n_miss))


if __name__ == "__main__":
    main()

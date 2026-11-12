#!/usr/bin/env python3
"""Reconstruct the empty allele around every full copy, so copies can be
collapsed into insertion SITES.

THE THREE COUNTS, KEPT APART.
  1. copies                 -- how many places carry a homolog of the element.
                               That is what the blastn sweep gives, and it is
                               NOT a number of insertions.
  2. distinct insertion sites -- copies sitting at the same orthologous locus in
                               different genomes are ONE site.  This script
                               supplies the key that lets them be collapsed.
  3. independent transposition events -- how many times the element actually
                               moved in history.  Presence/absence cannot give
                               this; it needs polarity on a tree.  Nothing here
                               estimates it, and nothing downstream should
                               report a number as if it did.

THE KEY.  For a copy at [s0,s1) on a contig, the empty allele is the
concatenation of the left and right flanks with the element deleted.  Two
copies at the same locus in two genomes reconstruct nearly the same empty
allele; two copies at different loci do not.  The same string is also the probe
for stage 3: blast it against genomes that lack the element and a contiguous
full-length match IS the empty site, which is what turns a copy into an
insertion polymorphism.

Canonical orientation: min(seq, revcomp(seq)) so a copy found on either strand
gives the same key, the same rule the insertion dedup already uses.

Copies whose flank would run off the contig end are kept but marked truncated,
because dropping them would silently bias against copies near contig breaks --
exactly where assembly-based IS detection is weakest.
"""
import argparse, collections, csv, os, sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if name: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", required=True)
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--kinds", default="full_copy")
    a = ap.parse_args()
    kinds = set(a.kinds.split(","))

    path = {}
    for r in csv.DictReader(open(a.corpus), delimiter="\t"):
        path[r["assembly"]] = r["fna_path"]

    by_asm = collections.defaultdict(list)
    n = 0
    for r in csv.DictReader(open(a.copies), delimiter="\t"):
        if r["copy_kind"] not in kinds: continue
        by_asm[r["assembly"]].append(r); n += 1
    sys.stderr.write("copies to flank: %d over %d assemblies\n" % (n, len(by_asm)))

    out = open(a.out + "_flanks.tsv", "w")
    fa = open(a.out + "_empty_ctx.fna", "w")
    out.write("copy_id\tbag_id\tassembly\tspecies\tcontig\ts0\ts1\tstrand\t"
              "left_trunc\tright_trunc\tctx_len\tempty_ctx_canon_sha\n")
    import hashlib
    done = trunc = missing = 0
    for asm, rows in sorted(by_asm.items()):
        p = path.get(asm)
        if not p or not os.path.exists(p):
            missing += len(rows); continue
        seqs = read_fasta(p)
        for r in rows:
            ctg = seqs.get(r["contig"])
            if ctg is None: missing += 1; continue
            s0, s1 = int(r["s0"]), int(r["s1"])
            l0 = max(0, s0 - a.flank); r1 = min(len(ctg), s1 + a.flank)
            L, R = ctg[l0:s0], ctg[s1:r1]
            lt = int(len(L) < a.flank); rt = int(len(R) < a.flank)
            trunc += (lt or rt)
            ctx = (L + R).upper()
            if not ctx: continue
            canon = min(ctx, rc(ctx))
            cid = "%s|%s|%s|%d|%d" % (r["bag_id"], asm, r["contig"], s0, s1)
            out.write("\t".join([cid, r["bag_id"], asm, r["species"], r["contig"],
                                 str(s0), str(s1), r["strand"], str(lt), str(rt),
                                 str(len(ctx)),
                                 hashlib.sha256(canon.encode()).hexdigest()[:16]]) + "\n")
            fa.write(">%s\n%s\n" % (cid, canon))
            done += 1
    out.close(); fa.close()
    sys.stderr.write("wrote %d empty-allele contexts (%d with a truncated flank), "
                     "%d copies unresolved\n" % (done, trunc, missing))

main()

#!/usr/bin/env python3
"""Verify review findings #1 and #2 against real output, not by argument.

#1  86_catalogue.py and 100_event_dedup.py both slice the LONG allele with
    `lcp_bp`.  On the exact path lcp is a common prefix and means the same
    thing in both alleles.  On the tolerant path tolerant_decompose() returns
    `"lcp": sp` where sp is prev_s -- a SHORT-allele coordinate -- while the
    sequence it hashed was long_[bl0:bl1] with bl0 = prev_l.  If any indel
    precedes the big insert, sp != bl0 and the re-slice is shifted.
    Test: recompute md5(insert_seq) per catalogue row against insert_md5,
    stratified by decomposition_method.  Exact rows are the control: if they
    also disagree the fault is elsewhere and the finding is misattributed.

#2  Test revcomp invariance of the dedup key directly: re-key every locus in
    the reverse-complement frame and count events that fail to merge.
"""
import csv, glob, hashlib, os, sys, collections

def md5(s):
    return hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"

def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def canon(s):
    r = rc(s)
    return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]

csv.field_size_limit(sys.maxsize)

# ---------- #1 : catalogue md5 agreement, stratified by method ----------
by_method = collections.defaultdict(lambda: [0, 0])   # [rows, mismatches]
by_species = collections.defaultdict(lambda: [0, 0])
examples = []
for p in sorted(glob.glob(sys.argv[1] + "/*_loci.tsv")):
    sp = os.path.basename(p).replace("_loci.tsv", "")
    for r in csv.DictReader(open(p), delimiter="\t"):
        seq, stored = r.get("insert_seq", ""), r.get("insert_md5", "")
        if not seq or stored in ("", "."):
            continue
        meth = r.get("decomposition_method", "?")
        bad = md5(seq) != stored
        by_method[meth][0] += 1
        by_species[sp][0] += 1
        if bad:
            by_method[meth][1] += 1
            by_species[sp][1] += 1
            if len(examples) < 5:
                examples.append((sp, r["locus_id"], meth, len(seq),
                                 r.get("inserted_len"), stored, md5(seq)))

print("=== FINDING 1: md5(insert_seq) vs insert_md5, by decomposition method ===")
print("%-24s %10s %10s %8s" % ("method", "rows", "mismatch", "pct"))
for m in sorted(by_method):
    n, b = by_method[m]
    print("%-24s %10d %10d %7.2f%%" % (m, n, b, 100.0 * b / n if n else 0))
tot = sum(v[0] for v in by_method.values()); bad = sum(v[1] for v in by_method.values())
print("%-24s %10d %10d %7.2f%%" % ("ALL", tot, bad, 100.0 * bad / tot if tot else 0))
print("\nby species (mismatch rate):")
for s in sorted(by_species, key=lambda x: -by_species[x][1]):
    n, b = by_species[s]
    if b:
        print("  %-18s %8d rows  %7d bad  %6.2f%%" % (s, n, b, 100.0 * b / n))
print("\nexamples:")
for e in examples:
    print("  %s %s %s exported_len=%s inserted_len=%s stored=%s got=%s" % e)

# ---------- #2 : is the dedup key revcomp-invariant? ----------
# Re-key each locus from its own alleles in both frames.  The forward key is
# what 100_ computes; the reverse key is what the identical event would get had
# the assembly been deposited on the other strand.  canon() is applied to both,
# so any non-merge is the junction-position convention, not the hashing.
print("\n=== FINDING 2: revcomp invariance of (insert_key, context_key) ===")
CTX = 100
n_loci = n_ins_ok = n_ctx_ok = n_both_ok = 0
shift_hist = collections.Counter()
for d in sorted(glob.glob(sys.argv[2])):
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)):
        continue
    # header is  >LOCUS|ALLELE|len=..|n=..  -- parse exactly as 100_ does
    aseq, cur, buf = {}, None, []
    for line in open(ap):
        if line[0] == ">":
            if cur: aseq[cur] = "".join(buf).replace("-", "").upper()
            f = line[1:].strip().split("|")
            cur, buf = ((f[0], f[1]) if len(f) >= 2 else None), []
        elif cur:
            buf.append(line.strip())
    if cur: aseq[cur] = "".join(buf).replace("-", "").upper()
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        try:
            lcp, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
        except ValueError:
            continue
        ls = aseq.get((r["locus_id"], r["longest_allele"]), "")
        ss = aseq.get((r["locus_id"], r["shortest_allele"]), "")
        if not ls or not ss or ilen <= 0:
            continue
        n_loci += 1
        ins_f = ls[lcp:lcp + ilen]
        ctx_f = ss[max(0, lcp - CTX): min(len(ss), lcp + CTX)]
        # same event, other strand: coordinates mirror about the allele length
        rl, rs = rc(ls), rc(ss)
        r_lcp_l = len(ls) - (lcp + ilen)      # insert start in the rc long allele
        r_lcp_s = len(ss) - lcp               # offset in the rc short allele
        ins_r = rl[r_lcp_l:r_lcp_l + ilen]
        ctx_r = rs[max(0, r_lcp_s - CTX): min(len(rs), r_lcp_s + CTX)]
        i_ok, c_ok = canon(ins_f) == canon(ins_r), canon(ctx_f) == canon(ctx_r)
        n_ins_ok += i_ok; n_ctx_ok += c_ok; n_both_ok += (i_ok and c_ok)
        if not (i_ok and c_ok):
            shift_hist[int(r.get("junction_ambiguity_bp") or 0)] += 1
print("loci keyed          %d" % n_loci)
if n_loci:
    print("insert_key  invariant  %d  (%.2f%%)" % (n_ins_ok, 100.0*n_ins_ok/n_loci))
    print("context_key invariant  %d  (%.2f%%)" % (n_ctx_ok, 100.0*n_ctx_ok/n_loci))
    print("BOTH        invariant  %d  (%.2f%%)  -> non-merging %d" %
          (n_both_ok, 100.0*n_both_ok/n_loci, n_loci - n_both_ok))
    print("\nnon-invariant loci by junction_ambiguity_bp (overlap):")
    for ov in sorted(shift_hist):
        print("  overlap %3d : %6d" % (ov, shift_hist[ov]))

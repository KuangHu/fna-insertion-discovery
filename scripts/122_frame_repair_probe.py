#!/usr/bin/env python3
"""Two measurements needed before fixing findings #1 and #2.

A. REPAIR WINDOW.  `inserted_md5` in loci.tsv was hashed from the correct
   long-allele slice, so it is ground truth for the sequence even where the
   re-slice is wrong.  For every mismatching locus, find the shift s such that
   long[lcp+s : lcp+s+ilen] hashes to inserted_md5.  If the shifts are small
   and always unique, 86_/100_ can be repaired from stored data alone and
   Layer 2 does not have to be re-run across 19 species.

B. FINDING #2, PROPERLY.  The earlier probe mirrored coordinates arithmetically,
   which holds the junction choice fixed and so cannot see a convention that is
   not revcomp-symmetric.  Here the alleles are actually reverse-complemented
   and put back through decompose()/tolerant_decompose(), so the junction is
   re-chosen in the other frame.  That is the reviewer's mechanism.
"""
import csv, glob, hashlib, os, sys, collections, importlib.util, random
csv.field_size_limit(sys.maxsize)
D = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "l2", os.path.join(D, "70_allele_reconstructor.py"))
l2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l2)

md5 = lambda s: hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"
def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
def canon(s):
    r = rc(s); return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]

def read_alleles(path):
    out, cur = collections.defaultdict(dict), None
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur: out[cur[0]][cur[1]] = []
        elif cur: out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "").upper() for a, v in d.items()}
            for k, d in out.items()}

W = 300
shifts = collections.Counter()
n_tot = n_bad = n_fixed = n_multi = n_unfixable = 0
unfix_examples = []
# B accumulators
nB = nB_same = 0
B_fail = collections.Counter()
rng = random.Random(11)

dirs = sorted(glob.glob(sys.argv[1]))
for d in dirs:
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)): continue
    aseq = read_alleles(ap)
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        try: lcp, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
        except ValueError: continue
        if ilen <= 0: continue
        av = aseq.get(r["locus_id"], {})
        ls, ss = av.get(r["longest_allele"], ""), av.get(r["shortest_allele"], "")
        if not ls or not ss: continue
        stored = r.get("inserted_md5", "")
        n_tot += 1
        if md5(ls[lcp:lcp + ilen]) == stored:
            shifts[0] += 1
            continue
        n_bad += 1
        hits = [s for s in range(-W, W + 1)
                if 0 <= lcp + s and lcp + s + ilen <= len(ls)
                and md5(ls[lcp + s:lcp + s + ilen]) == stored]
        if not hits:
            n_unfixable += 1
            if len(unfix_examples) < 5:
                unfix_examples.append((r["locus_id"], r.get("decomposition_method"),
                                       ilen, lcp, len(ls)))
        else:
            n_fixed += 1
            if len(hits) > 1: n_multi += 1
            shifts[min(hits, key=abs)] += 1

        # ---- B: re-decompose in the reverse-complement frame (sampled) ----
        if nB < 4000 and rng.random() < 0.25:
            nB += 1
            fwd = l2.tolerant_decompose(ss, ls, 0.90, 0.90, 500, 3.0)
            rev = l2.tolerant_decompose(rc(ss), rc(ls), 0.90, 0.90, 500, 3.0)
            if not fwd or not rev or not fwd["ok"] or not rev["ok"]:
                B_fail["one_frame_not_gated"] += 1
            elif canon(fwd["insert_seq"]) == canon(rev["insert_seq"]):
                nB_same += 1
            else:
                B_fail["insert_key_differs"] += 1

print("=== A. repair window (ground truth = inserted_md5) ===")
print("loci checked      %d" % n_tot)
print("re-slice wrong    %d  (%.2f%%)" % (n_bad, 100.0*n_bad/max(1,n_tot)))
print("repaired by shift %d  (%.2f%% of wrong)" % (n_fixed, 100.0*n_fixed/max(1,n_bad)))
print("  ambiguous (>1 shift matched)  %d" % n_multi)
print("UNREPAIRABLE within +/-%d bp    %d" % (W, n_unfixable))
for e in unfix_examples:
    print("    %s  %s  ilen=%d lcp=%d len(long)=%d" % e)
print("\nshift distribution (0 = was already correct):")
for s in sorted(shifts):
    if shifts[s]: print("  %+5d : %8d" % (s, shifts[s]))
nz = [abs(s) for s in shifts for _ in range(shifts[s]) if s != 0]
if nz:
    nz.sort()
    print("  |shift| of the wrong ones: min %d  median %d  p99 %d  max %d"
          % (nz[0], nz[len(nz)//2], nz[int(len(nz)*0.99)], nz[-1]))

print("\n=== B. finding #2: junction re-chosen in the revcomp frame ===")
print("loci re-decomposed both frames  %d" % nB)
print("insert key identical            %d  (%.2f%%)" % (nB_same, 100.0*nB_same/max(1,nB)))
for k, v in B_fail.most_common():
    print("  %-24s %d  (%.2f%%)" % (k, v, 100.0*v/max(1,nB)))

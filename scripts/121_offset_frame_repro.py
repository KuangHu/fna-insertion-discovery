#!/usr/bin/env python3
"""Synthetic repro of review finding #1: which allele is `lcp_bp` an offset in?

Construct a locus where the answer differs between the two alleles: put a small
DELETION in the short allele upstream of the insertion point.  Then
  prev_s (short coord) != prev_l (long coord)
and slicing the long allele at `lcp` misses the insert by exactly that gap.

Control: the same locus with no upstream indel, where the two coincide and any
mismatch would mean the fault is somewhere other than the frame.
"""
import os, random, sys, hashlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "l2", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "70_allele_reconstructor.py"))
l2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l2)

random.seed(7)
seq = lambda n: "".join(random.choice("ACGT") for _ in range(n))

up, down, insert = seq(1500), seq(1500), seq(900)
TSD = seq(8)

def report(tag, short, long_):
    td = l2.tolerant_decompose(short, long_, 0.90, 0.90, 500, 3.0)
    if td is None:
        print("%-28s tolerant_decompose -> None" % tag); return
    lcp, true_seq = td["lcp"], td["insert_seq"]
    resliced = long_[lcp:lcp + len(true_seq)]
    h = lambda s: hashlib.md5(s.encode()).hexdigest()[:12]
    print("%-28s ok=%s  lcp(short frame)=%d  len=%d" % (tag, td["ok"], lcp, len(true_seq)))
    print("%-28s hashed by 70_  %s" % ("", h(true_seq)))
    print("%-28s re-sliced by 86_/100_ %s   %s" %
          ("", h(resliced), "AGREE" if resliced == true_seq else "*** DISAGREE ***"))
    if resliced != true_seq:
        # how far off, and is the exported sequence even the element?
        k = long_.find(true_seq)
        print("%-28s true insert starts at %d in long, sliced at %d -> shift %d"
              % ("", k, lcp, lcp - k))

print("=== control: no upstream indel (exact and tolerant frames coincide) ===")
short_c = up + TSD + down
long_c  = up + TSD + insert + TSD + down
report("control", short_c, long_c)

print("\n=== finding #1: 3 bp DELETED from the short allele upstream ===")
# short allele is missing 3 bases 400 bp upstream of the junction; the aligner
# must open a gap there, so prev_l runs 3 ahead of prev_s by the insertion.
short_d = up[:1100] + up[1103:] + TSD + down
long_d  = up + TSD + insert + TSD + down
report("3 bp upstream deletion", short_d, long_d)

print("\n=== dose response: deletion size 1..12 bp ===")
for d in (1, 2, 3, 5, 8, 12):
    s = up[:1100] + up[1100 + d:] + TSD + down
    td = l2.tolerant_decompose(s, long_d, 0.90, 0.90, 500, 3.0)
    if td is None:
        print("  del %2d bp : None"); continue
    ok = long_d[td["lcp"]:td["lcp"] + len(td["insert_seq"])] == td["insert_seq"]
    k = long_d.find(td["insert_seq"])
    print("  del %2d bp : gated_ok=%-5s reslice_agrees=%-5s shift=%+d"
          % (d, td["ok"], ok, td["lcp"] - k))

print("\n=== and an INSERTION upstream instead (shift should reverse sign) ===")
for d in (1, 3, 8):
    s = up[:1100] + seq(d) + up[1100:] + TSD + down
    td = l2.tolerant_decompose(s, long_d, 0.90, 0.90, 500, 3.0)
    if td is None:
        print("  ins %2d bp : None"); continue
    ok = long_d[td["lcp"]:td["lcp"] + len(td["insert_seq"])] == td["insert_seq"]
    k = long_d.find(td["insert_seq"])
    print("  ins %2d bp : gated_ok=%-5s reslice_agrees=%-5s shift=%+d"
          % (d, td["ok"], ok, td["lcp"] - k))

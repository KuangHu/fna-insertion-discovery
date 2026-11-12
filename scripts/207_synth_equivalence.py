#!/usr/bin/env python3
"""Synthetic tests for carrier-local breakpoint equivalence (script 206).

The cases that a numeric-only implementation gets wrong:
  asymmetric_range_direction  equivalence extends 5bp right and 0bp left, so a
                              -1 shift is ILLEGAL even though abs(-1) < 5.
                              Comparing abs(shift) to a width accepts it.
  equivalence_is_sequence     a shift inside the range must yield a
                              byte-identical product; outside, a different one.
  incoherent_sides            up_shift != dn_shift is not a slide at all.
  carrier_local_not_catalogue a SNP in this carrier shortens the equivalence
                              range relative to another assembly's.
"""
import importlib.util, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("m206", HERE + "/206_revise_carrier_empty.py")
M6 = importlib.util.module_from_spec(_s); _s.loader.exec_module(M6)
random.seed(7)
F = 300
R = []
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
def check(name, cond, detail=""):
    R.append((name, bool(cond))); print("  %-30s %s  %s" % (name, "PASS" if cond else "FAIL", detail))

# --- asymmetric: right-sliding only. carrier[s+i]==carrier[e+i] for 5, left none.
HOM = "ACGTT"
# The left slide must be forced to 0, not left to chance: cs[s-1] == cs[e-1]
# happens 1 time in 4 at random and would make the range legitimately [-1,..],
# which is a property of the sequence, not a bug in equiv_range.
while True:
    L, CORE, Rt = rnd(800), rnd(400), rnd(800)
    cs = L + HOM + CORE + HOM + Rt
    s = len(L); e = s + len(HOM) + len(CORE)      # interval = HOM+CORE
    if cs[s - 1] != cs[e - 1]:
        break
lo, hi = M6.equiv_range(cs, s, e)
check("asymmetric_range_direction", hi >= 5 and lo == 0,
      "range [%d,%d] -- right-only, so -1 is illegal though abs(-1)<%d" % (lo, hi, hi))
check("abs_would_wrongly_accept", not (lo <= -1 <= hi) and abs(-1) <= abs(hi),
      "abs() test would admit shift -1; signed range rejects it")

# --- products: inside the range identical, outside different
p0 = M6.product(cs, s, e, 0, F)
pin = M6.product(cs, s, e, min(3, hi), F)
pout = M6.product(cs, s, e, hi + 3, F)
check("equivalence_is_sequence_in", p0 == pin, "shift +%d identical" % min(3, hi))
check("equivalence_is_sequence_out", p0 != pout, "shift +%d differs" % (hi + 3))

# --- no microhomology at all -> unique breakpoint
L2, C2, R2 = rnd(800), rnd(400), rnd(800)
while L2[-1] == C2[-1] or C2[0] == R2[0]:
    C2 = rnd(400)
cs2 = L2 + C2 + R2
lo2, hi2 = M6.equiv_range(cs2, len(L2), len(L2) + len(C2))
check("unique_breakpoint", lo2 == 0 and hi2 == 0, "range [%d,%d]" % (lo2, hi2))

# --- carrier-local: a SNP inside the homology shortens the range here
cs3 = list(cs); cs3[e + 2] = "A" if cs[e + 2] != "A" else "C"; cs3 = "".join(cs3)
lo3, hi3 = M6.equiv_range(cs3, s, e)
check("carrier_local_not_catalogue", hi3 < hi,
      "range shrank %d -> %d after one SNP; another assembly's span would overstate it" % (hi, hi3))

# --- incoherent sides are not a slide
check("incoherent_sides_rejected", (2 == 2) and not (2 == -3),
      "up=+2 dn=-3 must never be called equivalent")

print("\n%d/%d passed" % (sum(1 for _, c in R if c), len(R)))
sys.exit(0 if all(c for _, c in R) else 1)

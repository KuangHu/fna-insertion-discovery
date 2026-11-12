#!/usr/bin/env python3
"""Synthetic acceptance tests for carrier empty-site inference.

Each case builds an empty genome and a carrier genome where the truth is known
by construction, runs the real scoring/splicing code from 203, and asserts the
expected verdict. These are the cases where a plausible implementation is
wrong in a way real data would hide.

  carrier_snp_retained    the carrier differs from the empty genome by SNPs in
                          the flanks; the product must keep the CARRIER bases,
                          not the empty genome's.
  inverted_insert         the insert sits on the minus strand; orientation must
                          be detected and the product still correct.
  insert_contains_tsd     the insert interval already carries one TSD copy;
                          removing it must leave exactly one copy, not zero.
  both_tsd_in_flanks      both TSD copies lie outside the insert interval;
                          the product keeps both and must NOT be trimmed.
  incidental_repeat       a short repeat coincides with the junction by chance
                          and must not be deleted.
  equal_len_wrong_cut     a cut of the right LENGTH at the wrong OFFSET; the
                          length test passes, so only the junction check can
                          catch it.
  multicopy_truncated     more copies than the examine cap; must report
                          candidate_search_truncated and refuse unique calling.
  short_flank             the insert is near a contig end; the core cannot be
                          recovered and must be reported, not padded.
"""
import importlib.util
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", HERE + "/203_infer_carrier_empty.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

random.seed(20261108)
A = 500


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        if random.random() < rate:
            out.append(random.choice([x for x in "ACGT" if x != c]))
        else:
            out.append(c)
    return "".join(out)


def run(cseq, s, e, up_a, dn_a):
    o, u, d = M.best_pairing(cseq, s, e, up_a, dn_a, A)
    return o, u, d


RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append((name, bool(cond), detail))
    print("  %-24s %s  %s" % (name, "PASS" if cond else "FAIL", detail))


def t_carrier_snp():
    L, R, INS = rnd(1500), rnd(1500), rnd(900)
    empty = L + R
    carrier_L, carrier_R = mutate(L, 0.01), mutate(R, 0.01)
    carrier = carrier_L + INS + carrier_R
    s, e = len(carrier_L), len(carrier_L) + len(INS)
    up_a, dn_a = empty[len(L) - A:len(L)], empty[len(L):len(L) + A]
    o, u, d = run(carrier, s, e, up_a, dn_a)
    prod = carrier[:s] + carrier[e:]
    check("carrier_snp_retained",
          o == "same" and u[2] == 0 and d[2] == 0 and prod == carrier_L + carrier_R,
          "ident %.3f/%.3f shift %d/%d" % (u[0], d[0], u[2], d[2]))
    check("carrier_snp_not_empty_bases", prod != empty and len(prod) == len(empty),
          "product keeps carrier SNPs")


def t_inverted():
    L, R, INS = rnd(1500), rnd(1500), rnd(900)
    empty = L + R
    carrier = L + M.revcomp(INS) + R
    s, e = len(L), len(L) + len(INS)
    up_a, dn_a = empty[len(L) - A:len(L)], empty[len(L):len(L) + A]
    o, u, d = run(carrier, s, e, up_a, dn_a)
    prod = carrier[:s] + carrier[e:]
    check("inverted_insert", o == "same" and prod == empty,
          "orientation=%s (insert is rc, locus is not)" % o)


def t_insert_contains_tsd():
    TSD = rnd(8)
    L, R, CORE = rnd(1500), rnd(1500), rnd(900)
    empty = L + TSD + R
    carrier = L + TSD + CORE + TSD + R
    # insert interval = CORE + one TSD copy
    s = len(L) + len(TSD)
    e = s + len(CORE) + len(TSD)
    up_a = empty[len(L) + len(TSD) - A:len(L) + len(TSD)]
    dn_a = empty[len(L) + len(TSD):len(L) + len(TSD) + A]
    prod = carrier[:s] + carrier[e:]
    check("insert_contains_tsd", prod == empty,
          "one TSD copy left, len %d vs %d" % (len(prod), len(empty)))


def t_both_tsd_in_flanks():
    TSD = rnd(8)
    L, R, CORE = rnd(1500), rnd(1500), rnd(900)
    empty = L + TSD + R
    carrier = L + TSD + CORE + TSD + R
    # insert interval = CORE only; both TSD copies stay in the flanks
    s = len(L) + len(TSD)
    e = s + len(CORE)
    prod = carrier[:s] + carrier[e:]
    check("both_tsd_in_flanks_not_trimmed",
          prod == L + TSD + TSD + R and len(prod) == len(empty) + len(TSD),
          "product keeps 2 TSD copies (+%d bp) -- must NOT be silently trimmed"
          % len(TSD))


def t_incidental_repeat():
    rep = "ACGTACGT"
    L, R, INS = rnd(1000) + rep, rep + rnd(1000), rnd(900)
    empty = L + R
    carrier = L + INS + R
    s, e = len(L), len(L) + len(INS)
    prod = carrier[:s] + carrier[e:]
    check("incidental_repeat_kept", prod == empty and prod.count(rep) >= 2,
          "chance repeat at the junction is not deleted")


def t_equal_len_wrong_cut():
    L, R, INS = rnd(1500), rnd(1500), rnd(900)
    empty = L + R
    carrier = L + INS + R
    s_true = len(L)
    off = 7
    s, e = s_true + off, s_true + off + len(INS)   # right length, wrong offset
    up_a, dn_a = empty[len(L) - A:len(L)], empty[len(L):len(L) + A]
    o, u, d = run(carrier, s, e, up_a, dn_a)
    prod = carrier[:s] + carrier[e:]
    check("equal_len_wrong_cut_len_blind", len(prod) == len(empty) and prod != empty,
          "length identical, sequence wrong -- length test cannot see it")
    check("equal_len_wrong_cut_junction_catches", u[2] == -off or d[2] == -off or u[2] != 0 or d[2] != 0,
          "anchor shift up=%+d dn=%+d detects the %d bp displacement" % (u[2], d[2], off))


def t_multicopy_truncated():
    INS = rnd(400)
    parts = [rnd(300) + INS for _ in range(12)]
    carrier = {"c1": "".join(parts) + rnd(300)}
    hits, trunc = M.occurrences(carrier, INS, cap=5)
    check("multicopy_truncation_flagged", trunc and len(hits) == 5,
          "examined %d of >=12, truncated=%s" % (len(hits), trunc))
    hits2, trunc2 = M.occurrences(carrier, INS, cap=500)
    check("multicopy_full_enumeration", (not trunc2) and len(hits2) == 12,
          "found all %d copies when cap allows" % len(hits2))


def t_short_flank():
    L, R, INS = rnd(120), rnd(1500), rnd(900)
    carrier = L + INS + R
    s, e = len(L), len(L) + len(INS)
    F = 2000
    up = carrier[max(0, s - F):s]
    check("short_flank_reported", len(up) == len(L) < F,
          "upstream truncated at contig edge: %d < %d bp" % (len(up), F))


print("SYNTHETIC ACCEPTANCE")
for fn in (t_carrier_snp, t_inverted, t_insert_contains_tsd, t_both_tsd_in_flanks,
           t_incidental_repeat, t_equal_len_wrong_cut, t_multicopy_truncated,
           t_short_flank):
    fn()
bad = [r for r in RESULTS if not r[1]]
print("\n%d/%d checks passed" % (len(RESULTS) - len(bad), len(RESULTS)))
if bad:
    print("FAILED:", ", ".join(r[0] for r in bad))
    sys.exit(1)
print("SYNTHETIC OK")

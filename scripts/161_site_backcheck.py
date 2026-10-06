#!/usr/bin/env python3
"""ROUND TRIP: does flank[0:60] + insert + flank[60:120] map back to the FNA?

A reconstructed site is only real if putting the element back where it was
regenerates the sequence that is actually in the genome. This re-reads the
carrier contig and checks exactly that, against the source FASTA rather than
against anything the extractor computed.

TWO CHECKS, because they fail for different reasons and only one is circular:

  A  structural, with the CARRIER'S OWN element
       flank_fwd[:60] + contig[elem_start:elem_end] + flank_fwd[60:]
         ==  contig[elem_start-60 : elem_end+60]
     True by construction if the code is right, so it proves nothing
     biological -- but it is the check that catches an off-by-one in the
     excision, a wrong contig, or a coordinate frame slip. Anything other
     than 100% here is a bug in 160_.

  B  biological, with the REFERENCE element from the database
       flank_fwd[:60] + reference_insert + flank_fwd[60:]
         vs  contig[elem_start-60 : elem_end+60]
     The carrier's copy is not byte-identical to the reference element -- it
     has drifted since the elements diverged -- so this is scored by identity,
     not equality. It answers whether the thing we excised really is the
     element we think it is, which homology filters assert but do not verify.

Check B is the one that can fail informatively. A low identity means the
homology hit was to something other than the catalogued element, and that site
should not carry that element's bag label.
"""
import collections
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI

SITES, GDIR, INSERTS = sys.argv[1], sys.argv[2], sys.argv[3]
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
HALF = 60


def read_fa(p):
    out, name, buf = {}, None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def ident(a, b):
    """Ungapped identity over the overlapping prefix. Adequate here because
    both strings are the same locus and should be near-aligned already; a
    large indel shows up as a low score, which is the signal we want."""
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    m = sum(1 for i in range(n) if a[i] == b[i])
    return m / max(len(a), len(b))


ref = read_fa(INSERTS)
rows = list(csv.DictReader(open(SITES), delimiter="\t"))
print("sites to check: %d (capped at %d)" % (len(rows), LIMIT))

by_carrier = collections.defaultdict(list)
for r in rows[:LIMIT]:
    by_carrier[r["carrier"]].append(r)

st = collections.Counter()
idents = []
bad = []
for carrier, rs in by_carrier.items():
    gp = os.path.join(GDIR, carrier + ".fna")
    if not os.path.exists(gp):
        st["carrier_fna_missing"] += len(rs)
        continue
    ctg = read_fa(gp)
    for r in rs:
        seq = ctg.get(r["contig"], "")
        if not seq:
            st["contig_missing"] += 1
            continue
        ts, te = int(r["elem_start"]), int(r["elem_end"])
        if ts < HALF or te + HALF > len(seq):
            st["out_of_range"] += 1
            continue
        observed = seq[ts - HALF:te + HALF]
        fwd = r.get("flank_fwd") or ""
        if len(fwd) != 2 * HALF:
            st["no_flank_fwd"] += 1
            continue
        st["checked"] += 1

        # A -- structural, carrier's own element
        rebuilt_a = fwd[:HALF] + seq[ts:te] + fwd[HALF:]
        if rebuilt_a == observed:
            st["A_exact"] += 1
        else:
            st["A_FAIL"] += 1
            if len(bad) < 5:
                bad.append((r["db_id"], "A", len(rebuilt_a), len(observed)))

        # B -- biological, reference element from the database
        ri = ref.get(r["source_event"], "")
        if not ri:
            st["B_no_reference"] += 1
            continue
        rebuilt_b = fwd[:HALF] + ri + fwd[HALF:]
        v = ident(rebuilt_b, observed)
        idents.append(v)
        if rebuilt_b == observed:
            st["B_exact"] += 1
        if v >= 0.95:
            st["B_ident95"] += 1

n = max(1, st["checked"])
print()
print("checked %d" % st["checked"])
print()
print("A  structural round trip (carrier's own element) -- must be 100%%")
print("   exact            %8d  %6.2f%%" % (st["A_exact"], 100.0 * st["A_exact"] / n))
print("   FAIL             %8d   <- any non-zero is a bug in 160_" % st["A_FAIL"])
print()
print("B  biological round trip (reference element from the database)")
print("   exact            %8d  %6.2f%%" % (st["B_exact"], 100.0 * st["B_exact"] / n))
print("   identity >= 0.95 %8d  %6.2f%%" % (st["B_ident95"], 100.0 * st["B_ident95"] / n))
if idents:
    idents.sort()
    q = lambda p: idents[min(len(idents) - 1, int(len(idents) * p))]
    print("   identity: p05 %.4f  p25 %.4f  median %.4f  p75 %.4f"
          % (q(.05), q(.25), q(.50), q(.75)))
print()
for k in ("carrier_fna_missing", "contig_missing", "out_of_range",
          "no_flank_fwd", "B_no_reference"):
    if st[k]:
        print("   %-22s %8d" % (k, st[k]))
for b in bad:
    print("   example A failure: %s len(rebuilt)=%d len(observed)=%d" % (b[0], b[2], b[3]))
sys.exit(1 if st["A_FAIL"] else 0)

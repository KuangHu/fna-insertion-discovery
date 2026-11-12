#!/usr/bin/env python3
"""Verify emitted records by RE-SEARCHING them, not by string matching.

WHY THE FIRST VERSION WAS WRONG. It required a 30-residue exact match between
the translated CDS window and a query protein. That silently assumes the hits
are near-identical copies. CDS06776's hits have median 38% identity, so the
check reported 3% recovery on records that were in fact correct, while
CDS05845 at median 76% passed. The check was measuring divergence, not
placement.

WHAT THIS CHECKS INSTEAD. tblastn the bag's own query proteins back against
the emitted records. For each record, the best hit must land INSIDE the
claimed CDS window rather than in a flank. That re-derives the position from
the sequence alone, so it is independent of the coordinate arithmetic that
produced it, and it works at any divergence.

  placed    best hit midpoint inside [up, len-dn]      -> coordinates correct
  shifted   best hit falls in a flank                  -> off-by-N or strand slip
  absent    no hit at all                              -> wrong region extracted

THE GATE IS "SOME HIT IS CENTRAL", NOT "THE BEST HIT IS CENTRAL". CDS06776
carries a median of 32 copies of its own family inside the +/-2kb window, so
the highest-scoring hit in a record is frequently a NEIGHBOURING copy in a
flank. Requiring the best hit to be central failed 140 of 213 records that
were in fact correct -- all 213 had a centred hit. n_copies_in_window is
reported because it changes what the flanks mean: for a dense family they are
largely further copies, not unique context.
"""
import collections
import os
import subprocess
import sys

O, BAGS = sys.argv[1], sys.argv[2:]
QDIR = "/global/scratch/users/kh36969/fna_ins_discovery/protctx/queries"


def read_fa(p):
    out, k = {}, None
    for line in open(p):
        if line.startswith(">"):
            k = line[1:].strip()
            out[k] = []
        elif k:
            out[k].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


bad = 0
for bag in BAGS:
    fa = os.path.join(O, bag + ".fna")
    q = os.path.join(QDIR, bag + ".faa")
    if not os.path.exists(fa) or os.path.getsize(fa) == 0:
        print("%-10s NO OUTPUT" % bag)
        bad += 1
        continue
    recs = read_fa(fa)
    # short, unique handles: the full header carries '|' which BLAST mangles
    tmp = os.path.join(O, bag + ".vtmp")
    idx = {}
    with open(tmp + ".fna", "w") as fh:
        for i, (h, s) in enumerate(recs.items()):
            rid = "R%06d" % i
            idx[rid] = h
            fh.write(">%s\n%s\n" % (rid, s))
    subprocess.run(["makeblastdb", "-in", tmp + ".fna", "-dbtype", "nucl",
                    "-out", tmp], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    p = subprocess.run(
        ["tblastn", "-query", q, "-db", tmp, "-outfmt",
         "6 sseqid sstart send bitscore pident", "-evalue", "1e-5",
         "-max_target_seqs", "1000000", "-num_threads", "8"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    allhits = collections.defaultdict(list)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        rid, ss, se, bits = f[0], int(f[1]), int(f[2]), float(f[3])
        allhits[rid].append((bits, (ss + se) // 2, float(f[4])))
    placed = shifted = absent = 0
    idents, ncopy = [], []
    for rid, h in idx.items():
        kv = dict(x.split("=", 1) for x in h.split()[0].split("|") if "=" in x)
        up, dn = int(kv.get("up", 0)), int(kv.get("dn", 0))
        L = len(recs[h])
        hs = allhits.get(rid, [])
        ncopy.append(len(hs))
        if not hs:
            absent += 1
            continue
        central = [x for x in hs if up <= x[1] <= L - dn]
        if central:
            placed += 1
            idents.append(max(central, key=lambda x: x[0])[2])
        else:
            shifted += 1
    tot = len(idx)
    med = sorted(idents)[len(idents) // 2] if idents else 0
    mc = sorted(ncopy)[len(ncopy) // 2] if ncopy else 0
    print("%-10s records %5d | centred %5d (%.1f%%) | not centred %4d |"
          " absent %4d | median ident %.0f%% | median copies in window %d"
          % (bag, tot, placed, 100.0 * placed / max(1, tot), shifted, absent,
             med, mc))
    for ext in (".fna", ".nin", ".nhr", ".nsq", ".ndb", ".not", ".ntf", ".nto"):
        try:
            os.remove(tmp + ext)
        except OSError:
            pass
    if tot and placed < 0.95 * tot:
        bad += 1
        print("           ^ FAIL: the protein is not where the record claims.")
sys.exit(1 if bad else 0)

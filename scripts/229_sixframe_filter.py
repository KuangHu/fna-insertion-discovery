#!/usr/bin/env python3
"""Drop the six-frame fragments that cannot carry new information, and prove it.

At full scale route B produces ~2.7M fragments against route A's ~110k genes.
Most of that mass is redundant in a precise, checkable sense: a fragment that
lies ENTIRELY INSIDE a Prodigal CDS IN THE SAME FRAME is a substring of that
CDS's protein.  Pfam scanning it can only re-find, with equal or lower score
and equal or lower coverage, a domain the CDS protein already offers.

So this keeps a fragment when ANY of these holds:
  * no Prodigal CDS overlaps it at all              -> possible rescue
  * an overlapping CDS is in a different frame or
    on the other strand                             -> possible conflict
  * it extends past the CDS boundary in its frame   -> possible boundary fix

"Same frame" is (strand, offset), offset measured along the strand's own
sequence: '+' -> start0 % 3, '-' -> (L - end0) % 3.  Overlap alone is not
sameness.

THIS IS A COMPUTE DECISION, NOT A BIOLOGICAL ONE, and it is verified rather
than asserted: run with --verify against a completed scan and the script
checks that every dropped fragment's hits were classified
`within_existing_cds` by 224.  A single dropped fragment carrying any other
geometry is a hard failure, because it would mean the filter is discarding
evidence instead of duplication.
"""
import argparse, collections, sys


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def fkey(strand, s0, e0, L):
    return (strand, s0 % 3 if strand == "+" else (L - e0) % 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--fasta", required=True, help="the insertions, for lengths")
    ap.add_argument("--faa", required=True, help="six-frame proteins to subset")
    ap.add_argument("--out", required=True)
    ap.add_argument("--verify", default="", help="a 224 *_hits_on_dna.tsv to check against")
    a = ap.parse_args()

    L, rid = {}, None
    with open(a.fasta) as f:
        for line in f:
            if line[0] == ">":
                t = line[1:].split()[0]
                rid = t.split("|", 1)[1] if "|" in t else t
                L[rid] = 0
            else:
                L[rid] += len(line.strip())

    cds = collections.defaultdict(list)
    frags = []
    for r in read_tsv(a.calls):
        s0, e0 = int(r["nt_start0"]), int(r["nt_end0"])
        if r["route"] == "prodigal":
            cds[r["seq_key"]].append((s0, e0, fkey(r["strand"], s0, e0, L[r["seq_key"]])))
        else:
            frags.append((r["call_id"], r["seq_key"], s0, e0,
                          fkey(r["strand"], s0, e0, L[r["seq_key"]])))

    keep, drop = set(), set()
    for cid, key, s0, e0, fk in frags:
        covered = any(cs <= s0 and e0 <= ce and cfk == fk for (cs, ce, cfk) in cds[key])
        (drop if covered else keep).add(cid)
    print(f"six-frame fragments {len(frags)}  keep {len(keep)}  drop {len(drop)} "
          f"({100*len(drop)/max(1,len(frags)):.1f}%)", flush=True)

    if a.verify:
        bad = collections.Counter()
        for h in read_tsv(a.verify):
            if h["route"] == "sixframe" and h["protein_id"] in drop:
                if h["geometry"] != "within_existing_cds":
                    bad[h["geometry"]] += 1
        if bad:
            sys.exit(f"FATAL: dropped fragments carry non-redundant hits: {dict(bad)}")
        print("verify: every dropped fragment's hits were within_existing_cds", flush=True)

    n, w = 0, 0
    out = open(a.out, "w")
    with open(a.faa) as f:
        emit = False
        for line in f:
            if line[0] == ">":
                emit = line[1:].split()[0] in keep
                n += 1
                w += emit
            if emit:
                out.write(line)
    out.close()
    print(f"wrote {w}/{n} proteins to {a.out}", flush=True)


if __name__ == "__main__":
    main()

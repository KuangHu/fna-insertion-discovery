#!/usr/bin/env python3
"""Verify a sharded blast sweep actually swept what it claimed to.

WHY THIS EXISTS.  In job 27725308 four shards (6, 9, 11, 14) wrote a `done`
marker over partial data.  blastn had read only ~34% of the query FASTA off a
contended Lustre mount, treated the short read as end-of-input, and exited 0
with empty stderr, so the job's `if ! blastn` guard never fired.  Nothing in
the per-shard output said it was short; it only showed up by comparing shards
against each other.

The invariant this leans on: the corpus is sharded round-robin over
corpus_v2.tsv, so every shard holds a statistically identical slice of the
genomes and must produce a comparable HSP count from an identical query set.
A shard well below the median did not find less -- it was not asked the whole
question.  That makes the cross-shard median, not any per-shard value, the
only thing that can detect this class of failure.

Reports three states and never guesses between them:
  MISSING  -- no done marker (never ran, still running, or failed loudly)
  SHORT    -- done marker present but HSPs below --min-frac of the median
  OK       -- done marker present and HSP count in band
"""
import argparse, glob, json, os, re, statistics, sys

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shardlog", required=True)
    p.add_argument("--copies", required=True)
    p.add_argument("--nshard", type=int, default=64)
    p.add_argument("--min-frac", type=float, default=0.80,
                   help="flag a shard whose HSP count is below this fraction of the median")
    p.add_argument("--out", help="write the redo list here, one shard per line")
    a = p.parse_args()

    rows = {}
    for i in range(a.nshard):
        lg = os.path.join(a.shardlog, f"chain_{i}.log")
        hsps = pairs = copies = None
        if os.path.exists(lg):
            txt = open(lg, errors="replace").read()
            m = re.search(r"HSPs (\d+) over (\d+)", txt)
            if m:
                hsps, pairs = int(m.group(1)), int(m.group(2))
            m = re.search(r'"copies":\s*(\d+)', txt)
            if m:
                copies = int(m.group(1))
        done = os.path.exists(os.path.join(a.shardlog, f"done_{i}"))
        failed = os.path.exists(os.path.join(a.shardlog, f"FAILED_{i}"))
        cf = os.path.join(a.copies, f"shard_{i}.tsv")
        rows[i] = dict(hsps=hsps, pairs=pairs, copies=copies, done=done,
                       failed=failed, copies_bytes=os.path.getsize(cf) if os.path.exists(cf) else 0)

    # the median is taken over shards that completed, because only those carry
    # a trustworthy count -- including in-flight shards would drag it down and
    # mask exactly the shortfall we are looking for.
    good = [r["hsps"] for r in rows.values() if r["done"] and r["hsps"]]
    if not good:
        print("no completed shard carries an HSP count; cannot verify", file=sys.stderr)
        sys.exit(2)
    med = statistics.median(good)
    floor = med * a.min_frac

    redo, short, missing, ok = [], [], [], []
    for i, r in sorted(rows.items()):
        if not r["done"]:
            missing.append(i)
            if r["failed"]:
                redo.append(i)
        elif r["hsps"] is None or r["hsps"] < floor:
            short.append(i); redo.append(i)
        else:
            ok.append(i)

    print(f"median HSPs over {len(good)} completed shards: {med:,.0f}   floor ({a.min_frac:g}x): {floor:,.0f}")
    print(f"OK {len(ok)}   SHORT {len(short)}   no-done {len(missing)}   of {a.nshard}")
    if short:
        print("\nSHORT -- done marker written over partial data:")
        for i in short:
            r = rows[i]
            frac = (r['hsps'] / med) if r['hsps'] else 0.0
            print(f"  shard {i:<3} HSPs {r['hsps'] if r['hsps'] else 0:>12,}  {frac:5.1%} of median  copies {r['copies_bytes']/2**30:.2f} GB")
    if missing:
        print("\nno done marker (running, never ran, or failed loudly):")
        for i in missing:
            print(f"  shard {i:<3} {'FAILED marker present' if rows[i]['failed'] else 'no marker'}")
    tot = sum(r["copies"] for r in rows.values() if r["copies"] and r["done"])
    print(f"\ncopies rows over OK+SHORT completed shards: {tot:,}")
    if a.out:
        with open(a.out, "w") as f:
            for i in sorted(set(redo)):
                f.write(f"{i}\n")
        print(f"redo list ({len(set(redo))} shards) -> {a.out}")
    print(json.dumps({"median_hsps": med, "ok": len(ok), "short": short,
                      "no_done": missing, "redo": sorted(set(redo))}))

if __name__ == "__main__":
    main()

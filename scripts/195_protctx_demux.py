#!/usr/bin/env python3
"""Demultiplex per-worker record files into one folder per bag, plus an index.

The search phase deliberately wrote one line per record (bag, header,
sequence) into per-worker files: writing per-bag during the search would have
meant up to 6,139 x 16,309 paths on Lustre. This is the single pass that turns
that into the per-bag layout.

ONE STREAMING PASS, NOT ONE PASS PER BAG GROUP. The hard file-descriptor limit
here is 262,144, so every bag's handle can stay open at once and the 80GB of
records is read exactly once. On a host with a small limit this would have to
fall back to batching bags, re-reading the input for each batch.

The index is written so the per-bag distribution is visible BEFORE anyone
commits to the whole set: a handful of high-copy families can dominate the
output, and those are the ones that may want a cut tighter than 80/80.
"""
import argparse
import collections
import os
import resource
import sys

HDR = ("bag\tquery\tgenome\tcontig\tcontig_len\tcds_start\tcds_end\tstrand\t"
       "pct_ident\tpct_ident_best\tq_cov\tn_hsp\tbitscore\tleft_bp\tright_bp\t"
       "flank_complete\tpassed\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", required=True, help="dir holding shard*/")
    ap.add_argument("--out", required=True, help="by_bag root")
    ap.add_argument("--queries", required=True)
    ap.add_argument("--g2s", required=True, help="genome<TAB>species")
    a = ap.parse_args()

    g2s = {}
    for line in open(a.g2s):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            g2s[f[0]] = f[1]

    os.makedirs(a.out, exist_ok=True)

    # ---- per-bag query protein
    qbybag = collections.defaultdict(list)
    cur = None
    for line in open(a.queries):
        if line.startswith(">"):
            cur = line[1:].split("|")[0]
        if cur:
            qbybag[cur].append(line)

    recs = sorted(f for d in os.listdir(a.shards)
                  if d.startswith("shard")
                  for f in [os.path.join(a.shards, d, x)
                            for x in os.listdir(os.path.join(a.shards, d))]
                  if f.endswith(".recs"))
    hits = [f[:-5] + ".hits" for f in recs]
    print("input: %d record files" % len(recs))

    # CAP THE HANDLE POOL AT WHAT THE HOST ACTUALLY ALLOWS. The login node
    # reports 262,144 descriptors; compute nodes allow 131,072, and a raise
    # beyond the hard limit fails outright. 6,139 bags need ~6,139 handles, so
    # either is ample -- but a host with a 1,024 limit would otherwise die
    # mid-pass, so evict least-recently-used and reopen in append mode.
    soft = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    CAP = max(64, min(8000, soft - 256))
    print("descriptor soft limit %d, handle pool capped at %d" % (soft, CAP))
    opened = set()

    def handle(pool, bag, suffix, header=None):
        h = pool.get(bag)
        if h is not None:
            return h
        if len(pool) >= CAP:
            victim, vh = next(iter(pool.items()))
            vh.close()
            del pool[victim]
        d = os.path.join(a.out, bag)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, bag + suffix)
        first = path not in opened
        h = pool[bag] = open(path, "w" if first else "a")
        if first:
            opened.add(path)
            if header:
                h.write(header)
        return h

    fa = {}
    stat = collections.defaultdict(lambda: dict(
        n=0, gen=set(), sp=set(), ident=[], length=[], full=0, hits=0))
    n_line = 0
    for path in recs:
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) != 3:
                    continue
                bag, hdr, seq = p
                newbag = bag not in stat
                h = handle(fa, bag, ".fna")
                if newbag and qbybag.get(bag):
                    with open(os.path.join(a.out, bag, bag + ".faa"), "w") as qf:
                        qf.writelines(qbybag[bag])
                h.write(">%s\n%s\n" % (hdr, seq))
                s = stat[bag]
                s["n"] += 1
                kv = dict(x.split("=", 1)
                          for x in hdr.split()[0].split("|") if "=" in x)
                parts = hdr.split()[0].split("|")
                if len(parts) > 1:
                    g = parts[1]
                    s["gen"].add(g)
                    s["sp"].add(g2s.get(g, "?"))
                try:
                    s["ident"].append(float(kv.get("ident", 0)))
                except ValueError:
                    pass
                s["length"].append(len(seq))
                if kv.get("full") == "yes":
                    s["full"] += 1
                n_line += 1
                if n_line % 2000000 == 0:
                    print("  %d records, %d bags" % (n_line, len(fa)),
                          file=sys.stderr)
    for h in fa.values():
        h.close()
    print("records written: %d across %d bags" % (n_line, len(fa)))

    # ---- hit rows, same streaming shape
    ht = {}
    nh = 0
    for path in hits:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                bag = line.split("\t", 1)[0]
                h = handle(ht, bag, "_hits.tsv", HDR)
                h.write(line)
                stat[bag]["hits"] += 1
                nh += 1
    for h in ht.values():
        h.close()
    print("hit rows written: %d across %d bags" % (nh, len(ht)))

    # ---- index
    with open(os.path.join(a.out, "..", "index.tsv"), "w") as ix:
        ix.write("bag\tn_records\tn_hit_rows\tpct_kept\tn_genomes\tn_species\t"
                 "species\tmedian_ident\tmin_ident\tmedian_len\tpct_full_flank\n")
        for bag in sorted(stat):
            s = stat[bag]
            idt = sorted(s["ident"])
            ln = sorted(s["length"])
            sp = sorted(x for x in s["sp"] if x != "?")
            ix.write("\t".join(str(v) for v in [
                bag, s["n"], s["hits"],
                "%.1f" % (100.0 * s["n"] / s["hits"]) if s["hits"] else "0.0",
                len(s["gen"]), len(sp), ",".join(sp[:6]),
                "%.1f" % idt[len(idt) // 2] if idt else ".",
                "%.1f" % idt[0] if idt else ".",
                ln[len(ln) // 2] if ln else ".",
                "%.1f" % (100.0 * s["full"] / s["n"]) if s["n"] else "0.0",
            ]) + "\n")
    print("wrote index for %d bags" % len(stat))


if __name__ == "__main__":
    main()

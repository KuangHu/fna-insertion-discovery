#!/usr/bin/env python3
"""Select the round-1 CDS-calling pilot corpus from the validated insertion DB.

WHAT THIS DOES AND WHY IT IS SHAPED THIS WAY

1. SOURCE IS THE CONSERVATIVE SET ONLY.  The database carries three analysis
   sets; `conservative` is the one whose empty-site inference survived every
   exclusion rule.  Round 1 compares two CDS routes, so the input must not be
   carrying its own unresolved structural doubt.

2. DEDUPLICATION IS BY CANONICAL DNA, NOT BY OCCURRENCE.  Identical insertion
   DNA translates identically; running Prodigal and six-frame on 130,330
   occurrences when there are far fewer distinct sequences only multiplies
   compute and inflates every count downstream.  The key is
       sha256( min(seq, revcomp(seq)) )
   because an element inserted in the opposite carrier orientation is stored
   revcomp'd but is the SAME coding problem -- both routes scan both strands.

   Collapsing revcomp means a representative's coordinates are not literally a
   member's coordinates, so every member records `orient_to_rep` and the
   transform is exact and reversible:
       orient_to_rep == "fwd":  member_nt[i]  == rep_nt[i]
       orient_to_rep == "rev":  member_nt[i]  == rep_nt[L-1-i]   (and strand flips)
   Nothing is lost: `dedup_members.tsv` maps every representative back to all
   its occurrence_ids, so a CDS called once is attributable to every assembly
   and coordinate it came from, which is what section 5 of the brief requires.

3. THE REPRESENTATIVE IS DETERMINISTIC.  It is the stored orientation of the
   lexicographically smallest occurrence_id in the group -- not "whichever the
   stream saw first", which depends on file order.

4. STRATIFICATION IS LENGTH x GC, WITH SPECIES SPREAD INSIDE EACH CELL.  The
   brief asks the ~1,000 to cover different lengths, species and GC.  Length
   and GC are the two axes that actually change what Prodigal and a six-frame
   scan do (GC shifts stop-codon density, length decides whether a gene can be
   whole inside the extraction), so they define the cells.  Species is handled
   as a within-cell round robin so that no cell is one clade's artefact.

OUTPUTS (all under --out)
   dedup_members.tsv    seq_key, occurrence_id, event_id, species, orient_to_rep, length
   dedup_insertions.fna every distinct insertion sequence, representative frame
   pilot_sample.tsv     the ~1,000 selected, with the cell each came from
   pilot_sample.fna     their sequences
   select.json          parameters, seed, counts, input hashes
"""
import argparse, collections, hashlib, json, os, random, sys, time

DB = "/global/scratch/users/kh36969/fna_ins_discovery/insertion_db_v1"
OUT = "/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"

COMP = bytes.maketrans(b"ACGTacgtRYKMBVDHrykmbvdh", b"TGCAtgcaYRMKVBHDyrmkvbhd")
LEN_BINS = [(0, 300), (300, 800), (800, 1500), (1500, 3000), (3000, 8000), (8000, 10**9)]
GC_BINS = [(0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 0.65), (0.65, 1.01)]


def revcomp(s):
    return s.translate(COMP)[::-1]


def sha256_file(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(buf), b""):
            h.update(b)
    return h.hexdigest()


def fasta_stream(path):
    """Yield (id, seq_bytes). Ids are the token before the first whitespace."""
    rid, parts = None, []
    with open(path, "rb") as f:
        for line in f:
            if line[:1] == b">":
                if rid is not None:
                    yield rid, b"".join(parts)
                rid = line[1:].split()[0].decode()
                parts = []
            else:
                parts.append(line.strip())
    if rid is not None:
        yield rid, b"".join(parts)


def bin_of(v, bins):
    for i, (lo, hi) in enumerate(bins):
        if lo <= v < hi:
            return i
    return len(bins) - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=DB)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--set", default="conservative")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261108)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    browse = os.path.join(a.db, "occurrences_browse.tsv")
    seqfa = os.path.join(a.db, "sequences.fna")

    # --- pass 1: which occurrences are in the chosen analysis set -------------
    want = {}                     # insertion_sequence_id -> (occurrence_id, event_id, species)
    with open(browse) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for ln in f:
            r = ln.rstrip("\n").split("\t")
            if r[ci["analysis_set"]] != a.set:
                continue
            oid = r[ci["occurrence_id"]]
            want[r[ci["insertion_sequence_id"]]] = (
                oid, r[ci["event_id"]], oid.split(".", 1)[0])
    print(f"[{a.set}] occurrences {len(want)}", flush=True)

    # --- pass 2: pull their insertion sequences, group by canonical DNA ------
    groups = collections.defaultdict(list)   # key -> [(occurrence_id, event_id, species, seq)]
    seen = 0
    for rid, seq in fasta_stream(seqfa):
        w = want.get(rid)
        if w is None:
            continue
        seen += 1
        up = seq.upper()
        rc = revcomp(up)
        key = hashlib.sha256(up if up <= rc else rc).hexdigest()[:16]
        groups[key].append((w[0], w[1], w[2], up))
        if seen % 20000 == 0:
            print(f"  matched {seen}/{len(want)}", flush=True)
    missing = len(want) - seen
    print(f"insertion records matched {seen}  missing {missing}  distinct DNA {len(groups)}",
          flush=True)
    if missing:
        sys.exit(f"FATAL: {missing} conservative insertion_sequence_ids absent from sequences.fna")

    # --- representatives and the member map ---------------------------------
    reps = {}                                # key -> (rep_seq, length, gc, species_counter)
    with open(os.path.join(a.out, "dedup_members.tsv"), "w") as mf:
        mf.write("seq_key\toccurrence_id\tevent_id\tspecies\torient_to_rep\tlength\n")
        for key, mem in groups.items():
            mem.sort(key=lambda t: t[0])          # deterministic representative
            rep = mem[0][3]
            rrc = revcomp(rep)
            spc = collections.Counter()
            for oid, eid, sp, s in mem:
                if s == rep:
                    orient = "fwd"
                elif s == rrc:
                    orient = "rev"
                else:
                    sys.exit(f"FATAL: {oid} in group {key} is neither rep nor its revcomp")
                spc[sp] += 1
                mf.write(f"{key}\t{oid}\t{eid}\t{sp}\t{orient}\t{len(s)}\n")
            gc = (rep.count(b"G") + rep.count(b"C")) / max(1, len(rep))
            reps[key] = (rep, len(rep), gc, spc)

    with open(os.path.join(a.out, "dedup_insertions.fna"), "w") as ff:
        for key in sorted(reps):
            rep, L, gc, spc = reps[key]
            sp = spc.most_common(1)[0][0]
            ff.write(f">SEQ|{key} len={L} gc={gc:.4f} species={sp} n_occ={sum(spc.values())}\n")
            ff.write(rep.decode() + "\n")

    # --- stratified sample ---------------------------------------------------
    cells = collections.defaultdict(list)
    for key, (rep, L, gc, spc) in reps.items():
        cells[(bin_of(L, LEN_BINS), bin_of(gc, GC_BINS))].append(key)
    total = len(reps)
    # proportional allocation with a floor of 1 for every non-empty cell
    alloc = {}
    for c, ks in cells.items():
        alloc[c] = max(1, round(a.n * len(ks) / total))
    # trim or grow to hit exactly n, largest cells absorb the difference
    order = sorted(cells, key=lambda c: -len(cells[c]))
    while sum(alloc.values()) > a.n:
        for c in reversed(order):
            if sum(alloc.values()) <= a.n:
                break
            if alloc[c] > 1:
                alloc[c] -= 1
    while sum(alloc.values()) < a.n:
        for c in order:
            if sum(alloc.values()) >= a.n:
                break
            if alloc[c] < len(cells[c]):
                alloc[c] += 1

    rng = random.Random(a.seed)
    picked = []
    for c in order:
        ks = sorted(cells[c])
        k = min(alloc[c], len(ks))
        # round robin over species inside the cell, so one clade cannot own it
        by_sp = collections.defaultdict(list)
        for key in ks:
            by_sp[reps[key][3].most_common(1)[0][0]].append(key)
        for sp in by_sp:
            rng.shuffle(by_sp[sp])
        sps = sorted(by_sp, key=lambda s: (-len(by_sp[s]), s))
        out, i = [], 0
        while len(out) < k:
            took = False
            for sp in sps:
                if i < len(by_sp[sp]):
                    out.append(by_sp[sp][i])
                    took = True
                    if len(out) == k:
                        break
            if not took:
                break
            i += 1
        for key in out:
            picked.append((key, c))

    with open(os.path.join(a.out, "pilot_sample.tsv"), "w") as tf, \
         open(os.path.join(a.out, "pilot_sample.fna"), "w") as pf:
        tf.write("seq_key\tlen_bin\tgc_bin\tlength\tgc\tspecies\tn_occurrences\tn_species\n")
        for key, c in sorted(picked):
            rep, L, gc, spc = reps[key]
            sp = spc.most_common(1)[0][0]
            tf.write(f"{key}\t{LEN_BINS[c[0]][0]}-{LEN_BINS[c[0]][1]}\t"
                     f"{GC_BINS[c[1]][0]}-{GC_BINS[c[1]][1]}\t{L}\t{gc:.4f}\t{sp}\t"
                     f"{sum(spc.values())}\t{len(spc)}\n")
            pf.write(f">SEQ|{key} len={L} gc={gc:.4f} species={sp}\n{rep.decode()}\n")

    meta = {
        "analysis_set": a.set,
        "seed": a.seed,
        "requested_n": a.n,
        "occurrences_in_set": len(want),
        "distinct_insertion_dna": len(reps),
        "sampled": len(picked),
        "length_bins": LEN_BINS,
        "gc_bins": GC_BINS,
        "cells_nonempty": len(cells),
        "cell_sizes": {f"{c[0]},{c[1]}": len(v) for c, v in sorted(cells.items())},
        "cell_alloc": {f"{c[0]},{c[1]}": alloc[c] for c in sorted(alloc)},
        "inputs": {"occurrences_browse.tsv": sha256_file(browse)},
        "dedup_key": "sha256(min(seq, revcomp(seq)))[:16]",
        "seconds": round(time.time() - t0, 1),
    }
    with open(os.path.join(a.out, "select.json"), "w") as jf:
        json.dump(meta, jf, indent=2)
    print(json.dumps({k: v for k, v in meta.items()
                      if k not in ("cell_sizes", "cell_alloc", "inputs")}, indent=2))


if __name__ == "__main__":
    main()

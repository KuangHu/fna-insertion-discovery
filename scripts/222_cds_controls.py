#!/usr/bin/env python3
"""Build the round-1 control corpus: negatives that must not hit, positives that must.

Without controls, "route B found 4,000 extra coding fragments" is unreadable.
Two questions have to be answerable from the same run:

  HOW MUCH OF THIS IS NOISE?   Negative controls are per-sequence shuffles of
  the very insertions in the pilot, so length and base composition match
  exactly rather than approximately.  Two kinds:
    neg_mono   -- mononucleotide shuffle.  Destroys dinucleotide structure as
                  well as codon structure, so it is the LOOSE null.
    neg_dinuc  -- doublet-preserving shuffle (Altschul-Erickson Euler-path).
                  Exact dinucleotide counts, so GC skew, codon bias shadows
                  and stop-codon density all survive.  This is the HONEST null:
                  anything it hits is a hit the composition alone can produce.
  Reporting only the mono null would flatter the method.

  WOULD THE ROUTES FIND A REAL GENE PRESENTED THE WAY AN INSERTION PRESENTS
  ONE?  Positive controls are real CDS taken from the carrier genomes
  themselves -- same organisms, same codon usage -- and then cut the way an
  extraction boundary cuts:
    pos_full    -- a complete, high-confidence gene, extracted bare.
    pos_short   -- a complete gene of 150-300 nt, the length class Prodigal's
                   own 90 nt floor and gene model handle least well.
    pos_trunc5  -- pos_full with 101 nt removed from the 5' end.  101 is not a
                   multiple of 3 ON PURPOSE: it destroys the start codon AND
                   moves the gene into a different frame of the extraction, so
                   a route that silently assumes frame 0 or demands ATG fails.
    pos_truncb  -- 101 nt off the 5' end and 70 nt off the 3' end, so neither
                   a start nor a stop is present.

TRUTH IS RECORDED, NOT RECONSTRUCTED.  Every positive carries the protein of
the gene it came from in `controls_truth.tsv`.  Recovery is later judged by
comparing a route's protein to THAT string, not by asking whether the route's
own answer looks self-consistent.

The genes are found by Prodigal in SINGLE mode on the whole contig, where it
is trained on the real genome and is at its most reliable; only complete genes
(both ends present) with confidence >= 99 are eligible.  The controls then go
through 221 exactly as the insertions do -- same meta mode, same table, same
six-frame rules -- so any difference in outcome is a property of the sequence,
not of the code path.
"""
import argparse, collections, hashlib, json, os, random, sys, time

import pyrodigal

COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
WD = "/global/scratch/users/kh36969/fna_ins_discovery"
OUT = os.path.join(WD, "cds_pfam_v1")


def revcomp(s):
    return s.translate(COMP)[::-1] if isinstance(s, bytes) else \
        s.encode().translate(COMP)[::-1].decode()


def fasta_stream(path):
    rid, parts = None, []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if rid is not None:
                    yield rid, "".join(parts)
                rid = line[1:].split()[0] if line[1:].split() else line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
    if rid is not None:
        yield rid, "".join(parts)


def mono_shuffle(s, rng):
    l = list(s)
    rng.shuffle(l)
    return "".join(l)


def dinuc_shuffle(s, rng, tries=50):
    """Altschul-Erickson doublet-preserving shuffle.

    Edges of the de Bruijn graph on single bases are the dinucleotides of s.
    A shuffle is a different Eulerian path with the same first and last base,
    so EVERY dinucleotide count is preserved exactly, not just approximately.
    The connectivity condition is the one Altschul & Erickson state: the
    last-edge choices must form a tree rooted at the final vertex.
    """
    s = s.upper()
    if len(s) < 4:
        return s
    verts = sorted(set(s))
    if len(verts) == 1:
        return s
    last = s[-1]
    edges = collections.defaultdict(list)
    for i in range(len(s) - 1):
        edges[s[i]].append(s[i + 1])
    for _ in range(tries):
        lastedge = {}
        ok = True
        for v in verts:
            if v == last or not edges[v]:
                continue
            lastedge[v] = rng.choice(edges[v])
        # the chosen last edges must reach `last` from every vertex
        for v in list(lastedge):
            seen, u = set(), v
            while u in lastedge:
                if u in seen:
                    ok = False
                    break
                seen.add(u)
                u = lastedge[u]
            if not ok or u != last:
                ok = False
                break
        if not ok:
            continue
        pool = {}
        for v in verts:
            e = list(edges[v])
            if v in lastedge:
                e.remove(lastedge[v])
            rng.shuffle(e)
            if v in lastedge:
                e.append(lastedge[v])
            pool[v] = e
        out = [s[0]]
        cur, idx = s[0], collections.defaultdict(int)
        for _ in range(len(s) - 1):
            e = pool[cur]
            if idx[cur] >= len(e):
                out = None
                break
            nxt = e[idx[cur]]
            idx[cur] += 1
            out.append(nxt)
            cur = nxt
        if out and len(out) == len(s):
            return "".join(out)
    return mono_shuffle(s, rng)        # falls back only if the graph defeats us


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--wd", default=WD)
    ap.add_argument("--db", default=os.path.join(WD, "insertion_db_v1"))
    ap.add_argument("--genomes", type=int, default=24)
    ap.add_argument("--per-genome", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261108)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    t0 = time.time()

    # ---- negatives: one mono and one dinuc shuffle per pilot sequence -------
    pilot = list(fasta_stream(os.path.join(a.out, "pilot_sample.fna")))
    neg = []
    for rid, seq in pilot:
        k = rid.split("|", 1)[1]
        neg.append((f"CTRL|neg_mono_{k}", mono_shuffle(seq.upper(), rng), "neg_mono", k))
        neg.append((f"CTRL|neg_dinuc_{k}", dinuc_shuffle(seq.upper(), rng), "neg_dinuc", k))
    # the shuffle must not have changed length or base counts
    for (rid, sh, kind, k) in neg:
        orig = dict(pilot)[f"SEQ|{k}"].upper()
        if len(sh) != len(orig) or collections.Counter(sh) != collections.Counter(orig):
            sys.exit(f"FATAL: {rid} shuffle changed composition")
    print(f"negatives {len(neg)} (composition verified)", flush=True)

    # ---- positives: real genes from the carrier genomes themselves ---------
    want_sp = collections.Counter()
    with open(os.path.join(a.out, "pilot_sample.tsv")) as f:
        f.readline()
        for ln in f:
            want_sp[ln.split("\t")[5]] += 1
    src = {}
    with open(os.path.join(a.db, "source_manifest.jsonl")) as f:
        for ln in f:
            j = json.loads(ln)
            src.setdefault(j["fna_path"], j["assembly_accession_version"])
    # one genome per species, round robin, until --genomes picked
    by_sp = collections.defaultdict(list)
    for p in sorted(src):
        by_sp[p.split("/")[0].replace("genomes_", "").replace("_done_", "")].append(p)
    sps = sorted(by_sp, key=lambda s: (-want_sp.get(s, 0), s))
    chosen, i = [], 0
    while len(chosen) < a.genomes and i < 50:
        for s in sps:
            if len(chosen) >= a.genomes:
                break
            if i < len(by_sp[s]):
                chosen.append(by_sp[s][i])
        i += 1
    print(f"control genomes {len(chosen)} from {len(sps)} species groups", flush=True)

    pos, truth = [], []
    for p in chosen:
        full = os.path.join(a.wd, p)
        if not os.path.exists(full):
            print(f"  skip missing {p}", flush=True)
            continue
        recs = [(r, s) for r, s in fasta_stream(full) if len(s) > 200000]
        if not recs:
            continue
        rid, contig = max(recs, key=lambda t: len(t[1]))
        gf = pyrodigal.GeneFinder(closed=True)
        try:
            gf.train(contig)
        except Exception as e:
            print(f"  train failed {p}: {e}", flush=True)
            continue
        genes = gf.find_genes(contig)
        comp = [g for g in genes
                if not g.partial_begin and not g.partial_end and g.confidence() >= 99.0]
        shortg = [g for g in comp if 150 <= (g.end - g.begin + 1) <= 300]
        midg = [g for g in comp if 600 <= (g.end - g.begin + 1) <= 1500]
        rng.shuffle(shortg)
        rng.shuffle(midg)
        acc = src[p]
        n_short = max(1, a.per_genome // 4)
        n_mid = a.per_genome - n_short
        for tag, gl, n in (("short", shortg, n_short), ("mid", midg, n_mid)):
            for j, g in enumerate(gl[:n]):
                dna = contig[g.begin - 1:g.end]
                if g.strand < 0:
                    dna = revcomp(dna)
                prot = g.translate(include_stop=False)
                base = f"{acc}_{rid}_{g.begin-1}_{g.end}_{'+' if g.strand>0 else '-'}"
                if tag == "short":
                    cid = f"CTRL|pos_short_{base}"
                    pos.append((cid, dna, "pos_short"))
                    truth.append((cid, "pos_short", acc, rid, g.begin - 1, g.end,
                                  "+" if g.strand > 0 else "-", 0, 0, len(dna), prot))
                else:
                    cid = f"CTRL|pos_full_{base}"
                    pos.append((cid, dna, "pos_full"))
                    truth.append((cid, "pos_full", acc, rid, g.begin - 1, g.end,
                                  "+" if g.strand > 0 else "-", 0, 0, len(dna), prot))
                    if len(dna) > 101 + 90:
                        d5 = dna[101:]
                        cid5 = f"CTRL|pos_trunc5_{base}"
                        pos.append((cid5, d5, "pos_trunc5"))
                        truth.append((cid5, "pos_trunc5", acc, rid, g.begin - 1, g.end,
                                      "+" if g.strand > 0 else "-", 101, 0, len(d5), prot))
                    if len(dna) > 101 + 70 + 90:
                        db = dna[101:len(dna) - 70]
                        cidb = f"CTRL|pos_truncb_{base}"
                        pos.append((cidb, db, "pos_truncb"))
                        truth.append((cidb, "pos_truncb", acc, rid, g.begin - 1, g.end,
                                      "+" if g.strand > 0 else "-", 101, 70, len(db), prot))
        print(f"  {p} {rid}: complete {len(comp)}  short {len(shortg)}  mid {len(midg)}",
              flush=True)

    with open(os.path.join(a.out, "controls.fna"), "w") as f:
        for cid, seq, kind, _src in neg:
            f.write(f">{cid} kind={kind}\n{seq}\n")
        for cid, seq, kind in pos:
            f.write(f">{cid} kind={kind}\n{seq}\n")
    with open(os.path.join(a.out, "controls_truth.tsv"), "w") as f:
        f.write("control_id\tkind\tassembly\tcontig\tgene_start0\tgene_end0\tgene_strand\t"
                "trim5_nt\ttrim3_nt\textract_len_nt\ttrue_protein\n")
        for r in truth:
            f.write("\t".join(str(x) for x in r) + "\n")
        for cid, seq, kind, k in neg:
            f.write(f"{cid}\t{kind}\t\t\t\t\t\t0\t0\t{len(seq)}\t\n")

    meta = dict(negatives=len(neg), positives=len(pos),
                positive_kinds=dict(collections.Counter(k for _, _, k in pos)),
                genomes_used=len(chosen), seed=a.seed,
                seconds=round(time.time() - t0, 1))
    with open(os.path.join(a.out, "controls.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

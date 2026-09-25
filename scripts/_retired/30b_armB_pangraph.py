#!/usr/bin/env python3
"""ARM B, PanGraph engine -- same schema as 30_armB_graph.py, different graph.

PanGraph is the better conceptual fit for this problem: it groups homologous
sequence by nucleotide alignment ALONE, with no gene annotation, so every genome
is a path through blocks and an insertion is a path difference. The junction we
want is literally a place where paths diverge between two shared core blocks:

    core block A                    core block B
         │                               │
    g1:  A ─────────────────────────────► B      empty allele   = A|B
    g2:  A ─────► X ────────────────────► B      filled allele  = A|X|B
    g3:  A ─────► X ────────────────────► B
    g4:  A ─────────────────────────────► B

So: for every ADJACENT pair of core blocks, collect what each genome puts
between them. Genomes with nothing -> empty carriers. Genomes with block(s)
-> filled carriers. That is the whole algorithm.

Use this engine when the cluster is small enough (<= a few hundred genomes) and
you want annotation-independent block homology rather than minigraph's
reference-anchored bubbles. Run both and intersect for the highest confidence.

  pangraph build -o graph.json <fastas>
  pangraph export ...             (JSON is parsed directly here)
"""
import argparse
import json
import os
import statistics
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log, run                                   # noqa: E402


def build(paths, out_json, pangraph, threads, extra=None):
    cmd = [pangraph, "build", "-o", out_json, "-j", str(threads)]
    if extra:
        cmd += extra.split()
    cmd += paths
    run(cmd)
    return out_json


def load_paths(js):
    """Return {genome: [(block_id, strand, len), ...]} and {block_id: len}.

    PanGraph's JSON key names have shifted across versions, so each field is
    looked up through a few aliases and a clear error is raised rather than a
    silent empty result.
    """
    with open(js) as fh:
        g = json.load(fh)

    def pick(d, *names):
        for n in names:
            if n in d:
                return d[n]
        return None

    blocks = pick(g, "blocks", "Blocks") or []
    blen = {}
    if isinstance(blocks, dict):
        blocks = list(blocks.values())
    for b in blocks:
        bid = str(pick(b, "id", "ID", "name"))
        seq = pick(b, "sequence", "consensus", "Sequence") or ""
        blen[bid] = pick(b, "length", "len") or len(seq)

    raw = pick(g, "paths", "Paths") or []
    if isinstance(raw, dict):
        raw = [dict(v, name=k) for k, v in raw.items()]
    out = {}
    for p in raw:
        name = str(pick(p, "name", "Name", "id"))
        nodes = pick(p, "blocks", "nodes", "Blocks") or []
        seq = []
        for n in nodes:
            if isinstance(n, dict):
                bid = str(pick(n, "id", "block_id", "ID", "name"))
                strand = pick(n, "strand", "orientation")
                strand = "+" if strand in (True, "+", 1, None) else "-"
            else:
                bid, strand = str(n), "+"
            seq.append((bid, strand, blen.get(bid, 0)))
        out[name] = seq
    if not out:
        raise SystemExit("could not parse paths from %s -- run `pangraph export` "
                         "and check the key names against load_paths()" % js)
    return out, blen


def core_blocks(paths, min_frac, max_copies=1):
    """Single-copy blocks present in >= min_frac of genomes: the anchors."""
    n = len(paths)
    per = defaultdict(int)
    dup = set()
    for gname, seq in paths.items():
        c = Counter(b for b, _, _ in seq)
        for b, k in c.items():
            per[b] += 1
            if k > max_copies:
                dup.add(b)
    return {b for b, k in per.items() if k >= min_frac * n and b not in dup}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cluster-id", default=None)
    ap.add_argument("--pangraph", default="pangraph")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--build-args", default=None,
                    help="extra flags passed to `pangraph build`")
    ap.add_argument("--core-frac", type=float, default=0.90,
                    help="block must be single-copy in >= this fraction of genomes")
    ap.add_argument("--min-insert", type=int, default=300)
    ap.add_argument("--max-insert", type=int, default=5000)
    ap.add_argument("--min-empty", type=int, default=1)
    ap.add_argument("--min-filled", type=int, default=1)
    ap.add_argument("--reuse", action="store_true")
    args = ap.parse_args()

    cid = args.cluster_id or os.path.basename(args.manifest).split(".")[0]
    os.makedirs(args.outdir, exist_ok=True)
    paths_in = [l.strip() for l in open(args.manifest) if l.strip()]
    js = os.path.join(args.outdir, "pangraph.json")
    if not (args.reuse and os.path.exists(js)):
        log("%s: pangraph build over %d genomes" % (cid, len(paths_in)))
        build(paths_in, js, args.pangraph, args.threads, args.build_args)

    gpaths, blen = load_paths(js)
    log("%d genome paths, %d blocks" % (len(gpaths), len(blen)))
    core = core_blocks(gpaths, args.core_frac)
    log("%d core anchor blocks (>=%.0f%% single-copy)" % (len(core), 100 * args.core_frac))

    # for each genome, the reduced path over core blocks + what sits between
    between = defaultdict(dict)          # (A,B) -> {genome: [inner blocks]}
    for gname, seq in gpaths.items():
        anchors = [(i, b, s) for i, (b, s, _) in enumerate(seq) if b in core]
        for (i1, b1, s1), (i2, b2, s2) in zip(anchors, anchors[1:]):
            inner = seq[i1 + 1:i2]
            # orient the junction so A|B and B|A are the same junction
            key = (b1, b2) if b1 <= b2 else (b2, b1)
            between[key][gname] = inner

    ev = open(os.path.join(args.outdir, "armB_events.tsv"), "w")
    pr = open(os.path.join(args.outdir, "armB_pairs.tsv"), "w")
    ev.write("\t".join([
        "event_id", "cluster_id", "chrom", "bubble_start", "bubble_end",
        "is_inversion", "len_empty", "len_filled", "insert_size",
        "n_empty", "n_filled", "n_missing", "allele_state",
        "empty_genomes", "filled_genomes", "qc_tier", "flags"]) + "\n")
    pr.write("event_id\tcluster_id\tfilled_genome\tempty_genome\tchrom\t"
             "bubble_start\tbubble_end\tinsert_size\n")

    n_ev = 0
    for (a, b), per_g in sorted(between.items()):
        lens = {g: sum(blen.get(x[0], 0) for x in inner) for g, inner in per_g.items()}
        if len(lens) < 2:
            continue
        shortest = min(lens.values())
        empty = [g for g, L in lens.items() if L <= shortest + 50]
        filled = [g for g, L in lens.items() if L >= shortest + args.min_insert]
        if len(empty) < args.min_empty or len(filled) < args.min_filled:
            continue
        insert = int(statistics.median([lens[g] for g in filled]) -
                     statistics.median([lens[g] for g in empty]))
        if not (args.min_insert <= insert <= args.max_insert):
            continue
        eid = "%s.P%06d" % (cid, n_ev)
        n_missing = len(gpaths) - len(lens)
        ev.write("\t".join(map(str, [
            eid, cid, "%s|%s" % (a, b), 0, 0, 0,
            shortest, int(statistics.median([lens[g] for g in filled])), insert,
            len(empty), len(filled), n_missing, "PRESENCE_ABSENCE",
            ",".join(sorted(empty)[:20]), ",".join(sorted(filled)[:20]),
            "graph_only", "pangraph_core_junction"])) + "\n")
        for f in sorted(filled)[:1]:
            for e in sorted(empty)[:1]:
                pr.write("%s\t%s\t%s\t%s\t%s\t0\t0\t%d\n"
                         % (eid, cid, f, e, "%s|%s" % (a, b), insert))
        n_ev += 1

    ev.close()
    pr.close()
    log("%s: %d core-junction empty/filled alleles" % (cid, n_ev))
    log("coordinates are block-relative -- stage 31 supplies nucleotide positions")


if __name__ == "__main__":
    main()

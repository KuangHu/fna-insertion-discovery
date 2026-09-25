#!/usr/bin/env python3
"""Positional novelty test: which discovered inserts overlap a known IS call?

Replaces an earlier size-compatibility heuristic that was wrong twice over. It
asked only whether SOME ISEScan call in a filled carrier had a length within 20%
of the insert, and it `break`s at the first such call -- so a locus was credited
to whichever family the loop happened to reach first, and two unrelated 1.3 kb
elements at opposite ends of a chromosome counted as a match. Its family tally
put zero IS110 in the matches even though the pipeline's 1279 bp mode carries
DEDD_Tnp_IS110 at E=8e-28.

This version places the insert. Each inserted sequence is mapped back to a
genome that carries it, and the resulting interval is intersected with that
genome's ISEScan calls. An insert is `known` only when a call actually overlaps
its footprint by --min-overlap of the shorter interval.

The `novel` fraction is an UPPER BOUND on novelty, not a novelty rate:
ISEScan itself misses diverged and atypical elements, so an unmatched insert may
be a real element ISEScan cannot see, an assembly artefact, or a non-IS mobile
element. It bounds the pool worth annotating; it does not claim discovery.

    87_novelty_positional.py --inserts armB_inserts.fna --events armB_events.tsv \\
        --isescan-dir isescan_run/tsv --genome-dir isescan_run/fna --out novelty
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import paf as pafmod                                    # noqa: E402
from lib.util import log, overlap                                # noqa: E402


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inserts", nargs="+", required=True)
    ap.add_argument("--events", nargs="+", required=True)
    ap.add_argument("--isescan-dir", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.5)
    ap.add_argument("--min-ident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    ins = {}
    for p in args.inserts:
        ins.update(read_fasta(p))
    carriers = {}
    for p in args.events:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = [x for x in r.get("filled_genomes", "").split(",") if x]
                if g:
                    carriers[r["event_id"]] = g
    log("%d inserts, %d events with a filled carrier" % (len(ins), len(carriers)))

    isc = defaultdict(list)
    for p in glob.glob(os.path.join(args.isescan_dir, "*.tsv")):
        s = os.path.basename(p)[:-4]
        with open(p) as fh:
            # ISEScan writes COMMA-separated content into a .tsv filename.
            # Assuming tab silently parses each line into one bogus field, every
            # lookup KeyErrors, and the run reports a near-empty IS catalogue --
            # which inflates apparent novelty instead of failing. Sniff, as
            # 60_benchmark_armA_vs_isescan.py already does.
            head = fh.readline()
            delim = "," if head.count(",") > head.count("\t") else "\t"
            fh.seek(0)
            for r in csv.DictReader(fh, delimiter=delim):
                try:
                    a, b = int(r["isBegin"]), int(r["isEnd"])
                except (KeyError, ValueError, TypeError):
                    continue
                isc[(s, r.get("seqID", "."))].append(
                    (min(a, b), max(a, b), r.get("family", "."),
                     r.get("cluster", ".")))
    log("%d ISEScan intervals over %d contigs"
        % (sum(len(v) for v in isc.values()), len(isc)))

    # group inserts by the carrier genome they will be mapped into, so each
    # genome is indexed once rather than once per insert
    by_genome = defaultdict(list)
    for eid, seq in ins.items():
        for g in carriers.get(eid, [])[:1]:
            by_genome[g].append((eid, seq))

    tmp = tempfile.mkdtemp(prefix="novel.")
    rows, stat = [], Counter()
    for g, items in sorted(by_genome.items()):
        gp = os.path.join(args.genome_dir, g + ".fna")
        if not os.path.exists(gp):
            stat["genome_missing"] += len(items)
            continue
        q = os.path.join(tmp, g + ".fna")
        with open(q, "w") as fh:
            for eid, seq in items:
                fh.write(">%s\n%s\n" % (eid, seq))
        paf = os.path.join(tmp, g + ".paf")
        with open(paf, "w") as fh:
            subprocess.run([args.minimap2, "-c", "-x", "asm10", "-N", "50",
                            "-p", "0.2", "--secondary=yes", "-t",
                            str(args.threads), gp, q],
                           stdout=fh, stderr=subprocess.DEVNULL, check=False)
        placed = defaultdict(list)
        for r in pafmod.parse(paf):
            if r.identity < args.min_ident or r.qspan < args.min_cov * r.qlen:
                continue
            placed[r.qname].append((r.tname, r.ts, r.te))
        for eid, seq in items:
            hits = placed.get(eid, [])
            if not hits:
                stat["insert_not_placed"] += 1
                rows.append([eid, g, 0, "NOT_PLACED", ".", ".", 0.0])
                continue
            best = None
            for tname, ts, te in hits:
                for a, b, fam, clu in isc.get((g, tname), ()):
                    ov = overlap(ts, te, a, b)
                    if ov <= 0:
                        continue
                    frac = ov / max(1, min(te - ts, b - a))
                    if best is None or frac > best[0]:
                        best = (frac, fam, clu)
            if best and best[0] >= args.min_overlap:
                stat["known"] += 1
                rows.append([eid, g, len(hits), "KNOWN", best[1], best[2],
                             round(best[0], 3)])
            else:
                stat["unmatched"] += 1
                rows.append([eid, g, len(hits), "UNMATCHED",
                             best[1] if best else ".", best[2] if best else ".",
                             round(best[0], 3) if best else 0.0])

    with open(args.out + "_per_insert.tsv", "w") as fh:
        fh.write("event_id\tcarrier\tn_placements\tclass\tisescan_family\t"
                 "isescan_cluster\toverlap_fraction\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    scored = stat["known"] + stat["unmatched"]
    fam = Counter(r[4] for r in rows if r[3] == "KNOWN")
    log("=" * 72)
    log("POSITIONAL NOVELTY (overlap >= %.0f%% of the shorter interval)"
        % (100 * args.min_overlap))
    log("  inserts placed and scored     %4d" % scored)
    log("  overlapping a known IS call   %4d  %5.1f%%"
        % (stat["known"], 100.0 * stat["known"] / scored if scored else 0))
    log("  no overlapping call           %4d  %5.1f%%   <- upper bound on novelty"
        % (stat["unmatched"], 100.0 * stat["unmatched"] / scored if scored else 0))
    if stat["insert_not_placed"]:
        log("  insert would not place        %4d  (excluded from the denominator)"
            % stat["insert_not_placed"])
    log("")
    log("  families among the matches:")
    for k, v in fam.most_common(10):
        log("    %-12s %4d" % (k, v))
    log("=" * 72)
    log("wrote %s_per_insert.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

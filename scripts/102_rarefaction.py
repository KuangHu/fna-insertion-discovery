#!/usr/bin/env python3
"""Does adding genomes still pay? Retrospective rarefaction on panels already run.

No download, no re-run: subsample the panels that exist and watch what the
curves do. The criterion is pre-registered in docs/SIMPLIFICATION.md §14 --
marginal yield in the 300->400 band against the 100->200 band, >70% keep going,
<40% stop.

THREE CURVES, NOT ONE. `event = (insert_md5, target_context_md5)`. A new genome
brings new flanking context, which can split one event in two purely because the
context string differs -- so part of any rise is KEY FRAGMENTATION, not
discovery. Separating them is the point:

  distinct insert_key   -> is the ELEMENT repertoire complete?
  distinct context_key  -> is the TARGET-SITE repertoire complete?
  distinct event        -> the number currently reported

If insert saturates while context keeps climbing, the reading is that the same
elements keep appearing at new target sites -- which is target-site preference,
measured with no annotation, and worth more than the headline count.

SECOND TEST -- are the unused genomes new territory or clones? For every genome
not in any panel, its maximum ANI to the used set, streamed from the raw skani
shards in one pass. Mostly >=99.9 means clonal near-relatives and more genomes
buy redundancy; a substantial share <99.5 means unsampled lineages exist.

    102_rarefaction.py --loci-to-event dedup_loci_to_event.tsv \\
        --used armA_genomes.txt --all-genomes all.txt \\
        --shards 'shards/ani_*.tsv' --species ecoli --out rarefaction/ecoli
"""
import argparse
import collections
import csv
import glob
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci-to-event", required=True)
    ap.add_argument("--used", default=None)
    ap.add_argument("--all-genomes", default=None)
    ap.add_argument("--shard-dir", default=None,
                    help="directory of skani shards; globbed as ani_*.tsv "
                         "INSIDE this script. Preferred over --shards: passing "
                         "a pattern through sbatch --export requires escaping "
                         "that Python's glob cannot then match (a literal "
                         "backslash), which silently read 0 shards.")
    ap.add_argument("--shards", nargs="+", default=None,
                    help="skani shard files. Accepts either a quoted glob or an "
                         "already-expanded list -- the submitting shell expands "
                         "the pattern before sbatch sees it, so a single-value "
                         "flag received 40 positional args and argparse exited 2.")
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", type=int, nargs="+",
                    default=[25, 50, 100, 200, 300, 400])
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    random.seed(args.seed)
    by_panel = collections.defaultdict(list)
    for r in csv.DictReader(open(args.loci_to_event), delimiter="\t"):
        by_panel[r["panel"]].append((r["insert_key"], r["context_key"]))
    panels = sorted(by_panel)
    log("%d panels, %d keyed loci" % (len(panels), sum(len(v) for v in by_panel.values())))
    if len(panels) < 2:
        log("FATAL: %d distinct panel labels -- the panel field is mis-derived "
            "and every subsample would be identical" % len(panels))
        return 2

    steps = [s for s in args.steps if s <= len(panels)]
    if len(panels) not in steps:
        steps.append(len(panels))
    curves = {}
    for n in steps:
        ins, ctx, ev = [], [], []
        for _ in range(args.reps if n < len(panels) else 1):
            pick = random.sample(panels, n)
            I, C, E = set(), set(), set()
            for p in pick:
                for a, b in by_panel[p]:
                    I.add(a)
                    C.add(b)
                    E.add((a, b))
            ins.append(len(I))
            ctx.append(len(C))
            ev.append(len(E))
        curves[n] = (statistics.mean(ins), statistics.mean(ctx), statistics.mean(ev))
        log("  n=%-4d  insert %8.0f   context %8.0f   event %8.0f"
            % (n, curves[n][0], curves[n][1], curves[n][2]))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_rarefaction.tsv", "w") as fh:
        fh.write("species\tn_panels\tdistinct_insert\tdistinct_context\tdistinct_event\n")
        for n in steps:
            fh.write("%s\t%d\t%.1f\t%.1f\t%.1f\n" % (args.species, n, *curves[n]))

    top = max(curves)
    prev = max((x for x in curves if x < top), default=None)

    def marg(lo, hi, idx):
        if lo is None or hi is None or lo not in curves or hi not in curves:
            return float("nan")
        return (curves[hi][idx] - curves[lo][idx]) / (hi - lo)

    log("=" * 76)
    log("RAREFACTION -- %s   (criterion fixed in SIMPLIFICATION.md §14)" % args.species)
    log("")
    log("  %-16s %12s %12s %10s"
        % ("curve", "100->200/pan", "%d->%d/pan" % (prev or 0, top), "ratio"))
    log("  " + "-" * 54)
    verdicts = {}
    for idx, name in ((0, "insert_key"), (1, "context_key"), (2, "event")):
        a, b = marg(100, 200, idx), marg(prev, top, idx)
        rt = b / a if a else float("nan")
        v = ("KEEP GOING (>70%)" if rt > 0.70 else
             "SATURATING (<40%)" if rt < 0.40 else "MARGINAL (40-70%)")
        verdicts[name] = (rt, v)
        log("  %-16s %12.1f %12.1f %9.1f%%  %s" % (name, a, b, 100 * rt, v))
    log("")
    ri, rc_ = verdicts["insert_key"][0], verdicts["context_key"][0]
    if ri == ri and rc_ == rc_ and ri < 0.40 <= rc_:
        log("  *** ELEMENT REPERTOIRE SATURATED, TARGET SITES STILL CLIMBING ***")
        log("  The same elements keep appearing at new target sites. That is")
        log("  target-site preference measured with no annotation, and it is a")
        log("  stronger result than any total.")
    log("=" * 76)

    # ---- unused genomes: new territory or clones? -------------------------
    shard_files = []
    if args.shard_dir:
        shard_files = sorted(glob.glob(os.path.join(args.shard_dir, "ani_*.tsv")))
    elif args.shards:
        for pat in args.shards:
            shard_files.extend(sorted(glob.glob(pat))
                               if any(c in pat for c in "*?[") else [pat])
        shard_files = sorted(set(shard_files))
    if args.used and args.all_genomes and (args.shard_dir or args.shards):
        # ASSERTION: zero shards is a path/escaping failure, not an empty matrix.
        # Every unused genome then defaults to max-ANI 0.0 and lands in the
        # "<99.0" band, producing a clean 100%/0%/0%/0% split that reads as a
        # real finding.
        if not shard_files:
            log("FATAL: 0 shard files matched. Use --shard-dir; a glob passed "
                "through sbatch --export cannot survive to glob.glob().")
            return 2
        used = {os.path.basename(l.strip())[:-4] for l in open(args.used) if l.strip()}
        allg = {os.path.basename(l.strip())[:-4] for l in open(args.all_genomes) if l.strip()}
        # matched/total for the GENOME-NAME JOIN. If `used` and `all_genomes`
        # are basenamed differently, `used - allg` is non-empty and `unused`
        # silently becomes the whole pool, making every band meaningless.
        stray = used - allg
        log("  genome-name join: %d/%d used genomes present in the full list "
            "(%.1f%%)" % (len(used) - len(stray), len(used),
                          100.0 * (len(used) - len(stray)) / len(used) if used else 0))
        if stray:
            log("FATAL: %d used genomes are absent from --all-genomes, e.g. %s. "
                "The two lists are named differently, so the unused set is wrong."
                % (len(stray), sorted(stray)[:3]))
            return 2
        unused = allg - used
        log("")
        log("  used %d, total %d, UNUSED %d" % (len(used), len(allg), len(unused)))
        best = collections.defaultdict(float)
        files = shard_files
        for f in files:
            for i, line in enumerate(open(f)):
                if i == 0:
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                a = os.path.basename(p[0])[:-4]
                b = os.path.basename(p[1])[:-4]
                if a == b:
                    continue
                try:
                    v = float(p[2])
                except ValueError:
                    continue
                if a in unused and b in used and v > best[a]:
                    best[a] = v
                if b in unused and a in used and v > best[b]:
                    best[b] = v
        log("  shards read: %d ; unused genomes with any measured pair: %d/%d"
            % (len(files), len(best), len(unused)))
        h = collections.Counter()
        for g in unused:
            v = best.get(g, 0.0)
            h[">=99.9" if v >= 99.9 else "99.5-99.9" if v >= 99.5 else
              "99.0-99.5" if v >= 99.0 else "<99.0 or unreported"] += 1
        log("")
        log("  MAX ANI OF EACH UNUSED GENOME TO THE USED SET")
        for k in (">=99.9", "99.5-99.9", "99.0-99.5", "<99.0 or unreported"):
            log("    %-22s %6d   %5.1f%%" % (k, h[k], 100.0 * h[k] / len(unused)))
        lo = h["99.0-99.5"] + h["<99.0 or unreported"]
        log("")
        log("    below 99.5: %d (%.1f%%) -- %s"
            % (lo, 100.0 * lo / len(unused),
               "unsampled lineages exist, more genomes is new ground"
               if lo > 0.25 * len(unused) else
               "mostly clonal near-relatives; more genomes buys redundancy"))
        with open(args.out + "_unused_ani.tsv", "w") as fh:
            fh.write("genome\tmax_ani_to_used\n")
            for g in sorted(unused):
                fh.write("%s\t%.3f\n" % (g, best.get(g, 0.0)))
    log("wrote %s_rarefaction.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

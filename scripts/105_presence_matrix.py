#!/usr/bin/env python3
"""LINE M3 -- how many genomes of this species carry this insert? (presence/absence)

M1 asks "is it >=3 copies inside ONE genome". M2 asks "does it occur at another
event". Neither can answer "how widespread is this sequence in the species" --
and a single-copy element can be present in most genomes. That gap is what this
fills.

The gate is identical to M1's (identity >= 0.95, query coverage >= 0.90) so the
two lines are directly comparable rather than measuring different things with
different tolerances.

WHY IT MIGHT BE REDUNDANT, stated up front: if M1-positive inserts are also
M3-high and M1-negative inserts are M3-low, M3 adds nothing. **The number that
justifies this module is how many M1-NEGATIVE inserts are M3-HIGH** -- widely
distributed single-copy elements that the mobility line cannot see. That is
reported first.

M3 is emitted as a FIELD alongside M1/M2. It does NOT change the definition of
`usable`, S1, or M1's gate.

    105_presence_matrix.py --inserts <sp>_inserts.fna --genomes <dir> \\
        --mobility mobility_mobility.tsv --events <sp>_events.tsv --out x
"""
import argparse, collections, csv, glob, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--mobility", default=None)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--par", type=int, default=8)
    args = ap.parse_args()

    genomes = sorted(glob.glob(os.path.join(args.genomes, "*.fna")))
    n_ins = sum(1 for l in open(args.inserts) if l.startswith(">"))
    log("%d inserts, %d genomes" % (n_ins, len(genomes)))
    if not genomes or not n_ins:
        log("FATAL: empty input")
        return 2

    tmp = tempfile.mkdtemp(prefix="m3.")
    present = collections.Counter()
    sh = os.path.join(tmp, "one.sh")
    with open(sh, "w") as fh:
        fh.write('#!/bin/bash\nminimap2 -c -x asm10 -t %d "$1" "%s" 2>/dev/null | '
                 'awk -F"\\t" \'$11>0 && $10/$11>=%f && $11/$2>=%f {print $1}\' | '
                 'sort -u > "$2"\n' % (args.threads, args.inserts,
                                       args.min_ident, args.min_cov))
    os.chmod(sh, 0o755)
    jobs = [(g, os.path.join(tmp, os.path.basename(g) + ".hit")) for g in genomes]
    with open(os.path.join(tmp, "list"), "w") as fh:
        for g, o in jobs:
            fh.write("%s\t%s\n" % (g, o))
    subprocess.run("awk -F'\\t' '{print $1\"\\n\"$2}' %s | xargs -P %d -n 2 %s"
                   % (os.path.join(tmp, "list"), args.par, sh),
                   shell=True, check=False)
    done = 0
    for g, o in jobs:
        if not os.path.exists(o):
            continue
        done += 1
        for line in open(o):
            # The FASTA record is ">E000001|<species> len=... md5=..." so the
            # minimap2 query name carries the species suffix, while the events
            # table id does not. Strip it. EIGHTH wrong-key failure in this
            # project: the earlier "93.9% of inserts detected" check PASSED
            # because it compared present-vs-FASTA, and the break was
            # present-vs-EVENTS -- a different join at a different level.
            present[line.strip().split("|")[0]] += 1
    # standing practice 6: matched/total before any rate derived from the join
    log("genome scan: %d/%d genomes produced output (%.1f%%)"
        % (done, len(genomes), 100.0 * done / len(genomes)))
    if done < 0.9 * len(genomes):
        log("FATAL: %d of %d genomes produced no hit file" % (len(genomes) - done, len(genomes)))
        return 2
    hit_any = len(present)
    log("inserts detected in >=1 genome: %d/%d (%.1f%%)"
        % (hit_any, n_ins, 100.0 * hit_any / n_ins))
    if hit_any == 0:
        log("FATAL: no insert matched any genome -- a name or format mismatch")
        return 2

    ev = {r["event_id"]: r for r in csv.DictReader(open(args.events), delimiter="\t")}
    # ASSERTION on the join that actually feeds the output: present keys must
    # resolve to event ids. Checking "did the inserts hit any genome" is a
    # DIFFERENT join and passed while this one was broken.
    joined = sum(1 for k in present if k in ev)
    log("present->events join: %d/%d detected inserts resolve to an event id "
        "(%.1f%%)" % (joined, len(present), 100.0 * joined / len(present) if present else 0))
    if present and joined < 0.5 * len(present):
        log("FATAL: fewer than half of the detected inserts match an event id -- "
            "a key-format mismatch between the FASTA record names and the "
            "events table")
        return 2
    rows = []
    for eid, r in ev.items():
        n = present.get(eid, 0)
        rows.append({"event_id": eid, "species": args.species,
                     "insert_md5": r.get("insert_md5", "."),
                     "inserted_len": r.get("inserted_len", "."),
                     "M1_within_genome": r.get("M1_within_genome", "."),
                     "M2_cross_event": r.get("M2_cross_event", "."),
                     "M3_n_genomes_present": n,
                     "M3_frac_genomes_present": "%.4f" % (n / done if done else 0)})
    cols = list(rows[0].keys())
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_presence.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    def band(fr):
        return ("0" if fr == 0 else "<1%" if fr < 0.01 else "1-10%" if fr < 0.10
                else "10-50%" if fr < 0.50 else ">=50%")
    log("=" * 74)
    log("LINE M3 -- %s   (%d genomes scanned)" % (args.species, done))
    log("")
    log("  THE NUMBER THAT JUSTIFIES THIS MODULE:")
    neg = [r for r in rows if r["M1_within_genome"] != "True"]
    negwide = [r for r in neg if float(r["M3_frac_genomes_present"]) >= 0.10]
    log("    M1-NEGATIVE but present in >=10%% of genomes: %d of %d M1-negatives"
        % (len(negwide), len(neg)))
    log("    = %.1f%% of all events" % (100.0 * len(negwide) / len(rows)))
    log("    If small, M1 already covers the ground and M3 is redundant.")
    log("")
    log("  M3 distribution, stratified by M1:")
    log("  %-16s %10s %10s" % ("present in", "M1 pos", "M1 neg"))
    log("  " + "-" * 38)
    for b in ("0", "<1%", "1-10%", "10-50%", ">=50%"):
        p = sum(1 for r in rows if r["M1_within_genome"] == "True"
                and band(float(r["M3_frac_genomes_present"])) == b)
        q = sum(1 for r in rows if r["M1_within_genome"] != "True"
                and band(float(r["M3_frac_genomes_present"])) == b)
        log("  %-16s %10d %10d" % (b, p, q))
    log("=" * 74)
    log("wrote %s_presence.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

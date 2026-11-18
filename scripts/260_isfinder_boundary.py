#!/usr/bin/env python3
"""Boundary test against curated IS elements (ISfinder). A TEST MODULE: it reads
the catalogue and changes nothing in it.

Question: where a catalogue insert IS a known, curated IS element, did the
pipeline cut it at the element's ends?

    260_isfinder_boundary.py --isfinder IS.fna --isfinder-sha256 <hex> \\
        --events 'catalogue_v3/*_events.tsv' --panel config/isfinder_panel.tsv \\
        --out <new dir> [--threads 8]

Method. Every event's `insert_seq` is aligned to every ISfinder sequence
(minimap2 -c --eqx). For each insert the best IS (most matching bases) is kept,
and the IS ends are PROJECTED into the insert: IS bases the aligner left
unaligned at either end are assumed colinear, so a clipped end still gives a
position. Then

    left_extra  = insert bases before the projected IS start
    right_extra = insert bases after the projected IS end
    extra       = left_extra + right_extra
    delta       = extra - junction_overlap_bp

The catalogue's insert is the segment between the two flanks of the empty
allele; where the junction is a direct repeat (the TSD) the segment holds the
element plus ONE copy of the repeat, at either end or split between them. The
catalogue reports the repeat length itself (`junction_overlap_bp`). So a
base-exact call of a complete element has both extras >= 0 and delta == 0.

Classes (only inserts whose best hit has identity >= --min-ident and covers
>= --min-cov of the IS are judged; the rest are `no_is`):

    exact          delta == 0, both extras >= 0
    exact_shifted  -near <= delta < 0 and neither extra below delta: the
                   junction repeat runs on by chance into the element (its
                   first base equals the base after the repeat, p ~ 1/4 per
                   end), so `junction_overlap_bp` is longer than the repeat
                   and the cut is the same allele written one of several
                   equivalent ways. Not a boundary error.
    near           |delta| <= --near bp, otherwise
    over           delta >  --near : the insert holds more than the IS and its repeat
    under          delta < -near   : the insert ends inside the IS

`over` and `under` are NOT by themselves errors: an insert may be an IS plus
cargo, two elements, or a genuinely truncated copy. `exact` and
`exact_shifted` are positive evidence of a correct boundary; the others are
the cases to inspect. The ISfinder record itself can be off by a base or two
at an end (see docs/ISFINDER_BOUNDARY.md), which shows up as `near`.
"""
import argparse
import collections
import csv
import glob
import hashlib
import json
import os
import subprocess
import sys
import time

csv.field_size_limit(1 << 30)


CLASSES = ("exact", "exact_shifted", "near", "over", "under")


class Fail(Exception):
    pass


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def place(qlen, qs, qe, strand, tlen, ts, te):
    """Project the IS ends into the insert. Returns (left_extra, right_extra):
    insert bases outside the projected element, in the insert's orientation.
    Negative means the insert ends inside the element on that side."""
    head, tail = ts, tlen - te          # IS bases unaligned at its start / end
    if strand == "+":
        return qs - head, (qlen - qe) - tail
    if strand == "-":
        return qs - tail, (qlen - qe) - head
    raise Fail("bad strand %r" % strand)


def classify(left, right, overlap, near):
    delta = left + right - overlap
    if delta == 0 and left >= 0 and right >= 0:
        return "exact", delta
    if -near <= delta < 0 and min(left, right) >= delta:
        return "exact_shifted", delta
    if abs(delta) <= near:
        return "near", delta
    return ("over" if delta > 0 else "under"), delta


def read_events(patterns, fasta_out):
    """Write the inserts as FASTA; return {qid: row summary}."""
    need = ["event_id", "species", "insert_seq", "inserted_len",
            "junction_overlap_bp", "insert_frame_status",
            "decomposition_method", "S1_structurally_clean"]
    files = sorted(f for p in patterns for f in glob.glob(p))
    if not files:
        raise Fail("no event tables match %s" % patterns)
    ev, inputs = {}, []
    with open(fasta_out, "w") as fa:
        for f in files:
            inputs.append({"path": os.path.realpath(f), "sha256": sha256_file(f)})
            with open(f, newline="") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                miss = [c for c in need if c not in (rd.fieldnames or [])]
                if miss:
                    raise Fail("%s lacks columns %s" % (f, miss))
                for r in rd:
                    seq = r["insert_seq"]
                    if not seq:
                        continue
                    if len(seq) != int(r["inserted_len"]):
                        raise Fail("%s %s: insert_seq length != inserted_len"
                                   % (f, r["event_id"]))
                    qid = "%s|%s" % (r["species"], r["event_id"])
                    if qid in ev:
                        raise Fail("duplicate event %s" % qid)
                    ev[qid] = {k: r[k] for k in need if k != "insert_seq"}
                    ev[qid]["junction_overlap_bp"] = int(r["junction_overlap_bp"] or 0)
                    fa.write(">%s\n%s\n" % (qid, seq))
    return ev, inputs


def read_is(path, fasta_out):
    """{name: length}. ISfinder headers are name_family_group. The public copy
    repeats ~400 records verbatim: an identical repeat is dropped, the same name
    with another sequence is fatal. The de-duplicated set is written out and is
    what gets aligned."""
    seqs, cur, dropped = {}, None, 0
    recs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = [line[1:].split()[0], []]
                recs.append(cur)
            elif cur is not None:
                cur[1].append(line.strip().upper())
    for name, parts in recs:
        seq = "".join(parts)
        if not seq:
            raise Fail("empty IS record %s" % name)
        if name in seqs:
            if seqs[name] != seq:
                raise Fail("IS name %s is used for two different sequences" % name)
            dropped += 1
            continue
        seqs[name] = seq
    with open(fasta_out, "w") as fa:
        for name, seq in seqs.items():
            fa.write(">%s\n%s\n" % (name, seq))
    return {n: len(q) for n, q in seqs.items()}, dropped


def is_parts(tname):
    p = tname.split("_")
    if len(p) >= 3:
        return "_".join(p[:-2]), p[-2], p[-1]
    return tname, "", ""


def read_panel(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            out[r["is_name"]] = r
    return out


def best_hits(paf, ev, is_len, min_ident, min_cov):
    """Per insert: the best judged hit, and how many OTHER elements it holds."""
    hits = collections.defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            q, qlen, qs, qe, strand, t, tlen, ts, te, nm, alen = (
                c[0], int(c[1]), int(c[2]), int(c[3]), c[4], c[5],
                int(c[6]), int(c[7]), int(c[8]), int(c[9]), int(c[10]))
            if q not in ev or is_len.get(t) != tlen:
                raise Fail("PAF names a sequence that was not given: %s / %s" % (q, t))
            ident, cov = nm / alen, (te - ts) / tlen
            if ident >= min_ident and cov >= min_cov:
                hits[q].append((nm, ident, cov, qlen, qs, qe, strand, t, tlen, ts, te))
    out = {}
    for q, hs in hits.items():
        hs.sort(key=lambda h: (-h[0], -h[1], h[7]))
        b = hs[0]
        others = set()
        for h in hs[1:]:
            ov = min(h[5], b[5]) - max(h[4], b[4])
            if ov < 0.5 * (h[5] - h[4]):          # a different stretch of the insert
                others.add((h[4] // 200, h[5] // 200))
        out[q] = (b, len(others))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--isfinder", required=True)
    ap.add_argument("--isfinder-sha256", required=True)
    ap.add_argument("--events", required=True, nargs="+")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--near", type=int, default=3)
    ap.add_argument("--min-n", type=int, default=10,
                    help="per-IS table lists elements with at least this many judged inserts")
    a = ap.parse_args()
    try:
        run(a)
    except Fail as e:
        sys.exit("FATAL: %s" % e)


def run(a):
    if os.path.lexists(a.out):
        raise Fail("%s exists; every run gets a new directory" % a.out)
    got = sha256_file(a.isfinder)
    if got != a.isfinder_sha256:
        raise Fail("ISfinder file has sha256 %s, expected %s" % (got, a.isfinder_sha256))
    here = os.path.dirname(os.path.abspath(__file__))
    commit = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", here, "status", "--porcelain"],
                                capture_output=True, text=True).stdout.strip())
    os.makedirs(a.out)
    fa, paf = os.path.join(a.out, "inserts.fna"), os.path.join(a.out, "hits.paf")
    ev, inputs = read_events(a.events, fa)
    ref = os.path.join(a.out, "isfinder.dedup.fna")
    is_len, dropped = read_is(a.isfinder, ref)
    panel = read_panel(a.panel)
    unknown = sorted(n for n in panel if not any(is_parts(t)[0] == n for t in is_len))
    if unknown:
        raise Fail("panel names not in the ISfinder file: %s" % unknown)

    cmd = ["minimap2", "-c", "--eqx", "-x", "asm20", "-N", "50", "-p", "0.3",
           "-t", str(a.threads), ref, fa]
    # no --end-bonus: it drags the alignment to the insert's end through the
    # repeat bases (as an indel), which hides exactly the extra being measured
    ver = subprocess.run(["minimap2", "--version"], capture_output=True, text=True).stdout.strip()
    with open(paf, "w") as out, open(os.path.join(a.out, "minimap2.log"), "w") as err:
        rc = subprocess.run(cmd, stdout=out, stderr=err).returncode
    if rc != 0:
        raise Fail("minimap2 exit %d (see minimap2.log)" % rc)

    best = best_hits(paf, ev, is_len, a.min_ident, a.min_cov)
    cols = ["species", "event_id", "inserted_len", "junction_overlap_bp",
            "insert_frame_status", "decomposition_method", "S1_structurally_clean",
            "is_name", "is_family", "is_group", "is_len", "strand", "identity",
            "is_cov", "is_unaligned_head", "is_unaligned_tail", "left_extra",
            "right_extra", "extra", "delta", "class", "other_is_segments"]
    per = collections.defaultdict(lambda: collections.Counter())
    extras = collections.defaultdict(collections.Counter)
    species = collections.defaultdict(set)
    total = collections.Counter()
    with open(os.path.join(a.out, "per_event.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for q in sorted(ev):
            e = ev[q]
            if q not in best:
                total["no_is"] += 1
                continue
            (nm, ident, cov, qlen, qs, qe, strand, t, tlen, ts, te), others = best[q]
            left, right = place(qlen, qs, qe, strand, tlen, ts, te)
            cls, delta = classify(left, right, e["junction_overlap_bp"], a.near)
            name, fam, grp = is_parts(t)
            w.writerow([e["species"], e["event_id"], qlen, e["junction_overlap_bp"],
                        e["insert_frame_status"], e["decomposition_method"],
                        e["S1_structurally_clean"], name, fam, grp, tlen, strand,
                        "%.4f" % ident, "%.4f" % cov, ts, tlen - te, left, right,
                        left + right, delta, cls, others])
            total[cls] += 1
            per[name][cls] += 1
            if others:
                per[name]["with_other_is"] += 1
            if cls in ("exact", "exact_shifted", "near"):
                extras[name][left + right] += 1
            species[name].add(e["species"])

    pcols = ["is_name", "is_family", "is_len", "in_panel", "tsd_literature",
             "n_judged", "exact", "exact_shifted", "near", "over", "under",
             "pct_exact", "pct_exact_or_shifted", "modal_extra_bp", "modal_extra_share",
             "with_other_is", "species"]
    fam_len = {is_parts(t)[0]: (is_parts(t)[1], n) for t, n in is_len.items()}
    with open(os.path.join(a.out, "per_is.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(pcols)
        for name in sorted(per, key=lambda n: -sum(per[n][c] for c in CLASSES)):
            c = per[name]
            n = sum(c[k] for k in CLASSES)
            if n < a.min_n and name not in panel:
                continue
            ex = extras[name]
            mode, share = "", ""
            if ex:
                mode, k = max(ex.items(), key=lambda kv: (kv[1], -kv[0]))
                share = "%.3f" % (k / sum(ex.values()))
            w.writerow([name, fam_len[name][0], fam_len[name][1],
                        "yes" if name in panel else "no",
                        panel.get(name, {}).get("tsd_literature", ""),
                        n, c["exact"], c["exact_shifted"], c["near"], c["over"],
                        c["under"], "%.1f" % (100 * c["exact"] / n),
                        "%.1f" % (100 * (c["exact"] + c["exact_shifted"]) / n),
                        mode, share, c["with_other_is"],
                        ",".join(sorted(species[name]))])

    judged = sum(total[c] for c in CLASSES)
    res = {
        "status": "complete", "git_commit": commit, "git_dirty": dirty,
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "isfinder": {"path": os.path.realpath(a.isfinder), "sha256": got,
                     "sequences": len(is_len), "identical_repeats_dropped": dropped},
        "panel": {"path": os.path.realpath(a.panel), "sha256": sha256_file(a.panel)},
        "events": inputs, "minimap2": {"version": ver, "cmd": cmd},
        "params": {"min_ident": a.min_ident, "min_cov": a.min_cov, "near": a.near},
        "inserts": len(ev), "judged": judged, "classes": dict(total),
    }
    with open(os.path.join(a.out, "result.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("inserts", "judged", "classes")}))


if __name__ == "__main__":
    main()

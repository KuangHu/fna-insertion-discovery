#!/usr/bin/env python3
"""LINE M -- join Layer-2 inserts to Arm A mobility evidence. No new compute.

Arm A has produced 1,610 within-genome multi-copy families over 200 genomes
since the `-f 0` fix and has NEVER been joined to Layer 2. This is that join.

It answers the acceptance question directly: an insert is element-like if it is
multi-copy at DISPERSED loci inside one genome, or identical to an insert in
another genome. Both are family-agnostic and neither needs polarity, a panel or
a tree.

    usable = reconstruction clean  AND  mobility positive

WHAT THIS IS NOT. "Multi-copy at dispersed loci" is satisfied by rRNA operons
(7 copies in E. coli), prophage, REP elements and segmental duplications. This
module reports the evidence; it does NOT assert the insert is a mobile element.
The false-positive rate must be measured at Layer 5, after the fact, and if it
is high the criterion needs a filter. Until that test runs, `mobility_positive`
means "has the structural signature", not "is an IS".

Two independent sub-lines, reported separately because they can disagree:

  M1  within-genome  the insert maps to an Arm A family in its own carrier
                     with n_distinct_loci >= --min-loci. Arm A's own
                     flank_homology / internal identity are carried through so
                     a tandem array can be told from a dispersed family.
  M2  cross-EVENT    a byte-identical insert (inserted_md5) occurs at a
                     DIFFERENT EVENT (from 100_event_dedup.py), never merely a
                     different locus_id.

                     **locus_id IS PANEL-SCOPED, NOT A GLOBAL IDENTITY.** The
                     same genomic insertion seen by two panels carries two
                     locus_ids and one insert md5, so a locus-level test scores
                     it as "the same insert at another locus" -- which is panel
                     redundancy, not mobility. Measured on the K. pneumoniae
                     400-panel run: **85.7% by locus_id vs 35.2% by event_id**,
                     an inflation of 50.4 points over 15,444 loci. At 40 panels
                     the same bug read 42.8%. This is why the dedup must run
                     BEFORE the M join.

                     The earlier form of this same error was subtler: the very
                     first version tested "present in another GENOME", which at
                     >=99% ANI is nearly free -- 26 of 28 positives had the
                     insert only at their own locus. "Present in another genome" is NOT the
                     criterion: at >=99% ANI two genomes share almost
                     everything byte-identically, so same-locus identity is
                     shared ancestry, not mobility. Measured: of 28 loci
                     passing the naive genome-based test, 26 had the md5 only
                     at their own locus, and ALL 19 of its unique rescues were
                     same-locus. The locus constraint is what makes M2 mean
                     anything.

    99_mobility_join.py --recon <recon dirs> --arma armA_f0 --out mobility
"""
import argparse
import collections
import csv
import glob
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert                                                 # noqa: E402
from lib.util import log                                          # noqa: E402


def read_alleles(path):
    out, cur = collections.defaultdict(dict), None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "").upper() for a, v in d.items()}
            for k, d in out.items()}



def expand(paths):
    """Expand glob patterns that survived the shell.

    `slurm/analysis.sh` sets `set -f` so patterns reach the script intact
    (without it the submitting shell expands them and argparse sees dozens of
    positional args). The cost is that any script taking a directory list MUST
    glob for itself -- and when a chained job is submitted BEFORE its input
    directories exist, the pattern is all it ever gets. 100_event_dedup.py once
    read 0 panels this way and the whole A. baumannii chain stalled.
    """
    out = []
    for p in paths:
        out.extend(sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p])
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--arma", required=True)
    ap.add_argument("--events", default=None,
                    help="loci_to_event.tsv from 100_event_dedup.py. REQUIRED "
                         "for a correct M2: without it M2 falls back to "
                         "locus_id, which is panel-scoped and inflates the rate "
                         "by ~50 points on a 400-panel run.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-loci", type=int, default=3,
                    help="distinct loci an Arm A family must occupy for M1")
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    args.recon = expand(args.recon)

    loci, aseq, carriers = {}, {}, collections.defaultdict(dict)
    for d in args.recon:
        # SKIP a panel directory with no loci.tsv. The Layer 2 driver mkdir's
        # the output dir before checking that the panel has any seed loci, so a
        # panel with none leaves an EMPTY directory behind. 100_event_dedup
        # already skipped these; 99_ and 86_ did not, and both crashed with
        # FileNotFoundError after the dedup had already succeeded.
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            loci[r["locus_id"]] = r
        aseq.update(read_alleles(os.path.join(d, "alleles.fna")))
        p = os.path.join(d, "alleles.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                carriers[r["locus_id"]][r["allele_id"]] = \
                    [g for g in r.get("genomes", "").split(",") if g]
    log("%d Layer-2 loci" % len(loci))

    # ---- the decomposable set: the only loci with an insert to test --------
    dec = {}
    for lid, r in loci.items():
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        try:
            lcp, ilen = int(r["lcp_bp"]), int(r["inserted_len"])
        except ValueError:
            continue
        if ilen <= 0:
            continue
        ls = aseq.get(lid, {}).get(r["longest_allele"], "")
        if not ls:
            continue
        # NOT ls[lcp:lcp+ilen]: lcp_bp is a SHORT-allele offset and using it
        # here mis-sliced 31.1% of tolerant-path loci, which then made M2 ask
        # its question of the wrong sequence.
        i0, frame = lib_insert.locate_insert(r, ls)
        iseq = ls[i0:i0 + ilen] if i0 >= 0 else ""
        key = lib_insert.canonical_insert_key(ls, i0, ilen) if i0 >= 0 else ""
        dec[lid] = {"seq": iseq, "m2key": key, "frame": frame, "row": r,
                    "carriers": carriers.get(lid, {}).get(r["longest_allele"], [])}
    log("%d decomposable loci with an insert" % len(dec))

    # ---- M2: byte-identical insert in another genome -----------------------
    ev_of = {}
    if args.events and os.path.exists(args.events):
        for r in csv.DictReader(open(args.events), delimiter="\t"):
            ev_of[r["locus_id"]] = r["event_id"]
        n_map = sum(1 for l in dec if l in ev_of)
        log("event map: %d of %d decomposable loci carry an event_id"
            % (n_map, len(dec)))
        if dec and n_map < 0.5 * len(dec):
            log("FATAL: the event map covers %.1f%% of loci -- M2 would silently "
                "fall back to panel-scoped locus ids and inflate"
                % (100.0 * n_map / len(dec)))
            return 2
    else:
        log("WARNING: no --events file. M2 will use PANEL-SCOPED locus ids and "
            "will be INFLATED (measured: 85.7%% vs 35.2%% on 400 panels). Run "
            "100_event_dedup.py first and pass its loci_to_event.tsv.")

    def unit(lid):
        """Global identity of a locus: its event when known, else the locus id."""
        return ev_of.get(lid, lid)

    # KEY: canonical, not `inserted_md5`. inserted_md5 hashes the insert as
    # stored, so the same element captured on opposite strands in two
    # assemblies hashes to two values and M2 answers False for both. The
    # question M2 asks -- "is this element also somewhere else" -- has no
    # strand in it. lib_insert.canon_key hashes min(seq, revcomp(seq)).
    by_md5_unit = collections.defaultdict(set)
    by_md5_genome = collections.defaultdict(set)
    n_nokey = 0
    for lid, d in dec.items():
        key = d["m2key"]
        if key:
            by_md5_unit[key].add(unit(lid))
            for g in d["carriers"]:
                by_md5_genome[key].add(g)
        else:
            n_nokey += 1
    if n_nokey:
        log("  %d loci carry no verifiable insert sequence; M2 is False for "
            "them by absence of evidence, not by evidence of absence" % n_nokey)
    for lid, d in dec.items():
        key = d["m2key"]
        other = by_md5_unit.get(key, set()) - {unit(lid)} if key else set()
        d["m2_loci"] = len(by_md5_unit.get(key, set())) if key else 0
        d["m2_genomes"] = len(by_md5_genome.get(key, set())) if key else 0
        d["m2"] = len(other) > 0     # a DIFFERENT EVENT, not another locus_id


    # ---- M1: map each insert to its carrier's Arm A families ---------------
    fam = {}
    if not os.path.isdir(args.arma):
        # sa_k12 shares sa_k24's Arm A -- the k-experiment arms differ only in
        # panel size -- so a missing armA dir is a wiring mistake in the caller,
        # not a species with no multi-copy families. Say which, and stop: M1
        # silently False for every locus reads exactly like a real negative.
        log("FATAL: --arma %s does not exist. M1 cannot be computed and an "
            "all-False M1 column is indistinguishable from a true negative."
            % args.arma)
        return 5
    for f in os.listdir(args.arma):
        if not f.endswith("_families.tsv"):
            continue
        g = f[:-len("_families.tsv")]
        for r in csv.DictReader(open(os.path.join(args.arma, f)), delimiter="\t"):
            fam[(g, r["family_id"])] = r
    log("%d Arm A families loaded" % len(fam))

    by_genome = collections.defaultdict(list)
    for lid, d in dec.items():
        for g in d["carriers"]:
            by_genome[g].append(lid)

    tmp = tempfile.mkdtemp(prefix="mjoin.")
    hits = collections.defaultdict(list)
    for g, lids in sorted(by_genome.items()):
        el = os.path.join(args.arma, g + "_elements.fna")
        if not os.path.exists(el):
            continue
        q = os.path.join(tmp, g + ".q.fna")
        with open(q, "w") as fh:
            for lid in lids:
                fh.write(">%s\n%s\n" % (lid, dec[lid]["seq"]))
        paf = os.path.join(tmp, g + ".paf")
        with open(paf, "w") as fh:
            subprocess.run([args.minimap2, "-c", "-x", "asm10", "-t",
                            str(args.threads), el, q], stdout=fh,
                           stderr=subprocess.DEVNULL, check=False)
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qn, ql, tn = f[0], int(f[1]), f[5]
            nm, bl = int(f[9]), int(f[10])
            ident = nm / bl if bl else 0.0
            cov = bl / ql if ql else 0.0
            if ident >= args.min_ident and cov >= args.min_cov:
                # family_id in *_families.tsv is the FULL "<genome>.A0000"
                # string, and the elements.fna record name is that same string.
                # Splitting on "." to take the suffix produced a key that never
                # matched, so M1 silently returned 0 for every locus.
                hits[qn].append((ident, cov, g, tn))
    # matched/total for the ARM A FAMILY JOIN -- the join that once returned
    # M1 = 0 for every locus because `split(".")[-1]` produced a family key that
    # never matched. Seven wrong-key failures in this project, no two on the
    # same key: the defence cannot be per-key, it has to be that every join
    # reports its match rate by default.
    n_hit = sum(1 for lid in dec if hits.get(lid))
    n_fam = 0
    for lid in dec:
        for _, _, g, fid in hits.get(lid, []):
            if (g, fid) in fam:
                n_fam += 1
                break
    log("arm A join: %d/%d inserts aligned to an element (%.1f%%); "
        "%d/%d resolved to a family record (%.1f%%)"
        % (n_hit, len(dec), 100.0 * n_hit / len(dec) if dec else 0,
           n_fam, len(dec), 100.0 * n_fam / len(dec) if dec else 0))
    if n_hit and n_fam == 0:
        log("FATAL: %d inserts aligned to an Arm A element but NONE resolved to "
            "a family record -- the family key does not match "
            "*_families.tsv.family_id" % n_hit)
        return 2

    for lid, d in dec.items():
        best = sorted(hits.get(lid, []), reverse=True)[:1]
        d["m1"] = False
        d["m1_family"] = d["m1_copies"] = d["m1_loci"] = "."
        d["m1_flank_homology"] = d["m1_internal_id"] = "."
        if best:
            ident, cov, g, fid = best[0]
            r = fam.get((g, fid))
            if r:
                nl = int(r.get("n_distinct_loci") or 0)
                d["m1_family"] = "%s.%s" % (g, fid)
                d["m1_copies"] = r.get("n_copies")
                d["m1_loci"] = nl
                d["m1_flank_homology"] = r.get("flank_homology")
                d["m1_internal_id"] = r.get("internal_identity_median")
                d["m1"] = nl >= args.min_loci

    # ---- usable = reconstruction clean AND mobility positive ---------------
    rows = []
    for lid, d in sorted(dec.items()):
        r = d["row"]
        amb = int(r.get("junction_ambiguity_bp") or 0)
        lost = int(r.get("target_bases_lost") or 0)
        clean = (lost == 0 and amb <= 15
                 and r.get("placement_status", "ok") not in ("ambiguous",))
        mob = d["m1"] or d["m2"]
        rows.append({
            "locus_id": lid, "inserted_len": r["inserted_len"],
            "junction_overlap_bp": amb, "target_bases_lost": lost,
            "decomposition_method": r.get("decomposition_method"),
            "recon_clean": clean,
            "M1_within_genome": d["m1"], "M1_family": d["m1_family"],
            "M1_copies": d["m1_copies"], "M1_distinct_loci": d["m1_loci"],
            "M1_flank_homology": d["m1_flank_homology"],
            "M1_internal_identity": d["m1_internal_id"],
            "M2_cross_event": d["m2"], "M2_n_events_with_md5": d["m2_loci"],
            "M2_n_genomes": d["m2_genomes"],
            "mobility_positive": mob, "USABLE": clean and mob})
    cols = list(rows[0].keys()) if rows else []
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_mobility.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n = len(rows)
    c = collections.Counter()
    for r in rows:
        c["recon_clean"] += r["recon_clean"]
        c["M1"] += r["M1_within_genome"]
        c["M2"] += r["M2_cross_event"]
        c["M1_only"] += r["M1_within_genome"] and not r["M2_cross_event"]
        c["M2_only"] += r["M2_cross_event"] and not r["M1_within_genome"]
        c["both"] += r["M1_within_genome"] and r["M2_cross_event"]
        c["mob"] += r["mobility_positive"]
        c["usable"] += r["USABLE"]
    log("=" * 84)
    log("LINE M -- MOBILITY JOIN")
    log("")
    log("  decomposable loci                 %4d" % n)
    log("  reconstruction clean              %4d   %5.1f%%" % (c["recon_clean"], 100.0*c["recon_clean"]/n))
    log("")
    log("  M1 within-genome (>= %d loci)      %4d   %5.1f%%" % (args.min_loci, c["M1"], 100.0*c["M1"]/n))
    log("  M2 identical insert at ANOTHER LOCUS %4d   %5.1f%%" % (c["M2"], 100.0*c["M2"]/n))
    log("    M1 only %d   M2 only %d   both %d" % (c["M1_only"], c["M2_only"], c["both"]))
    log("  mobility positive (M1 or M2)      %4d   %5.1f%%" % (c["mob"], 100.0*c["mob"]/n))
    log("")
    log("  *** USABLE = clean AND mobile     %4d   %5.1f%% of decomposable ***"
        % (c["usable"], 100.0*c["usable"]/n))
    log("")
    log("  NOT VALIDATED: rRNA operons, prophage, REP elements and segmental")
    log("  duplications also satisfy multi-copy-at-dispersed-loci. The")
    log("  false-positive rate must be measured at Layer 5 before")
    log("  `mobility_positive` may be read as `is a mobile element`.")
    log("=" * 84)
    log("wrote %s_mobility.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

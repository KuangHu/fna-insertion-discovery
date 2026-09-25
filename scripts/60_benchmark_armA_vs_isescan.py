#!/usr/bin/env python3
"""Benchmark 1 -- Arm A vs ISEScan as a CONTINGENCY TABLE, not a recall number.

ISEScan is not truth. It finds IS by transposase pHMM + terminal repeats, so it
carries exactly the family bias this pipeline exists to avoid, and a single
"recall vs ISEScan" figure silently treats its blind spots as Arm A's errors.
The informative object is the 2x2:

                    ISEScan +                 ISEScan -
    Arm A +    known IS recovered      ANNOTATION-FREE CANDIDATES  <- the point
    Arm A -    missed known IS         background (not enumerable)

Read it in both directions:

  * ArmA+ / ISEScan+   Arm A found it on structure alone, with no HMM. This is
                       the sanity cell.
  * ArmA+ / ISEScan-   NOT automatically false positive. IS110/IS1111 has no
                       TIR and no TSD -- the two features ISEScan keys on --
                       so an annotation-free method is EXPECTED to win here.
                       This cell is the discovery pool and is characterised
                       separately in <out>_novel_pool.tsv.
  * ArmA- / ISEScan+   genuine sensitivity failure. Split by ISEScan copy
                       number, because Arm A never claimed singletons: it needs
                       >=3 copies at distinct loci. Singleton misses are a
                       declared blind spot, multi-copy misses are real bugs.
  * ArmA- / ISEScan-   the rest of the genome. Not countable without a unit of
                       "negative locus", so it is reported as "-", per design.

The novel pool is NOT graded here. Deciding whether an ArmA+/ISEScan- family is
real needs the stage 50 screens (rRNA via barrnap, ORF fraction, transposase
HMM), and running them here would reintroduce the annotation-first ordering the
pipeline is built to avoid. This script only isolates and describes the pool.

    60_benchmark_armA_vs_isescan.py \\
        --armA-copies-glob   'armA_5kb/*_copies.tsv' \\
        --armA-families-glob 'armA_5kb/*_families.tsv' \\
        --isescan-glob       'isescan_run/*/*.tsv' \\
        --out bench60
"""
import argparse
import csv
import glob
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log, overlap                               # noqa: E402

# ISEScan reports these as "partial"/"tir" etc.; a call with no transposase hit
# is not a confident IS and is excluded from the ISEScan+ column by default.
ISESCAN_MIN_LEN = 100


def sample_from_path(path):
    """GCA_000258145.1_copies.tsv -> GCA_000258145.1 ; one.fna.tsv -> one."""
    b = os.path.basename(path)
    for suf in ("_copies.tsv", "_families.tsv", ".fna.tsv", ".tsv",
                ".fna.csv", ".csv", ".fna"):
        if b.endswith(suf):
            b = b[: -len(suf)]
            break
    return b


def read_isescan(path, sample):
    """ISEScan tsv/csv: seqID, family, cluster, isBegin, isEnd, type, ..."""
    out = []
    with open(path) as fh:
        head = fh.readline()
        delim = "," if head.count(",") > head.count("\t") else "\t"
        fh.seek(0)
        for r in csv.DictReader(fh, delimiter=delim):
            try:
                s, e = int(r["isBegin"]), int(r["isEnd"])
            except (KeyError, ValueError, TypeError):
                continue
            if abs(e - s) < ISESCAN_MIN_LEN:
                continue
            out.append({"sample": sample, "contig": r.get("seqID", "."),
                        "start": min(s, e), "end": max(s, e),
                        "family": r.get("family", "."),
                        "cluster": r.get("cluster", "."),
                        "type": r.get("type", ".")})
    return out


def read_armA_copies(path):
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                yield {"sample": r["sample_id"], "contig": r["contig"],
                       "start": int(r["start"]), "end": int(r["end"]),
                       "family_id": r["family_id"],
                       "tier": r.get("qc_tier", ".")}
            except (KeyError, ValueError):
                continue


def read_armA_families(path):
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield r


def fnum(row, key, default=float("nan")):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return default


def med(vals):
    vals = [v for v in vals if v == v]                      # drop NaN
    return statistics.median(vals) if vals else float("nan")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--armA-copies-glob", required=True)
    ap.add_argument("--armA-families-glob", required=True)
    ap.add_argument("--isescan-glob", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="reciprocal overlap to count a locus match (default 0.5)")
    ap.add_argument("--multi-copy-min", type=int, default=3,
                    help="ISEScan copies at/above this count as multi-copy, "
                         "which is what Arm A actually targets (default 3)")
    ap.add_argument("--strong-only", action="store_true",
                    help="restrict Arm A to STRONG_MOBILE_CANDIDATE families")
    args = ap.parse_args()

    # ------------------------------------------------------------- load
    isc = []
    for p in sorted(glob.glob(args.isescan_glob)):
        isc.extend(read_isescan(p, sample_from_path(p)))
    fam = {}
    for p in sorted(glob.glob(args.armA_families_glob)):
        for r in read_armA_families(p):
            fam[r["family_id"]] = r
    copies = []
    for p in sorted(glob.glob(args.armA_copies_glob)):
        copies.extend(read_armA_copies(p))

    if args.strong_only:
        keep = {k for k, v in fam.items()
                if v.get("verdict") == "STRONG_MOBILE_CANDIDATE"}
        copies = [c for c in copies if c["family_id"] in keep]
        fam = {k: v for k, v in fam.items() if k in keep}

    # Only genomes with BOTH sides run can enter the table. Counting an Arm A
    # family from a genome ISEScan never saw would inflate the novel cell.
    s_isc = {h["sample"] for h in isc}
    s_arm = {c["sample"] for c in copies}
    shared = s_isc & s_arm
    if not shared:
        log("FATAL: no genome has both Arm A and ISEScan output.")
        log("  ISEScan samples: %s" % (sorted(s_isc)[:5] or "none"))
        log("  Arm A samples:   %s" % (sorted(s_arm)[:5] or "none"))
        return 1
    dropped = (s_isc | s_arm) - shared
    isc = [h for h in isc if h["sample"] in shared]
    copies = [c for c in copies if c["sample"] in shared]
    fam = {k: v for k, v in fam.items()
           if v.get("sample_id", k.split(".A")[0]) in shared}
    log("%d genomes with both sides (%d dropped); %d IS calls, %d Arm A copies, "
        "%d families" % (len(shared), len(dropped), len(isc), len(copies), len(fam)))

    # ------------------------------------------------------- locus matching
    per_family_n = defaultdict(int)
    for h in isc:
        per_family_n[(h["sample"], h["family"], h["cluster"])] += 1

    by_contig = defaultdict(list)
    for c in copies:
        by_contig[(c["sample"], c["contig"])].append(c)

    is_rows = []
    armA_hit = set()                       # family_ids overlapping some IS call
    hit = {"multi": 0, "singleton": 0}
    tot = {"multi": 0, "singleton": 0}
    for h in isc:
        best, best_r = None, 0.0
        for c in by_contig.get((h["sample"], h["contig"]), ()):
            ov = overlap(h["start"], h["end"], c["start"], c["end"])
            if ov == 0:
                continue
            r = min(ov / max(1, h["end"] - h["start"]),
                    ov / max(1, c["end"] - c["start"]))
            if r > best_r:
                best, best_r = c, r
        found = best is not None and best_r >= args.min_overlap
        n_cop = per_family_n[(h["sample"], h["family"], h["cluster"])]
        cls = "multi" if n_cop >= args.multi_copy_min else "singleton"
        tot[cls] += 1
        hit[cls] += found
        if found:
            armA_hit.add(best["family_id"])
        is_rows.append([h["sample"], h["contig"], h["start"], h["end"],
                        h["family"], h["type"], n_cop, cls, int(found),
                        "%.3f" % best_r, best["family_id"] if best else "."])

    novel = sorted(set(fam) - armA_hit)
    recovered = sorted(set(fam) & armA_hit)

    # ------------------------------------------------------------- 2x2
    cells = [
        ["ARM_A+", "ISESCAN+", len(recovered), "families",
         "known IS recovered with no HMM"],
        ["ARM_A+", "ISESCAN-", len(novel), "families",
         "ANNOTATION-FREE CANDIDATES -- the discovery pool"],
        ["ARM_A-", "ISESCAN+_multi", tot["multi"] - hit["multi"], "IS calls",
         "sensitivity failure (Arm A targets multi-copy)"],
        ["ARM_A-", "ISESCAN+_singleton", tot["singleton"] - hit["singleton"],
         "IS calls", "declared blind spot, not a grade"],
        ["ARM_A-", "ISESCAN-", "-", "-", "not enumerable by design"],
    ]
    with open(args.out + "_contingency.tsv", "w") as fh:
        fh.write("arm_a\tisescan\tn\tunit\tnote\n")
        for c in cells:
            fh.write("\t".join(map(str, c)) + "\n")

    with open(args.out + "_per_is.tsv", "w") as fh:
        fh.write("sample\tcontig\tis_start\tis_end\tis_family\tis_type\t"
                 "is_copy_number\tclass\tfound_by_armA\treciprocal_overlap\t"
                 "armA_family\n")
        for r in is_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    # -------------------------------------------- characterise the novel pool
    pool_cols = ["family_id", "sample_id", "n_copies", "n_distinct_loci",
                 "elem_len_median", "internal_identity_median",
                 "boundary_disp_L_bp", "boundary_disp_R_bp", "boundary_sharpness",
                 "S_mobility", "verdict", "elem_md5"]
    with open(args.out + "_novel_pool.tsv", "w") as fh:
        fh.write("\t".join(pool_cols) + "\n")
        for fid in novel:
            r = fam[fid]
            fh.write("\t".join(str(r.get(c, "")) for c in pool_cols) + "\n")

    def describe(ids, label):
        if not ids:
            return "%-28s n=0" % label
        rows = [fam[i] for i in ids]
        vd = defaultdict(int)
        for r in rows:
            vd[r.get("verdict", "?")] += 1
        top = ", ".join("%s=%d" % kv for kv in
                        sorted(vd.items(), key=lambda x: -x[1])[:3])
        return ("%-28s n=%-5d len=%-6.0f copies=%-4.1f loci=%-4.1f "
                "ident=%-6.2f dispL=%-5.1f S=%-5.2f  %s"
                % (label, len(rows),
                   med([fnum(r, "elem_len_median") for r in rows]),
                   med([fnum(r, "n_copies") for r in rows]),
                   med([fnum(r, "n_distinct_loci") for r in rows]),
                   med([fnum(r, "internal_identity_median") for r in rows]),
                   med([fnum(r, "boundary_disp_L_bp") for r in rows]),
                   med([fnum(r, "S_mobility") for r in rows]), top))

    def pct(a, b):
        return 100.0 * a / b if b else float("nan")

    log("=" * 78)
    log("CONTINGENCY  (Arm A families x ISEScan calls, %d shared genomes)"
        % len(shared))
    log("")
    log("                        ISEScan +            ISEScan -")
    log("  Arm A +      %8d families    %8d families   <- discovery pool"
        % (len(recovered), len(novel)))
    log("  Arm A -      %8d IS multi     %8s"
        % (tot["multi"] - hit["multi"], "-"))
    log("               %8d IS single    %8s"
        % (tot["singleton"] - hit["singleton"], "(not enumerable)"))
    log("")
    log("recall on multi-copy IS (what Arm A targets): %d/%d = %.1f%%"
        % (hit["multi"], tot["multi"], pct(hit["multi"], tot["multi"])))
    log("recall on singleton IS  (declared blind spot): %d/%d = %.1f%%"
        % (hit["singleton"], tot["singleton"], pct(hit["singleton"],
                                                   tot["singleton"])))
    log("-" * 78)
    log(describe(recovered, "ArmA+/ISEScan+ (known)"))
    log(describe(novel, "ArmA+/ISEScan- (NOVEL)"))
    log("-" * 78)
    log("The novel cell is NOT graded here. Run the stage 50 structural screens")
    log("(barrnap rRNA, ORF fraction, transposase HMM) on %s_novel_pool.tsv"
        % args.out)
    log("to separate rrn-like repeats from genuine annotation-free elements.")
    log("=" * 78)
    log("wrote %s_{contingency,per_is,novel_pool}.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

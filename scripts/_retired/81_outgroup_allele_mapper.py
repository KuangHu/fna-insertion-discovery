#!/usr/bin/env python3
"""Layer 3b -- observe which allele each outgroup carries. Still no polarity.

The outgroup's only job here is to be SEEN at the locus. It does not vote, it
is not compared to a reference, and no allele is called ancestral. The single
question is:

    the interval between the same two anchors, in this outgroup --
    is it one of the alleles the ingroup already showed, or a new one?

Matching is SEQUENCE-based, never length-based: two 40 bp intervals that differ
at every base are not the same allele, and a 1279 bp element with a handful of
SNPs is. Comparison is orientation-agnostic, because an outgroup contig can
carry the locus on either strand and a forward-only comparison returns ~0.50
identity for a perfect reverse-complement match -- an artefact that already
faked 6 "errors" once in this project.

Locus placement REUSES the locked Layer 2 `place_locus()` unchanged. Writing a
second locus mapper for outgroups would mean two code paths that could disagree
about where a locus is, and the outgroup path would be the unvalidated one.

allele_match_status, per outgroup per locus:
    exact_existing_allele   identical to a known ingroup allele
    near_existing_allele    >= --near-ident / --near-cov to a known allele
    novel_outgroup_allele   placed cleanly, matches nothing the ingroup has
    unresolved              anchors did not place, or placement was ambiguous

`unresolved` is a first-class outcome and is reported, never quietly dropped:
the fraction of loci with >=2 confidently placed outgroups is exactly the
"outgroup callability" that bounds how much of the catalogue Layer 3 can cover.

    81_outgroup_allele_mapper.py --loci allele_recon/armB_cl0000/loci.tsv \\
        --alleles allele_recon/armB_cl0000/alleles.fna \\
        --seed-loci allele_recon/loci_armB_cl0000.tsv \\
        --outgroups outgroup/cl0000/outgroups.manifest --outdir og/cl0000
"""
import argparse
import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp, write_fasta                # noqa: E402
from lib.util import log                                         # noqa: E402
from Bio.Align import PairwiseAligner                            # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_ar", os.path.join(_HERE, "70_allele_reconstructor.py"))
_ar = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ar)


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    return b


def read_alleles(path):
    """locus_id -> {allele_id: sequence}"""
    out, cur = defaultdict(dict), None
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "") for a, v in d.items()}
            for k, d in out.items()}


def _corr(a, b):
    if not a and not b:
        return 1.0, 1.0
    if not a or not b:
        return 0.0, 0.0
    al = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                         open_gap_score=-12, extend_gap_score=-0.4)
    try:
        aln = al.align(a, b)[0]
    except (ValueError, OverflowError, MemoryError, IndexError):
        return 0.0, 0.0
    sb, lb = aln.aligned
    cols = matches = 0
    for (s0, s1), (l0, l1) in zip(sb, lb):
        matches += sum(1 for x, y in zip(a[s0:s1], b[l0:l1]) if x == y)
        cols += s1 - s0
    if cols == 0:
        return 0.0, 0.0
    return matches / cols, cols / min(len(a), len(b))


def corr(a, b):
    """Orientation-agnostic: an outgroup contig may carry either strand."""
    f = _corr(a, b)
    r = _corr(a, revcomp(b))
    return f if f[0] >= r[0] else r


def _kmers(s, k):
    return {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else set()


def _containment(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def match_allele(seq, allele_items, k=15, top=2, min_len_ratio=0.34):
    """Best and runner-up ingroup allele for an outgroup interval.

    The naive version aligned every outgroup interval against EVERY allele in
    BOTH orientations, which is O(loci x outgroups x alleles x 2) full local
    alignments and made the sensitivity grid a 10-hour job. Three cheap stages
    remove almost all of that work without changing the answer:

      1. exact hash, forward and reverse-complement -- an exact allele match
         needs no alignment at all, and a third of placements are exact;
      2. length ratio -- a 40 bp interval cannot be a near-match to a 1319 bp
         allele, whatever the aligner says;
      3. k-mer containment to shortlist, then align only the top `top`
         candidates, and only in the orientation k-mers already chose.

    `allele_items` is a list of (allele_id, seq, kmers_fwd) prepared once per
    locus so the allele k-mer sets are not rebuilt for every outgroup.

    Returns [(identity, coverage, allele_id), ...] sorted best first.
    """
    if not allele_items:
        return []
    exact = {}
    rc = revcomp(seq) if seq else ""
    for aid, aseq, _ in allele_items:
        if seq == aseq or (rc and rc == aseq):
            exact[aid] = (1.0, 1.0, aid)
    ql = len(seq)
    qk_f, qk_r = _kmers(seq, k), _kmers(rc, k)
    scored = []
    for aid, aseq, ak in allele_items:
        if aid in exact:
            continue
        al = len(aseq)
        if max(ql, al) and min(ql, al) / max(ql, al) < min_len_ratio:
            continue
        if ak and (qk_f or qk_r):
            c = max(_containment(qk_f, ak), _containment(qk_r, ak))
        else:
            c = 1.0                       # too short to k-mer; let alignment decide
        scored.append((c, aid, aseq))
    scored.sort(key=lambda x: -x[0])
    out = list(exact.values())
    for c, aid, aseq in scored[:top]:
        i, cov = corr(seq, aseq)
        out.append((i, cov, aid))
    out.sort(key=lambda x: (-x[0], -x[1]))
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True, help="Layer 2 loci.tsv")
    ap.add_argument("--alleles", required=True, help="Layer 2 alleles.fna")
    ap.add_argument("--seed-loci", required=True,
                    help="the seed TSV Layer 2 was run on (for anchor coords)")
    ap.add_argument("--seed-genomes", required=True,
                    help="manifest of the ingroup, to fetch seed anchors from")
    ap.add_argument("--outgroups", required=True, help="outgroups.manifest")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--anchor", type=int, default=500)
    ap.add_argument("--max-interval", type=int, default=7000)
    ap.add_argument("--min-ident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--margin-ratio", type=float, default=0.05)
    ap.add_argument("--near-ident", type=float, default=0.95)
    ap.add_argument("--near-cov", type=float, default=0.90)
    ap.add_argument("--allele-margin", type=float, default=0.005,
                    help="best allele must beat the runner-up by this identity "
                         "margin. 61%% of outgroup placements are near-matches "
                         "rather than exact, and at 98%% ANI ordinary SNPs make "
                         "identity to two similar alleles nearly equal "
                         "(.991 vs .989). Without a margin those become firm "
                         "votes in a majority rule.")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    loci = {}
    with open(args.loci) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            loci[r["locus_id"]] = r
    seeds = {}
    with open(args.seed_loci) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seeds[r["locus_id"]] = r
    alleles = read_alleles(args.alleles)
    ing = {sample_name(p): p for p in
           (l.strip() for l in open(args.seed_genomes) if l.strip())}
    ogs = [l.strip() for l in open(args.outgroups) if l.strip()]
    log("%d loci, %d outgroups" % (len(loci), len(ogs)))
    if not ogs:
        log("no outgroups supplied -- nothing to map")
        return 1

    # rebuild the same anchors Layer 2 used, from the same seed genome
    tmp = tempfile.mkdtemp(prefix="ogmap.")
    qpath = os.path.join(tmp, "anchors.fna")
    fas, kept = {}, []
    with open(qpath, "w") as qf:
        for lid, r in sorted(loci.items()):
            sd = seeds.get(lid)
            if not sd:
                continue
            gp = ing.get(sd["seed_genome"])
            if not gp:
                continue
            fa = fas.setdefault(sd["seed_genome"], Fasta(gp))
            if sd["contig"] not in fa:
                continue
            a, b = int(sd["start"]), int(sd["end"])
            clen = fa.length(sd["contig"])
            if a - args.anchor < 0 or b + args.anchor > clen:
                continue
            write_fasta(qf, lid + "#L", fa.fetch(sd["contig"], a - args.anchor, a))
            write_fasta(qf, lid + "#R", fa.fetch(sd["contig"], b, b + args.anchor))
            kept.append(lid)
    log("%d loci with reconstructable anchors" % len(kept))

    obs = defaultdict(dict)
    for gp in ogs:
        g = sample_name(gp)
        paf = os.path.join(tmp, g + ".paf")
        try:
            _ar.map_anchors(qpath, gp, paf, args.threads, args.minimap2)
        except subprocess.CalledProcessError:
            log("  minimap2 failed on %s" % g)
            continue
        hits = _ar.best_hits(paf, args.min_ident, args.min_cov)
        fa = Fasta(gp)
        n = 0
        for lid in kept:
            lh, rh = hits.get(lid + "#L", []), hits.get(lid + "#R", [])
            if not lh or not rh:
                continue
            # no localisation prior: an outgroup's coordinates are unrelated to
            # the seed's, so only the pair-scoring and margin apply
            pl = _ar.place_locus(lh, rh, args.max_interval, None, 20000,
                                 args.margin_ratio)
            if pl is None:
                continue
            seq = (fa.fetch(pl["contig"], pl["start"], pl["end"])
                   if pl["end"] > pl["start"] else "")
            if pl["strand"] == "-":
                seq = revcomp(seq)
            obs[lid][g] = (seq, pl)
            n += 1
        log("  %-22s placed %d/%d loci" % (g, n, len(kept)))

    # allele k-mer sets built ONCE per locus, not once per outgroup
    kidx = {}
    for lid in kept:
        kidx[lid] = [(aid, aseq, _kmers(aseq, 15))
                     for aid, aseq in alleles.get(lid, {}).items()]

    rows = []
    per_locus = defaultdict(lambda: {"placed": 0, "conf": 0})
    for lid in kept:
        known = alleles.get(lid, {})
        for g, (seq, pl) in sorted(obs.get(lid, {}).items()):
            scored = match_allele(seq, kidx.get(lid, []))
            best_id, bi, bc = (scored[0][2], scored[0][0], scored[0][1]) \
                if scored else (".", 0.0, 0.0)
            sec_id, si = (scored[1][2], scored[1][0]) if len(scored) > 1 else (".", 0.0)
            marg = bi - si if len(scored) > 1 else 1.0
            exact = (bi >= 0.9999 and bc >= 0.9999
                     and len(seq) == len(known.get(best_id, "")))
            if pl["status"] == "ambiguous":
                st = "unresolved"
            elif exact:
                st = "exact_existing_allele"
            elif bi >= args.near_ident and bc >= args.near_cov:
                st = "near_existing_allele"
            else:
                st = "novel_outgroup_allele"
            # An outgroup that matches two alleles almost equally well has not
            # identified one; it must not cast a firm vote in Layer 3c.
            if st == "near_existing_allele" and marg < args.allele_margin:
                assign = "ambiguous_between_%s_%s" % tuple(sorted((best_id, sec_id)))
            elif st == "novel_outgroup_allele":
                assign = "NOVEL"
            elif st == "unresolved":
                assign = "."
            else:
                assign = best_id
            per_locus[lid]["placed"] += 1
            if st != "unresolved" and not str(assign).startswith("ambiguous_between"):
                per_locus[lid]["conf"] += 1
            rows.append([lid, g, len(seq), best_id, "%.4f" % bi, "%.4f" % bc,
                         sec_id, "%.4f" % si, "%.4f" % marg, assign,
                         pl["status"], "%.4f" % pl["margin"], st])
        if lid not in obs:
            rows.append([lid, ".", -1, ".", "NA", "NA", ".", "NA", "NA", ".",
                         "unplaced", "NA", "unresolved"])

    with open(os.path.join(args.outdir, "outgroup_alleles.tsv"), "w") as fh:
        fh.write("locus_id\toutgroup\tinterval_len\tbest_match_allele\t"
                 "allele_match_identity\tallele_match_coverage\t"
                 "second_best_allele\tsecond_best_identity\tallele_margin\t"
                 "allele_assignment\t"
                 "placement_status\tplacement_margin\tallele_match_status\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    call = Counter()
    for lid in kept:
        c = per_locus[lid]["conf"]
        call[">=2 confident outgroups" if c >= 2
             else "1 confident outgroup" if c == 1 else "0"] += 1
    st = Counter(r[12] for r in rows if r[1] != ".")
    amb_assign = sum(1 for r in rows
                     if r[1] != "." and str(r[9]).startswith("ambiguous_between"))
    n2 = call[">=2 confident outgroups"]
    with open(os.path.join(args.outdir, "callability.tsv"), "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in call.items():
            fh.write("%s\t%d\n" % (k, v))
        fh.write("loci_total\t%d\n" % len(kept))
        fh.write("outgroup_callability_pct\t%.1f\n"
                 % (100.0 * n2 / len(kept) if kept else 0))

    log("=" * 72)
    log("OUTGROUP ALLELE OBSERVATION  (no polarity assigned)")
    log("")
    log("  allele_match_status over %d placements:" % sum(st.values()))
    for k, v in st.most_common():
        log("    %-26s %5d" % (k, v))
    log("")
    log("  ambiguous allele assignments (excluded from confident votes): %d"
        % amb_assign)
    log("")
    log("  OUTGROUP CALLABILITY -- how much of the catalogue Layer 3 can reach:")
    for k in (">=2 confident outgroups", "1 confident outgroup", "0"):
        log("    %-26s %5d  %5.1f%%"
            % (k, call[k], 100.0 * call[k] / len(kept) if kept else 0))
    log("")
    log("  %.1f%% of loci have >=2 independent outgroups and are eligible for "
        "polarity" % (100.0 * n2 / len(kept) if kept else 0))
    log("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())

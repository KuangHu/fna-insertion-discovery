#!/usr/bin/env python3
"""ARM A -- within-genome multi-copy element discovery. Family-agnostic.

No transposase HMM, no IS library, no annotation. Discovery rests on one
structural signal only:

    copy1  AAAA---[====================]---CCCC
    copy2  GGGT---[====================]---TTAT
    copy3  CTAC---[====================]---AGGG
              ^     highly conserved     ^
           divergent                  divergent
            flank                       flank

An interval repeated at >=N distinct loci, with high internal identity and
ABRUPT loss of homology at both edges, is a mobile-element candidate. High
internal identity with homologous flanks is a segmental duplication instead,
and is scored down rather than discarded.

Boundaries are voted by the copies themselves: the element edge is where
cross-copy homology terminates, so alignment endpoints ARE the boundary
estimate, and >=3 copies vote on it.

Outputs (prefix given by --out):
  *_families.tsv   one row per candidate element  (+ S_mobility components)
  *_copies.tsv     one row per copy = one insertion-site context
  *_elements.fna   representative sequence per family
  *_sites.fna      per-copy flank context, element-oriented
  *_profiles.tsv   cross-copy conservation profile across each boundary
"""
import argparse
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import paf as pafmod                                  # noqa: E402
from lib.fasta import Fasta, revcomp, write_fasta               # noqa: E402
from lib.util import UnionFind, log, merge_intervals, overlap, run, seq_md5  # noqa: E402


# ---------------------------------------------------------------- self-align
def self_align(genome, out_paf, threads, minimap2="minimap2", k=19, w=19, m=200,
               f=0):
    """All-vs-all self alignment. -X = -DP --dual=no: skips the self diagonal
    and keeps every chain, so repeat copies are all reported.

    -f 0 DISABLES high-frequency minimizer masking, and it is not optional here.
    minimap2 masks the most frequent 0.02% of minimizers by default, which is
    right for read mapping and exactly backwards for this pipeline: the highest
    copy-number repeat in a genome is the one whose minimizers are most
    frequent, so the default silently deletes the very elements Arm A exists to
    find. Measured on GCA_011404755.1, which carries 60 IS1 copies: the default
    emitted 138 self-alignments and ZERO covering any IS1 locus; -f 0 emitted
    3,442, with 57 covering one locus. Across the 20-genome ISEScan benchmark
    this accounted for 113 of 502 multi-copy misses (22.5%), and every miss with
    ISEScan copy number >=27 was in that group.

    Cost is 0.66 s -> 4.3 s per 5 Mb genome, which is affordable.
    """
    cmd = [minimap2, "-X", "-c", "-k", str(k), "-w", str(w), "-m", str(m),
           "-N", "50", "-p", "0.05", "-f", str(f), "-t", str(threads),
           genome, genome]
    with open(out_paf, "w") as fh:
        run(cmd, stdout=fh, capture_output=False)
    return out_paf


# ------------------------------------------------------------ copy discovery
def collect_blocks(paf_path, min_id, min_len, max_len, min_sep):
    """Homology blocks between two DISTINCT loci."""
    blocks = []
    for r in pafmod.parse(paf_path):
        if r.identity < min_id:
            continue
        if not (min_len <= r.qspan <= max_len and min_len <= r.tspan <= max_len):
            continue
        # reject same-locus / tandem-adjacent hits: those are one copy, not two
        if r.qname == r.tname:
            if overlap(r.qs, r.qe, r.ts, r.te) > 0:
                continue
            if abs(r.qs - r.ts) < min_sep:
                continue
        blocks.append(r)
    return blocks


def canonical_copies(blocks, slack):
    """Merge every observed interval into canonical copy intervals per contig."""
    per_contig = defaultdict(list)
    for r in blocks:
        per_contig[r.qname].append((r.qs, r.qe))
        per_contig[r.tname].append((r.ts, r.te))
    copies, cid = {}, 0
    for contig, ivs in per_contig.items():
        for s, e in merge_intervals(ivs, slack=slack):
            copies[cid] = (contig, s, e)
            cid += 1
    # contig -> [(cid, s, e)] for lookup
    index = defaultdict(list)
    for c, (contig, s, e) in copies.items():
        index[contig].append((c, s, e))
    for v in index.values():
        v.sort(key=lambda x: x[1])
    return copies, index


def assign(index, contig, s, e, min_rec=0.5):
    """Canonical copy id with the largest reciprocal overlap."""
    best, best_ov = None, 0.0
    for cid, cs, ce in index.get(contig, ()):
        ov = overlap(s, e, cs, ce)
        if ov == 0:
            continue
        rec = min(ov / max(1, e - s), ov / max(1, ce - cs))
        if rec >= min_rec and ov > best_ov:
            best, best_ov = cid, ov
    return best


def build_families(blocks, index):
    """Union-find over canonical copies; carry the edges for orientation."""
    uf, edges = UnionFind(), []
    for r in blocks:
        a = assign(index, r.qname, r.qs, r.qe)
        b = assign(index, r.tname, r.ts, r.te)
        if a is None or b is None or a == b:
            continue
        uf.union(a, b)
        edges.append((a, b, r.strand, r.identity, r.qs, r.qe, r.ts, r.te))
    return uf.groups(), edges


def orient_family(members, edges):
    """BFS strand assignment so every copy is read in element orientation."""
    mset = set(members)
    adj = defaultdict(list)
    for a, b, strand, *_ in edges:
        if a in mset and b in mset:
            sign = 1 if strand == "+" else -1
            adj[a].append((b, sign))
            adj[b].append((a, sign))
    orient, conflict = {}, False
    seed = sorted(members)[0]
    orient[seed] = 1
    stack = [seed]
    while stack:
        u = stack.pop()
        for v, sign in adj[u]:
            want = orient[u] * sign
            if v not in orient:
                orient[v] = want
                stack.append(v)
            elif orient[v] != want:
                conflict = True
    for m in members:                      # unreachable copies default to +
        orient.setdefault(m, 1)
    return orient, conflict



# ------------------------------------------------- boundary refinement (vote)
def get_context(fa, copies, m, orient, pad):
    """Element-oriented context. Returns (seq, anchorL, anchorR)."""
    contig, s, e = copies[m]
    clen = fa.length(contig)
    lo, hi = max(0, s - pad), min(clen, e + pad)
    seq = fa.fetch(contig, lo, hi)
    if orient[m] == -1:
        return revcomp(seq), hi - e, hi - s
    return seq, s - lo, e - lo


def refine_boundaries(ctxs, rep_i, tmp, minimap2, threads, min_cov=0.5):
    """Vote the element boundary in ONE common frame: the representative copy.

    Every copy's context is aligned to the representative's context. Because
    the flanks are unrelated, each alignment TERMINATES at the element edge --
    so its endpoints in representative coordinates are that copy's vote.

    ctxs: list of (name, seq, aL, aR), element-oriented.
    Returns (refined_anchors, stats) where refined_anchors maps index ->
    (aL, aR) corrected into each copy's own coordinates.
    """
    rep_name, rep_seq, rL, rR = ctxs[rep_i]
    ref_fa, qry_fa, out_paf = tmp + ".rep.fa", tmp + ".qry.fa", tmp + ".ref.paf"
    with open(ref_fa, "w") as fh:
        write_fasta(fh, "REP", rep_seq)
    with open(qry_fa, "w") as fh:
        for i, (nm, sq, _, _) in enumerate(ctxs):
            if i != rep_i:
                write_fasta(fh, "Q%d" % i, sq)
    stats = {"votes_L": [], "votes_R": [], "n_votes": 0, "n_strand_conflict": 0}
    anchors = {rep_i: (rL, rR)}
    hits = {}
    if os.path.getsize(qry_fa) > 0:
        cmd = [minimap2, "-c", "-x", "asm20", "-k", "15", "-w", "10",
               "--secondary=no", "-t", str(threads), ref_fa, qry_fa]
        with open(out_paf, "w") as fh:
            run(cmd, stdout=fh, capture_output=False)
        span = max(1, rR - rL)
        for r in pafmod.parse(out_paf):
            i = int(r.qname[1:])
            # keep the alignment that covers most of the element in rep frame
            cov = overlap(r.ts, r.te, rL, rR) / span
            if cov < min_cov:
                continue
            if r.strand == "-":
                stats["n_strand_conflict"] += 1
                continue
            if i not in hits or cov > hits[i][0]:
                hits[i] = (cov, r)
    for i, (cov, r) in hits.items():
        stats["votes_L"].append(r.ts)
        stats["votes_R"].append(r.te)
    stats["n_votes"] = len(hits)
    if stats["votes_L"]:
        cL = statistics.median(stats["votes_L"])
        cR = statistics.median(stats["votes_R"])
        stats["disp_L"] = statistics.median([abs(v - cL) for v in stats["votes_L"]])
        stats["disp_R"] = statistics.median([abs(v - cR) for v in stats["votes_R"]])
    else:
        cL, cR = rL, rR
        stats["disp_L"] = stats["disp_R"] = float("nan")
    stats["cons_L"], stats["cons_R"] = cL, cR
    stats["elem_len_consensus"] = int(cR - cL)
    anchors[rep_i] = (int(cL), int(cR))
    # map the consensus rep-frame boundary back into each copy's own frame
    for i, (cov, r) in hits.items():
        anchors[i] = (int(r.qs + (cL - r.ts)), int(r.qs + (cR - r.ts)))
    for i in range(len(ctxs)):                      # unaligned copies keep nominal
        anchors.setdefault(i, (ctxs[i][2], ctxs[i][3]))
    return anchors, stats


# ---------------------------------------------------- conservation profiling
def kmers(seq, k):
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if "N" in km:
            continue
        out.add(min(km, revcomp(km)))
    return out


def jacc(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def median_pair_jaccard(seqs, k):
    """Median pairwise k-mer Jaccard across a set of window sequences."""
    ks = [kmers(s, k) for s in seqs if len(s) >= k]
    if len(ks) < 2:
        return float("nan")
    vals = [jacc(ks[i], ks[j]) for i in range(len(ks)) for j in range(i + 1, len(ks))]
    return statistics.median(vals) if vals else float("nan")


def profile(contexts, flank, inside, step, k):
    """Conservation across one boundary.

    contexts: list of element-oriented strings, all anchored so that index
    `flank` is the element start (edge='L') or the element end (edge='R').
    Returns [(offset_bp, median_jaccard, n_copies)]; offset < 0 = outside.
    """
    rows = []
    for off in range(-flank, inside, step):
        wins = []
        for ctx, anchor in contexts:
            a = anchor + off
            if a < 0 or a + step > len(ctx):
                continue
            wins.append(ctx[a:a + step])
        rows.append((off, median_pair_jaccard(wins, k), len(wins)))
    return rows


# ------------------------------------------------------------------- scoring
def norm(x, lo, hi):
    if x != x:                                     # NaN
        return 0.0
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def s_mobility(internal_id, flank_hom, sharpness, n_loci, w):
    comp = {
        "s_internal": norm(internal_id, 85.0, 99.5),
        "s_flank_divergence": 1.0 - norm(flank_hom, 0.0, 0.6),
        "s_boundary_sharpness": norm(sharpness, 0.0, 0.7),
        "s_multiplicity": norm(n_loci, 2.0, 8.0),
    }
    comp["S_mobility"] = sum(w[k] * comp[k] for k in w)
    return comp


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome", help="one assembly FASTA (uncompressed)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--sample-id", default=None, help="defaults to genome basename")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--paf", default=None, help="reuse an existing self-align PAF")
    # discovery thresholds
    ap.add_argument("--min-identity", type=float, default=90.0,
                    help="min internal identity %% between copies (default 90)")
    ap.add_argument("--min-len", type=int, default=300)
    ap.add_argument("--max-len", type=int, default=5000,
                    help="insertion elements only; >5 kb is out of scope")
    ap.add_argument("--min-copies", type=int, default=3,
                    help=">=3 distinct loci: 2 copies are often segmental dup")
    ap.add_argument("--min-locus-sep", type=int, default=10000,
                    help="min separation on the same contig to count as a distinct locus")
    ap.add_argument("--merge-slack", type=int, default=50)
    # flank / profile
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--profile-inside", type=int, default=300)
    ap.add_argument("--profile-step", type=int, default=50)
    ap.add_argument("--kmer", type=int, default=15)
    ap.add_argument("--max-copies-profiled", type=int, default=30)
    ap.add_argument("--min-sharpness", type=float, default=0.15,
                    help="min inside-minus-outside conservation drop (default 0.15)")
    ap.add_argument("--max-flank-homology", type=float, default=0.40,
                    help="flank Jaccard at/above this = segmental duplication")
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()

    sample = args.sample_id or os.path.basename(args.genome).split(".fna")[0]
    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)

    tmpdir = args.out + "_tmp"
    os.makedirs(tmpdir, exist_ok=True)

    fa = Fasta(args.genome, samtools=args.samtools)
    log("%s: %d contigs, %.2f Mb" % (sample, len(fa.names()), fa.total() / 1e6))

    paf_path = args.paf or (args.out + "_self.paf")
    if not args.paf or not os.path.exists(paf_path):
        log("self-aligning ...")
        self_align(args.genome, paf_path, args.threads, args.minimap2)

    blocks = collect_blocks(paf_path, args.min_identity, args.min_len,
                            args.max_len, args.min_locus_sep)
    log("%d homology blocks pass filters" % len(blocks))
    if not blocks:
        for suf in ("_families.tsv", "_copies.tsv", "_profiles.tsv"):
            open(args.out + suf, "w").close()
        open(args.out + "_elements.fna", "w").close()
        open(args.out + "_sites.fna", "w").close()
        log("no candidates")
        return

    copies, index = canonical_copies(blocks, args.merge_slack)
    families, edges = build_families(blocks, index)
    log("%d canonical copies -> %d raw families" % (len(copies), len(families)))

    # identity per copy pair, for internal-identity summary
    pair_id = defaultdict(list)
    for a, b, strand, ident, *_ in edges:
        pair_id[tuple(sorted((a, b)))].append(ident)

    fam_fh = open(args.out + "_families.tsv", "w")
    cp_fh = open(args.out + "_copies.tsv", "w")
    pr_fh = open(args.out + "_profiles.tsv", "w")
    el_fh = open(args.out + "_elements.fna", "w")
    st_fh = open(args.out + "_sites.fna", "w")

    fam_cols = ["family_id", "sample_id", "n_copies", "n_distinct_loci",
                "elem_len_median", "elem_len_mad", "elem_len_consensus",
                "boundary_disp_L_bp", "boundary_disp_R_bp", "n_boundary_votes",
                "internal_identity_median",
                "left_flank_jaccard", "right_flank_jaccard", "flank_homology",
                "boundary_sharpness_L", "boundary_sharpness_R", "boundary_sharpness",
                "s_internal", "s_flank_divergence", "s_boundary_sharpness",
                "s_multiplicity", "S_mobility", "verdict",
                "orientation_conflict", "n_copies_profiled", "rep_copy",
                "elem_md5", "flags"]
    cp_cols = ["family_id", "copy_id", "sample_id", "contig", "start", "end",
               "strand", "length", "contig_len", "dist_to_contig_start",
               "dist_to_contig_end", "both_flanks_in_contig", "n_frac_elem",
               "n_frac_flank", "boundary_refined", "left_flank_200",
               "right_flank_200", "qc_tier"]
    fam_fh.write("\t".join(fam_cols) + "\n")
    cp_fh.write("\t".join(cp_cols) + "\n")
    pr_fh.write("family_id\tedge\toffset_bp\tmedian_jaccard\tn_copies\n")

    weights = {"s_internal": 0.25, "s_flank_divergence": 0.30,
               "s_boundary_sharpness": 0.25, "s_multiplicity": 0.20}
    n_kept = 0

    for fi, members in enumerate(sorted(families, key=lambda m: -len(m))):
        members = sorted(members)
        if len(members) < args.min_copies:
            continue
        # distinct loci: different contig, or far apart on the same contig
        loci = []
        for m in members:
            contig, s, e = copies[m]
            if not any(c == contig and abs(s - s2) < args.min_locus_sep
                       for c, s2, _ in loci):
                loci.append((contig, s, e))
        if len(loci) < args.min_copies:
            continue

        fam_id = "%s.A%04d" % (sample, fi)
        orient, conflict = orient_family(members, edges)
        lens = [copies[m][2] - copies[m][1] for m in members]
        len_med = statistics.median(lens)
        len_mad = statistics.median([abs(x - len_med) for x in lens]) if len(lens) > 1 else 0.0

        idents = [statistics.median(v) for kk, v in pair_id.items()
                  if kk[0] in set(members) and kk[1] in set(members)]
        int_id = statistics.median(idents) if idents else float("nan")

        # element-oriented contexts, then let the copies VOTE the boundary
        prof_members = members[:args.max_copies_profiled]
        truncated = len(members) > len(prof_members)
        pad = args.flank * 2
        ctxs = [(m,) + get_context(fa, copies, m, orient, pad) for m in prof_members]
        ctxs = [(str(m), sq, aL, aR) for m, sq, aL, aR in ctxs]
        # representative = the copy whose length is closest to the median
        rep_i = min(range(len(ctxs)),
                    key=lambda i: abs((ctxs[i][3] - ctxs[i][2]) - len_med))
        anchors, bstat = refine_boundaries(
            ctxs, rep_i, os.path.join(tmpdir, fam_id), args.minimap2, args.threads)

        ctxL, ctxR, flanksL, flanksR = [], [], [], []
        for i, (nm, sq, _, _) in enumerate(ctxs):
            aL, aR = anchors[i]
            if not (0 <= aL < aR <= len(sq)):
                continue
            ctxL.append((sq, aL))
            ctxR.append((sq, aR))
            flanksL.append(sq[max(0, aL - args.flank):aL])
            flanksR.append(sq[aR:aR + args.flank])

        # refined element coordinates, mapped back to genome space
        refined = {}
        for i, m in enumerate(prof_members):
            aL, aR = anchors[i]
            contig, s0, e0 = copies[m]
            clen = fa.length(contig)
            lo, hi = max(0, s0 - pad), min(clen, e0 + pad)
            gs, ge = ((lo + aL, lo + aR) if orient[m] == 1 else (hi - aR, hi - aL))
            if 0 <= gs < ge <= clen:
                refined[m] = (gs, ge)

        jL = median_pair_jaccard([s for s in flanksL if s], args.kmer)
        jR = median_pair_jaccard([s for s in flanksR if s], args.kmer)
        flank_hom = statistics.mean([v for v in (jL, jR) if v == v]) \
            if any(v == v for v in (jL, jR)) else float("nan")

        rowsL = profile(ctxL, args.flank, args.profile_inside, args.profile_step, args.kmer)
        rowsR = profile([(s, a - args.profile_inside) for s, a in ctxR],
                        0, args.profile_inside + args.flank, args.profile_step, args.kmer)
        for off, j, n in rowsL:
            pr_fh.write("%s\tL\t%d\t%.4f\t%d\n" % (fam_id, off, j if j == j else -1, n))
        for off, j, n in rowsR:
            pr_fh.write("%s\tR\t%d\t%.4f\t%d\n" % (fam_id, off - args.profile_inside,
                                                   j if j == j else -1, n))

        def sharp(rows, inside_positive):
            """Inside-minus-outside conservation drop across the boundary."""
            ins = [j for o, j, _ in rows if j == j and (o >= 0 if inside_positive else o < 0)]
            out = [j for o, j, _ in rows if j == j and (o < 0 if inside_positive else o >= 0)]
            if not ins or not out:
                return float("nan")
            return max(0.0, statistics.median(ins) - statistics.median(out))

        shL = sharp(rowsL, True)
        shR = sharp([(o - args.profile_inside, j, n) for o, j, n in rowsR], False)
        sharpness = statistics.mean([v for v in (shL, shR) if v == v]) \
            if any(v == v for v in (shL, shR)) else float("nan")

        comp = s_mobility(int_id, flank_hom, sharpness, len(loci), weights)
        sh_ok = sharpness == sharpness and sharpness >= args.min_sharpness
        # a mobile element must show a real edge; rrn / phage repeats do not
        if flank_hom == flank_hom and flank_hom >= args.max_flank_homology:
            verdict = "SEGMENTAL_DUP_LIKE"
        elif not sh_ok and len(loci) >= 3:
            verdict = "REPEAT_NO_BOUNDARY"
        elif comp["S_mobility"] >= 0.70 and sharpness >= 0.30 and len(loci) >= 3:
            verdict = "STRONG_MOBILE_CANDIDATE"
        elif comp["S_mobility"] >= 0.50:
            verdict = "MOBILE_CANDIDATE"
        else:
            verdict = "WEAK"

        flags = []
        if len_med > args.max_len:
            flags.append("merged_span_exceeds_max")
        if truncated:
            flags.append("profiled_%d_of_%d_copies" % (len(prof_members), len(members)))
        if conflict:
            flags.append("orientation_conflict")

        rep = members[0]
        rcontig, rs, re_ = copies[rep]
        rep_seq = fa.fetch(rcontig, rs, re_)
        if orient[rep] == -1:
            rep_seq = revcomp(rep_seq)
        write_fasta(el_fh, "%s len=%d copies=%d S=%.3f" %
                    (fam_id, len(rep_seq), len(loci), comp["S_mobility"]), rep_seq)

        fam_fh.write("\t".join(map(str, [
            fam_id, sample, len(members), len(loci), int(len_med), int(len_mad),
            bstat["elem_len_consensus"],
            "%.1f" % bstat["disp_L"] if bstat["disp_L"] == bstat["disp_L"] else "NA",
            "%.1f" % bstat["disp_R"] if bstat["disp_R"] == bstat["disp_R"] else "NA",
            bstat["n_votes"],
            "%.2f" % int_id if int_id == int_id else "NA",
            "%.4f" % jL if jL == jL else "NA",
            "%.4f" % jR if jR == jR else "NA",
            "%.4f" % flank_hom if flank_hom == flank_hom else "NA",
            "%.4f" % shL if shL == shL else "NA",
            "%.4f" % shR if shR == shR else "NA",
            "%.4f" % sharpness if sharpness == sharpness else "NA",
            "%.4f" % comp["s_internal"], "%.4f" % comp["s_flank_divergence"],
            "%.4f" % comp["s_boundary_sharpness"], "%.4f" % comp["s_multiplicity"],
            "%.4f" % comp["S_mobility"], verdict, int(conflict), len(prof_members),
            "%s:%d-%d" % (rcontig, rs, re_), seq_md5(rep_seq),
            ";".join(flags) or ".",
        ])) + "\n")

        # per-copy site rows -- this is the multi-site bag for downstream modelling
        for m in members:
            contig, s, e = copies[m]
            if m in refined:
                s, e = refined[m]
            clen = fa.length(contig)
            elem = fa.fetch(contig, s, e)
            lf = fa.fetch(contig, max(0, s - args.flank), s)
            rf = fa.fetch(contig, e, min(clen, e + args.flank))
            if orient[m] == -1:
                elem, lf, rf = revcomp(elem), revcomp(rf), revcomp(lf)
            both_in = (s - args.flank >= 0) and (e + args.flank <= clen)
            nfe = elem.upper().count("N") / max(1, len(elem))
            fl = lf + rf
            nff = fl.upper().count("N") / max(1, len(fl))
            if both_in and nff < 0.01 and min(s, clen - e) >= args.flank:
                tier = "high"
            elif min(s, clen - e) >= 100:
                tier = "medium"
            else:
                tier = "unresolvable_contig_break"
            cp_fh.write("\t".join(map(str, [
                fam_id, "%s.c%d" % (fam_id, m), sample, contig, s, e,
                "+" if orient[m] == 1 else "-", e - s, clen, s, clen - e,
                int(both_in), "%.4f" % nfe, "%.4f" % nff, int(m in refined),
                lf[-200:], rf[:200], tier,
            ])) + "\n")
            write_fasta(st_fh, "%s.c%d %s:%d-%d strand=%s tier=%s" %
                        (fam_id, m, contig, s, e,
                         "+" if orient[m] == 1 else "-", tier),
                        lf + elem + rf)
        n_kept += 1

    for fh in (fam_fh, cp_fh, pr_fh, el_fh, st_fh):
        fh.close()
    if not args.keep_tmp:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
    log("%s: %d families kept (>=%d distinct loci)" % (sample, n_kept, args.min_copies))


if __name__ == "__main__":
    main()

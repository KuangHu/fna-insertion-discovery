#!/usr/bin/env python3
"""Scan Layer-2 loci for TWO SUBSTANTIAL ALLELES AT ONE ANCHORED SITE.

The question this answers is the project's founding one, asked without any
annotation: do independent elements arrive at the same target site? A locus
qualifies when its two longest alleles share an exact prefix and suffix but carry
different substantial middles -- two things sitting at one anchored position.

EVERY GUARD HERE EXISTS BECAUSE AN EARLIER AD-HOC VERSION GOT IT WRONG:

  1. ORIENTATION. IS elements insert either way. Every middle pair is scored
     against revcomp as well and the better orientation is taken. A forward-only
     comparison once reported ~0.50 identity for pairs that matched at 1.000
     reverse-complemented.

  2. THE NOISE FLOOR IS MEASURED, NOT CHOSEN. An 0.90 identity cutoff was
     uncalibrated and wrong. revcomp preserves base composition exactly, so the
     revcomp scores of these same middles are a composition-matched null. The
     floor is their median, and classification is relative to it.

  3. NO ALIGNED COLUMNS IS NOT 0% IDENTITY. When the aligner returns no aligned
     block the identity is UNDEFINED and written as "NA". Scoring it 0.0 once
     inflated a category by six loci.

  4. THE FLANK CONTROL MUST COME FROM THE GENOME, NEVER THE ALLELE. "Flank
     identity" read off a Layer-2 allele is a TAUTOLOGY: the flanks ARE the
     longest common prefix/suffix, so decomposition defines them as exactly
     matching and the answer is 100.00 on every locus, always. The control here
     extracts --flank bp from each carrier genome OUTSIDE the allele boundary.
     Validity limit: flank matching is an inclusion criterion, so this control is
     sound for the loci it returns and cannot estimate how often they occur.

  5. COVERAGE IS REPORTED AND THRESHOLDED. High identity over a partial
     alignment is not the same finding. --min-cov is applied uniformly, with no
     per-locus exceptions, even when it costs a large-gap locus.

  6. EFFECTIVE n IS GENOME PAIRS AND CONNECTED COMPONENTS, NOT LOCI. One
     divergent genome pair produces elevated middle divergence at every site it
     shares. The summary reports loci, distinct carrier pairs, and the connected
     components of the pair graph -- the last is the honest ceiling.

    96_shared_target_scan.py --recon sweep/*/recon --genome-dir <dir> \\
        --out shared_target/ecoli_sweep
"""
import argparse
import collections
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
NA = float("nan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_alleles_fna(path):
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


def load_genome(path):
    seqs, cur = {}, None
    for line in open(path):
        if line[0] == ">":
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--genome-dir", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-middle", type=int, default=300,
                    help="a middle shorter than this is not 'substantial'")
    ap.add_argument("--min-cov", type=float, default=0.75,
                    help="alignment coverage below this = not well-covered. "
                         "Applied uniformly, no per-locus exceptions.")
    ap.add_argument("--same-element", type=float, default=0.90,
                    help="identity at or above this reads as one element, "
                         "lengths differing -- consistent with a single arrival "
                         "plus a deletion, so it is the WEAKER class")
    ap.add_argument("--floor-margin", type=float, default=0.05,
                    help="a pair within this of the measured noise floor is "
                         "called unrelated")
    ap.add_argument("--offset-tol", type=int, default=20,
                    help="two long alleles count as being at the SAME site only "
                         "if their insertion offsets against the site's shortest "
                         "allele differ by at most this. Junction ambiguity is a "
                         "few bp; 150 bp apart is two different sites.")
    ap.add_argument("--contain-cov", type=float, default=0.90,
                    help="coverage of the shorter middle inside the longer, above "
                         "which the locus reads as NESTED (one element inside "
                         "another) rather than two arrivals at one site")
    ap.add_argument("--flank", type=int, default=2000)
    args = ap.parse_args()

    from Bio.Align import PairwiseAligner
    aln = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                          open_gap_score=-12, extend_gap_score=-0.4)

    def idcov(a, b):
        """Returns (identity, coverage) or (NA, NA) when there is NO alignment.
        NA is not 0.0 -- an undefined comparison must never score as maximally
        dissimilar."""
        if len(a) < 50 or len(b) < 50:
            return NA, NA
        try:
            al = aln.align(a, b)[0]
        except Exception:
            return NA, NA
        A, B = al.aligned
        cols = m = 0
        for (s0, s1), (l0, l1) in zip(A, B):
            m += sum(1 for x, y in zip(a[s0:s1], b[l0:l1]) if x == y)
            cols += s1 - s0
        if cols == 0:
            return NA, NA
        return m / cols, cols / min(len(a), len(b))

    # ---- gather loci -----------------------------------------------------
    gpath = {}
    for d in args.genome_dir:
        for f in os.listdir(d):
            if f.endswith(".fna"):
                gpath.setdefault(f[:-4], os.path.join(d, f))
    log("%d genome files visible" % len(gpath))

    cand = []
    n_loci = n_shortlong = 0
    for d in args.recon:
        tag = os.path.basename(os.path.dirname(os.path.abspath(d))) or "recon"
        aseq = read_alleles_fna(os.path.join(d, "alleles.fna"))
        carriers = collections.defaultdict(dict)
        ap_ = os.path.join(d, "alleles.tsv")
        if os.path.exists(ap_):
            for r in csv.DictReader(open(ap_), delimiter="\t"):
                carriers[r["locus_id"]][r["allele_id"]] = \
                    r.get("genomes", "").split(",")
        for lid, av in aseq.items():
            n_loci += 1
            if len(av) < 2:
                continue
            ordered = sorted(av.items(), key=lambda kv: -len(kv[1]))
            two = ordered[:2]
            (a1, s1), (a2, s2) = two
            # BOTH alleles must be long RELATIVE TO THE SITE'S SHORTEST allele.
            # Without this a 2-allele locus (short + long) is admitted as if it
            # were two elements at one site: its "two middles" are just the
            # insert and the residual of a plain insertion. Two loci entered an
            # earlier hand-run set exactly this way. For a 2-allele locus the
            # shortest IS s2, so s2 - shortest = 0 and the locus is excluded.
            shortest = len(ordered[-1][1])
            if (len(s1) - shortest < args.min_middle
                    or len(s2) - shortest < args.min_middle):
                n_shortlong += 1
                continue
            n = min(len(s1), len(s2))
            p = 0
            while p < n and s1[p] == s2[p]:
                p += 1
            q = 0
            while q < n - p and s1[len(s1) - 1 - q] == s2[len(s2) - 1 - q]:
                q += 1
            m1, m2 = s1[p:len(s1) - q], s2[p:len(s2) - q]
            if len(m1) < args.min_middle or len(m2) < args.min_middle:
                continue
            c1 = (carriers.get(lid, {}).get(a1) or [""])[0]
            c2 = (carriers.get(lid, {}).get(a2) or [""])[0]
            cand.append(dict(panel=tag, locus=lid, a1=a1, a2=a2, s1=s1, s2=s2,
                             m1=m1, m2=m2, lcp=p, lcs=q, g1=c1, g2=c2,
                             s0=ordered[-1][1], a0=ordered[-1][0]))
    log("%d Layer-2 loci scanned, %d have TWO alleles both >= %d bp longer than the\n           site's shortest allele (%d rejected as short-vs-long)"
        % (n_loci, len(cand), args.min_middle, n_shortlong))
    if not cand:
        log("no candidates")
        return 0

    # ---- orientation-agnostic identity + composition-matched null --------
    for c in cand:
        fi, fc = idcov(c["m1"], c["m2"])
        ri, rcv = idcov(c["m1"], rc(c["m2"]))
        c["fwd_id"], c["fwd_cov"] = fi, fc
        c["rev_id"], c["rev_cov"] = ri, rcv
        pick_rev = (ri == ri) and (fi != fi or ri > fi)
        c["orientation"] = "revcomp" if pick_rev else "forward"
        c["best_id"] = ri if pick_rev else fi
        c["best_cov"] = rcv if pick_rev else fc
    nullv = [c["rev_id"] for c in cand
             if c["orientation"] == "forward" and c["rev_id"] == c["rev_id"]]
    floor = statistics.median(nullv) if nullv else NA
    n_rev = sum(1 for c in cand if c["orientation"] == "revcomp")
    log("noise floor (median revcomp identity, composition-matched): %s  from %d"
        % ("%.3f" % floor if floor == floor else "NA", len(nullv)))
    log("pairs where revcomp is the better orientation: %d" % n_rev)

    for c in cand:
        b = c["best_id"]
        if b != b:
            c["klass"] = "no_alignment_NA"
        elif b >= args.same_element:
            c["klass"] = "same_element_len_differs"
        elif floor == floor and b <= floor + args.floor_margin:
            c["klass"] = "unrelated_at_noise_floor"
        else:
            c["klass"] = "diverged_homologue"
        c["well_covered"] = (c["best_cov"] == c["best_cov"]
                             and c["best_cov"] >= args.min_cov)

    # ---- NULL 1: SAME SITE? offsets against the site's shortest allele ----
    # Pairwise LCP/LCS between the two LONG alleles guarantees only that the
    # divergent middles lie between shared flanks. It does NOT guarantee they
    # occupy the same position in the target. Two elements inserted 150 bp apart
    # in one anchored interval produce exactly the same blocky signature and are
    # NOT a shared target site. Each long allele is therefore decomposed against
    # the site's SHORTEST allele, and the two insertion offsets are compared.
    #
    # NULL 2: NESTED? If one middle is contained in the other, a single element
    # with a second inserted into it explains the locus without two arrivals at
    # one site.
    def decompose_vs(long_s, short_s):
        n = min(len(long_s), len(short_s))
        a = 0
        while a < n and long_s[a] == short_s[a]:
            a += 1
        b = 0
        while b < n - a and long_s[len(long_s) - 1 - b] == short_s[len(short_s) - 1 - b]:
            b += 1
        return a, len(long_s) - a - b

    for c in cand:
        o1, i1 = decompose_vs(c["s1"], c["s0"])
        o2, i2 = decompose_vs(c["s2"], c["s0"])
        c["offset1"], c["offset2"] = o1, o2
        c["ins1"], c["ins2"] = i1, i2
        c["offset_delta"] = abs(o1 - o2)
        c["same_site"] = c["offset_delta"] <= args.offset_tol
        # containment: how much of the SHORTER middle aligns inside the longer
        sm, lm = (c["m1"], c["m2"]) if len(c["m1"]) <= len(c["m2"]) else (c["m2"], c["m1"])
        ci, cc = idcov(sm, lm)
        cir, ccr = idcov(sm, rc(lm))
        if (cir == cir) and (ci != ci or cir > ci):
            ci, cc = cir, ccr
        c["contain_id"], c["contain_cov"] = ci, cc
        c["nested"] = (ci == ci and cc == cc
                       and ci >= args.same_element and cc >= args.contain_cov)

    # ---- external flank control, from the GENOME not the allele ----------
    gcache = {}

    def locate(gs, seq):
        for cid, s in gs.items():
            i = s.find(seq)
            if i >= 0:
                return cid, i, i + len(seq), 1
            i = s.find(rc(seq))
            if i >= 0:
                return cid, i, i + len(seq), -1
        return None

    for c in cand:
        c["flank_id"] = c["flank_cov"] = NA
        c["flank_status"] = "."
        if not c["g1"] or not c["g2"]:
            c["flank_status"] = "no_carrier"
            continue
        ext = {}
        ok = True
        for g, allele in ((c["g1"], c["s1"]), (c["g2"], c["s2"])):
            if g not in gpath:
                ok = False
                c["flank_status"] = "genome_missing"
                break
            if g not in gcache:
                if len(gcache) > 40:
                    gcache.clear()
                gcache[g] = load_genome(gpath[g])
            loc = locate(gcache[g], allele)
            if loc is None:
                ok = False
                c["flank_status"] = "allele_not_located"
                break
            cid, st, en, strand = loc
            ctg = gcache[g][cid]
            up = ctg[max(0, st - args.flank):st]
            dn = ctg[en:en + args.flank]
            if strand == -1:
                up, dn = rc(dn), rc(up)
            ext[g] = (up, dn)
        if not ok:
            continue
        res = [idcov(ext[c["g1"]][k], ext[c["g2"]][k]) for k in (0, 1)]
        res = [r for r in res if r[0] == r[0]]
        if not res:
            c["flank_status"] = "flank_too_short"
            continue
        c["flank_id"] = statistics.mean(r[0] for r in res)
        c["flank_cov"] = statistics.mean(r[1] for r in res)
        c["flank_status"] = "ok"

    cols = ["panel", "locus", "a1", "a2", "g1", "g2", "len_m1", "len_m2",
            "lcp", "lcs", "fwd_id", "fwd_cov", "rev_id", "rev_cov",
            "orientation", "best_id", "best_cov", "well_covered",
            "flank_id", "flank_cov", "flank_status", "within_locus_gap",
            "offset1", "offset2", "offset_delta", "same_site",
            "ins1", "ins2", "contain_id", "contain_cov", "nested", "klass"]
    for c in cand:
        c["len_m1"], c["len_m2"] = len(c["m1"]), len(c["m2"])
        c["within_locus_gap"] = (100 * (c["flank_id"] - c["best_id"])
                                 if c["flank_id"] == c["flank_id"]
                                 and c["best_id"] == c["best_id"] else NA)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + "_shared_target.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in sorted(cand, key=lambda x: (x["panel"], x["locus"])):
            fh.write("\t".join(
                ("%.4f" % c[k]) if isinstance(c.get(k), float) and c[k] == c[k]
                else ("NA" if isinstance(c.get(k), float) else str(c.get(k, ".")))
                for k in cols) + "\n")

    # ---- effective n: pairs, then connected components -------------------
    def report(sel, title):
        pairs = collections.Counter()
        for c in sel:
            if c["g1"] and c["g2"]:
                pairs[tuple(sorted((c["g1"], c["g2"])))] += 1
        par = {}

        def find(x):
            par.setdefault(x, x)
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for a, b in pairs:
            par.setdefault(a, a)
            par.setdefault(b, b)
            par[find(a)] = find(b)
        comp = collections.defaultdict(lambda: {"g": set(), "p": 0, "l": 0})
        for (a, b), nl in pairs.items():
            k = comp[find(a)]
            k["g"] |= {a, b}
            k["p"] += 1
            k["l"] += nl
        log("  %-34s %4d loci  %3d pairs  %3d genomes  %2d COMPONENTS"
            % (title, len(sel), len(pairs), len({g for p in pairs for g in p}),
               len(comp)))
        return comp

    log("=" * 96)
    log("SHARED TARGET SITE SCAN")
    log("")
    h = collections.Counter(c["klass"] for c in cand)
    log("  classification against the MEASURED floor (%s):"
        % ("%.3f" % floor if floor == floor else "NA"))
    for k in ("diverged_homologue", "same_element_len_differs",
              "unrelated_at_noise_floor", "no_alignment_NA"):
        log("    %-30s %4d" % (k, h[k]))
    log("")
    log("  well-covered means best_cov >= %.2f, applied with no exceptions"
        % args.min_cov)
    log("")
    log("  NULL 1 -- same site? (insertion offsets vs the site's shortest allele,")
    log("           tolerance %d bp):" % args.offset_tol)
    ss = sum(1 for c in cand if c["same_site"])
    log("    same site      %4d" % ss)
    log("    DIFFERENT site %4d   <- blocky signature, but two separate positions"
        % (len(cand) - ss))
    log("  NULL 2 -- nested? (shorter middle contained in the longer at >= %.2f "
        "identity and >= %.2f coverage):" % (args.same_element, args.contain_cov))
    log("    nested         %4d" % sum(1 for c in cand if c["nested"]))
    log("")
    dh = [c for c in cand if c["klass"] == "diverged_homologue"]
    un = [c for c in cand if c["klass"] == "unrelated_at_noise_floor"]
    report(cand, "ALL candidates")
    report(dh, "diverged homologue (all)")
    dhw = [c for c in dh if c["well_covered"]]
    report(dhw, "diverged homologue, WELL-COVERED")
    surv = [c for c in dhw if c["same_site"] and not c["nested"]]
    comp = report(surv, "^ AND same-site AND not nested")
    report([c for c in un if c["well_covered"]], "unrelated, WELL-COVERED")
    log("")
    if comp:
        log("  components of the well-covered diverged-homologue pair graph:")
        for v in sorted(comp.values(), key=lambda x: -x["l"]):
            log("    %d genomes, %d pairs, %d loci" % (len(v["g"]), v["p"], v["l"]))
    gaps = [c["within_locus_gap"] for c in surv
            if c["within_locus_gap"] == c["within_locus_gap"]]
    if gaps:
        log("")
        log("  WITHIN-LOCUS GAP (external flank identity minus middle identity),")
        log("  flanks taken from the carrier GENOMES, not from the alleles:")
        log("    median %.1f points   range %.1f to %.1f   n=%d"
            % (statistics.median(gaps), min(gaps), max(gaps), len(gaps)))
        log("    A gap this size between genomes whose local flanking sequence")
        log("    matches cannot arise by vertical inheritance in place.")
    sf = collections.Counter(c["flank_status"] for c in cand)
    log("")
    log("  flank control status: " + "  ".join("%s=%d" % kv for kv in sf.most_common()))
    log("=" * 96)
    log("wrote %s_shared_target.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

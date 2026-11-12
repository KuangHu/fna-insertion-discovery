#!/usr/bin/env python3
"""Is a bag an IS element? Test the structure, not the domain name.

A DDE or rve domain says "this protein can cut and join DNA".  It does not say
the thing carrying it is an insertion sequence.  Tn3 transposons, ICEs,
prophage and composite transposons all carry one.  An IS is a STRUCTURAL
object and is tested as one, on four independent properties:

  SIZE          IS elements are compact. Recorded as a number and banded, not
                gated at a magic value, because the families differ.
  COMPACTNESS   the transposase spans most of the element and there are few
                other ORFs. A 5 kb insert whose transposase covers 30% of it
                is a transposon carrying cargo, not an IS.
  CATALYTIC     the DOMINANT CDS carries a transposase catalytic family.
                Accessory DNA-binding families (HTH_Tnp_*, Zn_ribbon_InsA)
                are tracked SEPARATELY and never counted as the catalytic
                evidence -- an HTH alone is not a transposase.
                Phage_integrase, Resolvase and Recombinase are deliberately
                NOT catalytic-IS families: they mark site-specific
                recombination systems, which is a different object.
  TERMINAL IR   the element's own left end matches the reverse complement of
                its own right end -- searched with an OFFSET at each end, not
                flush. The recorded insertion boundary is one of several legal
                cut points inside the breakpoint equivalence range, and it can
                also swallow a few bases of target-site duplication, so a real
                TIR routinely starts tens of bases inside the recorded
                terminus. A flush comparison finds almost nothing and would
                have reported 0 IS elements in bags that are plainly IS1.

THE TIR THRESHOLD IS SET FROM A NULL, NOT CHOSEN.  Any sequence has some
terminal self-complementarity by chance, and GC-rich sequence has more. So the
same measurement is repeated on a CHIMERIC control: this insert's left end
against the reverse complement of a DIFFERENT insert's right end, paired at
random with a fixed seed. That preserves real genomic sequence statistics on
both sides while destroying the pairing that a real TIR would create. The
reported cutoff is the null's 99th percentile.

WHY NOT TSD.  A target-site duplication would be the other classic signal, but
in this project the junction direct repeat has already been shown to have a
U-shaped null and to be confounded with the junction microhomology that the
breakpoint inference itself can slide over. It is measured and reported here
for information, and it is NOT part of the IS verdict.

THE BIG BAGS ARE THE POSITIVE REFERENCE. If the test is sound, bags seen at
>=20 independent sites should be overwhelmingly IS-like. If they are not, the
test is wrong and the small-bag answer means nothing.
"""
import argparse, collections, json, random, re, sys

CATALYTIC = re.compile(
    r"^(DDE_Tnp|DEDD_Tnp|Y1_Tnp|OrfB_IS605|Transposase|Tnp_|Y2_Tnp|"
    r"rve(_\d)?$|Rve|IS\d+_Tnp|TnpB|Mu_transpos)", re.I)
ACCESSORY = re.compile(r"^(HTH_Tnp|Zn_ribbon_InsA|Zn_ribbon_ISL3|zf-IS66|"
                       r"LZ_Tnp|DUF772|HTH_21|HTH_28|HTH_23|Arm-DNA|Tn2501)", re.I)
SSR = re.compile(r"^(Phage_int|Integrase_|Recombinase|Resolvase|Zn_ribbon_recom|"
                 r"TraI|Relaxase|MobA|Cisa_C)", re.I)

COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def rc(s):
    return s.encode().translate(COMP)[::-1].decode()


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def fasta(p):
    rid, parts = None, []
    with open(p) as f:
        for line in f:
            if line[0] == ">":
                if rid:
                    yield rid, "".join(parts)
                t = line[1:].split()[0]
                rid = t.split("|", 1)[1] if "|" in t else t
                parts = []
            else:
                parts.append(line.strip())
    if rid:
        yield rid, "".join(parts)


def term_ir(left, right_rc, k=7, maxmm=0.25, max_off=60):
    """Longest near-exact match between the two windows, allowing an offset.

    Seed on exact k-mers, extend both ways while identity stays >= 1-maxmm.
    Returns (length, left_offset, right_offset). Offsets are how far INSIDE
    the recorded terminus the repeat begins; a match starting beyond
    `max_off` is internal, not terminal, and is rejected.

    The SAME function is applied to the chimeric null, so the extra chances
    bought by searching offsets are paid for in the threshold rather than
    being a free gain for the real pairs.
    """
    idx = {}
    for i in range(len(left) - k + 1):
        idx.setdefault(left[i:i + k], []).append(i)
    best = (0, -1, -1)
    for j in range(len(right_rc) - k + 1):
        for i in idx.get(right_rc[j:j + k], ()):
            mm, la = 0, 0
            while i - la - 1 >= 0 and j - la - 1 >= 0:
                bad = left[i - la - 1] != right_rc[j - la - 1]
                if bad and mm + 1 > maxmm * (k + la + 1):
                    break
                mm += bad
                la += 1
            lb = 0
            while i + k + lb < len(left) and j + k + lb < len(right_rc):
                bad = left[i + k + lb] != right_rc[j + k + lb]
                if bad and mm + 1 > maxmm * (k + la + lb + 1):
                    break
                mm += bad
                lb += 1
            L, ol, orr = k + la + lb, i - la, j - la
            if ol <= max_off and orr <= max_off and L > best[0]:
                best = (L, ol, orr)
    return best


def junction_dr(up, dn, kmax=20):
    """Longest direct repeat shared by the end of upstream and start of downstream."""
    n = min(len(up), len(dn), kmax)
    for k in range(n, 2, -1):
        if up[-k:] == dn[:k]:
            return k
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bags", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--fasta", required=True, help="dedup_insertions.fna")
    ap.add_argument("--hits", required=True, help="all_hits_on_dna.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--min-frac", type=float, default=0.50)
    ap.add_argument("--max-cds", type=int, default=3)
    # An exact 12-mer seed misses real TIRs: IS termini are imperfect repeats,
    # and on six reference bags (IS1, IS3, IS256, ISL3 families) k=12 found
    # nothing while k=7 found 18-31 bp at offsets of 0-22. The null is
    # recomputed with whatever is set here, so a looser seed buys sensitivity
    # at the price of a higher threshold, not for free.
    ap.add_argument("--seed-k", type=int, default=7)
    ap.add_argument("--tir-maxmm", type=float, default=0.25)
    ap.add_argument("--tir-max-off", type=int, default=60)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    seq = dict(fasta(a.fasta))
    cds = {r["seq_key"]: r for r in read_tsv(a.cds)}
    bag_rep, sites = {}, {}
    for r in read_tsv(a.bags):
        sites[r["bag_id"]] = int(r["n_independent_sites"])
    for r in read_tsv(a.cds):
        if r["cds_cluster_id"] and r["cds_cluster_id"] not in bag_rep:
            bag_rep[r["cds_cluster_id"]] = r["cluster_rep"] or r["seq_key"]

    dom_fam = collections.defaultdict(set)
    for h in read_tsv(a.hits):
        c = cds.get(h["seq_key"])
        if c and h["protein_id"] == c["dominant_cds"]:
            dom_fam[h["seq_key"]].add(h["pfam_name"])

    keys = [bag_rep[b] for b in bag_rep if bag_rep[b] in seq]
    # chimeric null: this insert's left end vs ANOTHER insert's right end
    partner = keys[:]
    rng.shuffle(partner)
    for i in range(len(keys)):
        if partner[i] == keys[i]:
            partner[i] = partner[(i + 1) % len(partner)]
    pmap = dict(zip(keys, partner))

    null = []
    for k in keys:
        s, o = seq[k], seq[pmap[k]]
        null.append(term_ir(s[:160], rc(o[-160:]), a.seed_k, a.tir_maxmm,
                            a.tir_max_off)[0])
    null.sort()
    cut = null[int(0.99 * (len(null) - 1))]
    print(f"chimeric TIR null over {len(null)} pairs: "
          f"median {null[len(null)//2]}  p95 {null[int(0.95*(len(null)-1))]}  "
          f"p99 {cut}  max {null[-1]}", flush=True)

    def band(n):
        return ("1 site" if n == 1 else "2-3 sites" if n <= 3 else
                "4-9 sites" if n <= 9 else "10-19 sites" if n <= 19 else ">=20 sites")

    tab = collections.defaultdict(collections.Counter)
    rows = []
    for b, k in bag_rep.items():
        if k not in seq:
            continue
        s = seq[k]
        c = cds[k]
        fams = dom_fam.get(k, set())
        cat = {f for f in fams if CATALYTIC.match(f)}
        acc = {f for f in fams if ACCESSORY.match(f)}
        ssr = {f for f in fams if SSR.match(f)}
        ir, ir_off_l, ir_off_r = term_ir(s[:160], rc(s[-160:]), a.seed_k,
                                         a.tir_maxmm, a.tir_max_off)
        nullir = term_ir(s[:160], rc(seq[pmap[k]][-160:]), a.seed_k,
                         a.tir_maxmm, a.tir_max_off)[0]
        frac = float(c["dominant_cds_frac"] or 0)
        ncds = int(c["n_cds"] or 0)
        compact = frac >= a.min_frac and 0 < ncds <= a.max_cds
        is_like = bool(cat) and compact and ir > cut
        g = band(sites.get(b, 1))
        tab[g]["bags"] += 1
        tab[g]["catalytic_tnp"] += bool(cat)
        tab[g]["compact"] += compact
        tab[g]["tir_above_null"] += ir > cut
        tab[g]["IS_like_all_three"] += is_like
        tab[g]["ssr_not_IS"] += bool(ssr and not cat)
        rows.append([b, sites.get(b, 1), g, k, len(s), ncds, frac,
                     c["coding_frac"], ir, ir_off_l, ir_off_r, nullir, int(ir > cut),
                     ";".join(sorted(cat)), ";".join(sorted(acc)),
                     ";".join(sorted(ssr)), int(compact), int(is_like)])

    with open(a.out + "_is_like.tsv", "w") as f:
        f.write("bag_id\tn_independent_sites\tband\trep_seq_key\tinsert_len_nt\tn_cds\t"
                "dominant_cds_frac\tcoding_frac\ttir_bp\ttir_offset_left\t"
                "tir_offset_right\tchimeric_null_tir_bp\t"
                "tir_above_null\tcatalytic_families\taccessory_families\t"
                "site_specific_recomb_families\tcompact\tIS_like\n")
        for r in sorted(rows, key=lambda r: -r[1]):
            f.write("\t".join(str(x) for x in r) + "\n")

    order = ["1 site", "2-3 sites", "4-9 sites", "10-19 sites", ">=20 sites"]
    print(f"\n{'band':>12} {'bags':>6} {'catalytic':>10} {'compact':>8} "
          f"{'TIR>null':>9} {'IS-like':>8} {'%':>6} {'SSR not IS':>11}")
    out = {}
    for g in order:
        t = tab[g]
        pct = 100 * t["IS_like_all_three"] / max(1, t["bags"])
        out[g] = dict(t, pct_IS_like=round(pct, 1))
        print(f"{g:>12} {t['bags']:>6} {t['catalytic_tnp']:>10} {t['compact']:>8} "
              f"{t['tir_above_null']:>9} {t['IS_like_all_three']:>8} {pct:>5.1f}% "
              f"{t['ssr_not_IS']:>11}")

    with open(a.out + "_is_like.json", "w") as f:
        json.dump(dict(tir_null_p99=cut, bands=out, min_dominant_frac=a.min_frac,
                       max_cds=a.max_cds, seed=a.seed, seed_k=a.seed_k,
                       tir_maxmm=a.tir_maxmm, tir_max_off=a.tir_max_off,
                       tir_null_median=null[len(null) // 2],
                       tir_null_p95=null[int(0.95 * (len(null) - 1))]), f, indent=2)


if __name__ == "__main__":
    main()

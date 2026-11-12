#!/usr/bin/env python3
"""Why is a one-site bag only seen once? Compare its dominant CDS to the common bags.

A bag seen at a single independent site has three very different explanations,
and the Pfam family alone cannot tell them apart:

  it is a DECAYED OR PARTIAL COPY of a common element -- same transposase,
      shorter protein, so it misses the 80% coverage gate and clusters alone;
  it is a DIVERGENT MEMBER of a common family, below the 50% identity gate;
  it is genuinely UNRELATED to anything seen often, which for a mobility
      domain means a low-copy element class -- integrative, conjugative,
      Tn3-like -- rather than a new IS.

So each one-site bag's dominant protein is searched against the dominant
proteins of the bags seen at >= --big-sites independent sites, at high
sensitivity and WITHOUT a coverage requirement, because the coverage gate is
precisely the thing under test.  A gate cannot be used to audit itself.

THE VERDICTS, and what each rules out:
  SHOULD_HAVE_CLUSTERED            >=50% id AND >=80% coverage both ways.
                                   This must be zero; a nonzero count means
                                   the clustering is internally inconsistent
                                   and nothing downstream can be trusted. It
                                   is therefore reported as a hard check.
  same_family_shorter_protein      >=50% id, query < 0.8x the target length
  same_family_longer_protein       >=50% id, query > 1.2x -- fusion/composite
  same_family_similar_length_low_cov  >=50% id, length matches, alignment does not
  divergent_homolog                30-50% identity
  weak_similarity                  <30% identity
  no_homolog_in_any_big_bag        nothing at e <= --evalue

LENGTH IS COMPARED, NOT ASSUMED.  "Shorter" is measured against the matched
representative's length, not against some absolute idea of how long a
transposase should be, because the families here differ in length by an order
of magnitude.
"""
import argparse, collections, json, os, subprocess, sys


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def subset_faa(src, ids, dst):
    keep, n = set(ids), 0
    with open(src) as f, open(dst, "w") as o:
        emit = False
        for line in f:
            if line[0] == ">":
                emit = line[1:].split()[0] in keep
                n += emit
            if emit:
                o.write(line)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bags", required=True, help="<tag>_bags.tsv")
    ap.add_argument("--cds", required=True, help="insertion_cds_<tag>.tsv")
    ap.add_argument("--faa", required=True, help="dominant_<tag>.faa")
    ap.add_argument("--query-bags", required=True,
                    help="<tag>_singleton_mge.tsv -- the one-site bags to explain")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--big-sites", type=int, default=20)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--sens", type=float, default=7.5)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.tmp, exist_ok=True)

    qrows = list(read_tsv(a.query_bags))
    qids = [r["seq_key"] for r in qrows if r["seq_key"]]
    big = {r["bag_id"] for r in read_tsv(a.bags)
           if int(r["n_independent_sites"]) >= a.big_sites}
    reps = sorted({r["cluster_rep"] for r in read_tsv(a.cds)
                   if r["cds_cluster_id"] in big and r["cluster_rep"]})
    # a query must never be its own target, or every bag explains itself
    reps = [r for r in reps if r not in set(qids)]

    qfa, tfa = os.path.join(a.tmp, "q.faa"), os.path.join(a.tmp, "t.faa")
    nq = subset_faa(a.faa, qids, qfa)
    nt = subset_faa(a.faa, reps, tfa)
    print(f"queries {nq}/{len(qids)}   targets {nt}/{len(reps)} "
          f"(bags with >= {a.big_sites} sites)", flush=True)
    if nq == 0 or nt == 0:
        sys.exit("FATAL: empty query or target set")

    m8 = a.out + "_vs_bigbag.m8"
    cmd = ["mmseqs", "easy-search", qfa, tfa, m8, os.path.join(a.tmp, "mm"),
           "--threads", str(a.threads), "-s", str(a.sens), "-e", str(a.evalue),
           "--max-seqs", "300", "-v", "1", "--format-output",
           "query,target,fident,alnlen,qlen,tlen,qcov,tcov,evalue,bits"]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    best = {}
    with open(m8) as f:
        for ln in f:
            q, t, fid, aln, ql, tl, qc, tc, e, b = ln.split()
            r = dict(t=t, fid=float(fid), qc=float(qc), tc=float(tc),
                     ql=int(ql), tl=int(tl), e=e, b=float(b))
            if q not in best or r["b"] > best[q]["b"]:
                best[q] = r

    cnt = collections.Counter()
    hdr = list(qrows[0].keys())
    with open(a.out + "_mge_classified.tsv", "w") as o:
        o.write("\t".join(hdr) + "\tverdict\tbest_bigbag_rep\tbest_identity\t"
                "query_cov\ttarget_cov\tquery_aa\ttarget_aa\tevalue\n")
        for r in qrows:
            v = best.get(r["seq_key"])
            if v is None:
                verdict, extra = "no_homolog_in_any_big_bag", [""] * 7
            else:
                if v["fid"] >= 0.50 and v["qc"] >= 0.80 and v["tc"] >= 0.80:
                    verdict = "SHOULD_HAVE_CLUSTERED"
                elif v["fid"] >= 0.50:
                    verdict = ("same_family_shorter_protein" if v["ql"] < v["tl"] * 0.8
                               else "same_family_longer_protein" if v["ql"] > v["tl"] * 1.2
                               else "same_family_similar_length_low_cov")
                elif v["fid"] >= 0.30:
                    verdict = "divergent_homolog"
                else:
                    verdict = "weak_similarity"
                extra = [v["t"], f"{v['fid']:.3f}", f"{v['qc']:.3f}", f"{v['tc']:.3f}",
                         str(v["ql"]), str(v["tl"]), v["e"]]
            cnt[verdict] += 1
            o.write("\t".join([r[c] for c in hdr] + [verdict] + extra) + "\n")

    fams = collections.defaultdict(collections.Counter)
    for r in qrows:
        v = best.get(r["seq_key"])
        verdict = ("no_homolog_in_any_big_bag" if v is None else
                   "has_homolog")
        for x in (r.get("dominant_family_names") or "").split(";"):
            if x:
                fams[verdict][x] += 1
    with open(a.out + "_mge_families_by_verdict.tsv", "w") as o:
        o.write("verdict\tpfam_id\tn_bags\n")
        for k, c in fams.items():
            for nm, n in c.most_common():
                o.write(f"{k}\t{nm}\t{n}\n")

    ids = sorted(v["fid"] for v in best.values())
    summ = dict(queries=len(qrows), targets=nt, big_sites=a.big_sites,
                evalue=a.evalue, sensitivity=a.sens, verdicts=dict(cnt),
                median_best_identity_where_hit=round(ids[len(ids) // 2], 3) if ids else None)
    with open(a.out + "_singleton_vs_bigbags.json", "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))
    if cnt["SHOULD_HAVE_CLUSTERED"]:
        sys.exit(f"FATAL: {cnt['SHOULD_HAVE_CLUSTERED']} one-site bags pass the "
                 f"clustering gate against a big bag -- the clustering is inconsistent")


if __name__ == "__main__":
    main()

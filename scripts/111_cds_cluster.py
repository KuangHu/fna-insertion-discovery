#!/usr/bin/env python3
"""Predict CDS in every insert, cluster the proteins, carve the non-coding remainder.

This is the bridge to the Channel A/B bag format (docs/CANONICAL_BAG_SPEC.md in
the DL project). Three products, one pass:

  cds_cluster_id       mmseqs cluster of the insert's DOMINANT protein.

                       **Deliberately NOT named transposase_id: nothing in this
                       step identifies a transposase.** prodigal predicts ORFs ab
                       initio and mmseqs groups them by sequence. No HMM, no IS
                       library, no functional claim. A cluster is a cluster.

                       MULTI-CDS RULE: when an insert carries several ORFs, the
                       cluster of the LONGEST ORF defines its cds_cluster_id.
                       Rationale: where an element has a transposase it is
                       normally the longest ORF, and passenger genes
                       (resistance, toxin-antitoxin) are shorter. This is a
                       HEURISTIC: `dominant_orf_frac` reports how much of the
                       insert the chosen ORF covers so a caller can see where it
                       is weak, and `n_orfs` reports how often the question
                       arises at all.

                       A PROTEIN cluster, not a nucleotide one: elements of one
                       family diverge far more in DNA than in protein.
  noncoding_regions    the insert MINUS its CDS spans. In the bag spec this is
                       the ncRNA context, and it is a DEPLOY field.
  cds inventory        how many ORFs, and how much of the insert is coding.

WHY PROTEIN AND NOT NUCLEOTIDE. `106_element_families.py` already clusters the
inserts by DNA. That answers "is this the same element"; it does NOT answer "is
this the same transposase family", which is what `cds_cluster_id` must mean.
Both are kept -- they are different groupings and are expected to disagree.

DOMINANT ORF, not all ORFs. An insert may carry passenger genes (resistance,
toxin-antitoxin). The transposase is usually the longest ORF and usually spans
most of the element, so the dominant ORF defines identity and the rest are
recorded but unused. That is a HEURISTIC: `dominant_orf_frac` reports how much
of the insert it covers, so a caller can see where it is weak.

NO ANNOTATION IS USED. prodigal is ab initio and mmseqs is sequence-only.
Nothing is matched against a transposase HMM or an IS library, so a cluster is
a cluster, not a named family. Naming is Layer 5.

    111_cds_cluster.py --inserts database/insertions_inserts.fna \\
        --db database/insertions.tsv --out cds/
"""
import argparse
import collections
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def read_fa(p):
    out, cur = {}, None
    for line in open(p):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-seq-id", type=float, default=0.50,
                    help="mmseqs protein identity for a transposase cluster. 0.5 "
                         "is the usual family-level cut for transposases; it is "
                         "recorded so the choice is visible, not defended as "
                         "optimal.")
    ap.add_argument("--cov", type=float, default=0.80)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    faa = os.path.join(args.out, "orfs.faa")
    gff = os.path.join(args.out, "orfs.gff")

    if not (os.path.exists(faa) and os.path.getsize(faa) > 0):
        log("prodigal -p meta over the inserts ...")
        rc = subprocess.run(["prodigal", "-i", args.inserts, "-a", faa,
                             "-f", "gff", "-o", gff, "-p", "meta", "-q"],
                            check=False).returncode
        if rc != 0 or not os.path.exists(faa) or os.path.getsize(faa) == 0:
            log("FATAL: prodigal produced nothing (rc=%d)" % rc)
            return 2

    n_orf = sum(1 for l in open(faa) if l.startswith(">"))
    n_ins = sum(1 for l in open(args.inserts) if l.startswith(">"))
    log("%d ORFs predicted across %d inserts (%.2f per insert)"
        % (n_orf, n_ins, n_orf / n_ins if n_ins else 0))
    if n_orf == 0:
        log("FATAL: no ORFs predicted")
        return 2

    # prodigal names an ORF "<record>_<n>" and writes "# start # end # strand"
    orfs = collections.defaultdict(list)
    for line in open(faa):
        if not line.startswith(">"):
            continue
        p = line[1:].strip().split("#")
        name = p[0].strip()
        try:
            st, en = int(p[1]), int(p[2])
        except (IndexError, ValueError):
            continue
        rec = name.rsplit("_", 1)[0]
        orfs[rec].append((name, st, en, en - st + 1))

    ins_len = {k: len(v) for k, v in read_fa(args.inserts).items()}
    dom = {}
    for rec, lst in orfs.items():
        lst.sort(key=lambda x: -x[3])
        dom[rec] = lst[0]
    log("%d of %d inserts carry >=1 ORF (%.1f%%)"
        % (len(dom), n_ins, 100.0 * len(dom) / n_ins if n_ins else 0))
    if n_ins and len(dom) < 0.5 * n_ins:
        log("NOTE: fewer than half the inserts contain a predicted ORF, so "
            "cds_cluster_id will be absent for the rest")

    keep = {v[0] for v in dom.values()}
    domfaa = os.path.join(args.out, "dominant.faa")
    with open(domfaa, "w") as fh:
        w = False
        for line in open(faa):
            if line.startswith(">"):
                w = line[1:].strip().split("#")[0].strip() in keep
            if w:
                fh.write(line)

    log("clustering %d dominant ORFs: mmseqs id>=%.2f cov>=%.2f"
        % (len(keep), args.min_seq_id, args.cov))
    pre = os.path.join(args.out, "clu")
    rc = subprocess.run(["mmseqs", "easy-cluster", domfaa, pre,
                         os.path.join(args.out, "mmtmp"),
                         "--min-seq-id", str(args.min_seq_id),
                         "-c", str(args.cov), "--cov-mode", "0",
                         "--threads", str(args.threads), "-v", "1"],
                        check=False).returncode
    tsv = pre + "_cluster.tsv"
    if rc != 0 or not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        log("FATAL: mmseqs produced no cluster table (rc=%d)" % rc)
        return 2

    rep = {}
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            rep[f[1]] = f[0]
    log("mmseqs: %d members in %d clusters" % (len(rep), len(set(rep.values()))))
    # standing practice 6: report the join rate before anything derived from it
    matched = sum(1 for k in keep if k in rep)
    log("dominant-ORF -> cluster join: %d/%d (%.1f%%)"
        % (matched, len(keep), 100.0 * matched / len(keep) if keep else 0))
    if keep and matched < 0.5 * len(keep):
        log("FATAL: fewer than half the dominant ORFs appear in the cluster "
            "table -- a name mismatch between the FASTA and mmseqs output")
        return 2

    cid = {}
    for i, r in enumerate(sorted(set(rep.values())), 1):
        cid[r] = "CDS%05d" % i

    db = {r["db_id"]: r for r in csv.DictReader(open(args.db), delimiter="\t")}
    cols = ["db_id", "species", "phylum", "inserted_len", "n_orfs",
            "dominant_orf", "dominant_orf_len_aa", "dominant_orf_frac",
            "cds_cluster_id", "coding_frac", "noncoding_regions_bp"]
    out = os.path.join(args.out, "insert_cds.tsv")
    n_tnp = 0
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for rec, L in ins_len.items():
            lst = orfs.get(rec, [])
            d = dom.get(rec)
            cov_bp = 0
            for _, st, en, _ in lst:
                cov_bp += max(0, min(en, L) - max(st, 1) + 1)
            cov_bp = min(cov_bp, L)
            t = cid.get(rep.get(d[0], ""), ".") if d else "."
            if t != ".":
                n_tnp += 1
            row = {"db_id": rec,
                   "species": db.get(rec, {}).get("species", "."),
                   "phylum": db.get(rec, {}).get("phylum", "."),
                   "inserted_len": L,
                   "n_orfs": len(lst),
                   "dominant_orf": d[0] if d else ".",
                   "dominant_orf_len_aa": (d[3] // 3) if d else 0,
                   "dominant_orf_frac": "%.3f" % (d[3] / L) if d and L else "0",
                   "cds_cluster_id": t,
                   "coding_frac": "%.3f" % (cov_bp / L) if L else "0",
                   "noncoding_regions_bp": L - cov_bp}
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    sz = collections.Counter(cid[rep[k]] for k in rep if rep[k] in cid)
    bysp = collections.defaultdict(set)
    for rec, d in dom.items():
        t = cid.get(rep.get(d[0], ""))
        if t:
            bysp[t].add(db.get(rec, {}).get("species", "?"))

    log("=" * 76)
    log("CDS CLUSTERING")
    log("")
    log("  inserts                       %7d" % n_ins)
    log("  with >=1 ORF                  %7d  %5.1f%%"
        % (len(dom), 100.0 * len(dom) / n_ins if n_ins else 0))
    log("  assigned a cds_cluster_id     %7d  %5.1f%%"
        % (n_tnp, 100.0 * n_tnp / n_ins if n_ins else 0))
    log("  CDS clusters                  %7d" % len(cid))
    log("")
    log("  cluster sizes: singletons %d   >=10 %d   >=100 %d   max %d"
        % (sum(1 for v in sz.values() if v == 1),
           sum(1 for v in sz.values() if v >= 10),
           sum(1 for v in sz.values() if v >= 100),
           max(sz.values()) if sz else 0))
    log("")
    log("  CROSS-SPECIES clusters (structural evidence of transfer):")
    for n in (2, 3, 5):
        log("    in >=%d species   %6d" % (n, sum(1 for v in bysp.values() if len(v) >= n)))
    log("=" * 76)
    log("wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

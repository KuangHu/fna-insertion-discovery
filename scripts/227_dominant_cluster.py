#!/usr/bin/env python3
"""Pick each insertion's dominant CDS and cluster those proteins into bags.

THIS KEEPS THE RATIFIED CONVENTION AND DOES NOT RE-OPEN IT.  Identity of an
insertion is the cluster of its LONGEST Prodigal CDS, clustered as PROTEIN at
mmseqs --min-seq-id 0.50 -c 0.80 --cov-mode 0.  That is the rule 111_ set and
the user confirmed; nothing here changes it.  What changes is only the input:
the CDS now come from the validated occurrence database, deduplicated to
distinct insertion DNA, with translation table 11 fixed and Prodigal running
in meta mode with open ends.

A CLUSTER IS A CLUSTER.  No HMM, no IS library, no functional claim is used to
form it.  Pfam, where it has been run, is attached afterwards as a description
of a bag, never as a condition for belonging to one.

TWO BAG SETS ARE EMITTED, AND THE DIFFERENCE BETWEEN THEM IS THE POINT.
  strict     dominant CDS = longest Prodigal CDS.  Insertions where Prodigal
             called nothing have NO bag.
  augmented  the same, except that an insertion with no Prodigal CDS may take
             a context-rescued CDS (225, status complete_candidate, not
             likely_host_gene) as its dominant CDS.
Both are written.  Reporting only the augmented set would hide how much of it
rests on rescue rather than on prediction; reporting only the strict set would
throw away the insertions the extraction boundary broke.

TIES.  Longest wins; equal length is broken by Prodigal score, then by call_id,
so the choice is reproducible rather than dependent on iteration order.

PER-INSERTION FIELDS, as 111_ defined them:
  n_cds               how often the dominant-ORF question even arises
  dominant_cds_frac   how much of the insertion the chosen CDS covers -- the
                      honest weakness indicator for this heuristic
  coding_frac         union of all CDS spans over the insertion length
  noncoding_regions   the insertion minus its CDS spans
"""
import argparse, collections, json, os, subprocess, sys


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
                rid = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if rid:
        yield rid, "".join(parts)


def merge(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--faa", required=True, help="prodigal proteins")
    ap.add_argument("--fasta", required=True, help="dedup_insertions.fna")
    ap.add_argument("--rescue", default="", help="optional *_rescue.tsv")
    ap.add_argument("--rescue-faa", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-seq-id", type=float, default=0.50)
    ap.add_argument("--cov", type=float, default=0.80)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    L = {}
    for rid, s in fasta(a.fasta):
        L[rid.split("|", 1)[1] if "|" in rid else rid] = len(s)
    prot = {rid: s for rid, s in fasta(a.faa)}

    cds = collections.defaultdict(list)
    for r in read_tsv(a.calls):
        if r["route"] != "prodigal":
            continue
        cds[r["seq_key"]].append(r)

    dom, info = {}, {}
    for key, rows in cds.items():
        rows.sort(key=lambda r: (-int(r["length_nt"]), -float(r["score"] or 0), r["call_id"]))
        d = rows[0]
        dom[key] = (d["call_id"], prot.get(d["call_id"], ""), "prodigal")
        spans = merge([[int(r["nt_start0"]), int(r["nt_end0"])] for r in rows])
        cov = sum(e - s for s, e in spans)
        info[key] = dict(n_cds=len(rows), dominant_cds=d["call_id"],
                         dominant_cds_len_aa=len(prot.get(d["call_id"], "")),
                         dominant_cds_frac=round(int(d["length_nt"]) / L[key], 3),
                         coding_frac=round(cov / L[key], 3),
                         noncoding_regions=";".join(
                             f"{s}-{e}" for s, e in
                             [(0, spans[0][0])] * (spans[0][0] > 0) +
                             [(spans[i][1], spans[i + 1][0]) for i in range(len(spans) - 1)] +
                             [(spans[-1][1], L[key])] * (spans[-1][1] < L[key])),
                         dominant_cds_source="prodigal",
                         partial_begin=d["partial_begin"], partial_end=d["partial_end"])

    strict_keys = set(dom)
    aug = dict(dom)
    if a.rescue and os.path.exists(a.rescue):
        rprot = {rid: s for rid, s in fasta(a.rescue_faa)} if a.rescue_faa else {}
        best = {}
        for r in read_tsv(a.rescue):
            if r["seq_key"] in strict_keys:
                continue
            if r["status"] != "complete_candidate" or r["likely_host_gene"] == "1":
                continue
            pid = f"R|{r['seq_key']}|{r['pfam_acc']}|{r['domain_contig_start0']}-{r['domain_contig_end0']}"
            p = rprot.get(pid, "")
            if not p:
                continue
            cur = best.get(r["seq_key"])
            if cur is None or len(p) > len(cur[1]):
                best[r["seq_key"]] = (pid, p, "context_rescued")
        aug.update(best)
        for k, v in best.items():
            info[k] = dict(n_cds=0, dominant_cds=v[0], dominant_cds_len_aa=len(v[1]),
                           dominant_cds_frac=round(len(v[1]) * 3 / L[k], 3),
                           coding_frac=round(len(v[1]) * 3 / L[k], 3),
                           noncoding_regions="", dominant_cds_source="context_rescued",
                           partial_begin="0", partial_end="0")

    res = {}
    for tag, table in (("strict", dom), ("augmented", aug)):
        dfa = os.path.join(a.out, f"dominant_{tag}.faa")
        with open(dfa, "w") as f:
            for key in sorted(table):
                cid, p, src = table[key]
                if p:
                    f.write(f">{key} {cid} src={src}\n{p}\n")
        pre = os.path.join(a.out, f"clu_{tag}")
        tmp = os.path.join(a.out, f"mmtmp_{tag}")
        cmd = ["mmseqs", "easy-cluster", dfa, pre, tmp,
               "--min-seq-id", str(a.min_seq_id), "-c", str(a.cov),
               "--cov-mode", "0", "--threads", str(a.threads), "-v", "1"]
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        clu = {}
        with open(pre + "_cluster.tsv") as f:
            for ln in f:
                rep, mem = ln.rstrip("\n").split("\t")
                clu[mem] = rep
        joined = sum(1 for k in table if k in clu)
        if joined < 0.5 * len(table):
            sys.exit(f"FATAL: only {joined}/{len(table)} dominant CDS joined a cluster")
        bagid = {}
        for i, rep in enumerate(sorted(set(clu.values())), 1):
            bagid[rep] = f"BAG{tag[0].upper()}{i:06d}"
        with open(os.path.join(a.out, f"insertion_cds_{tag}.tsv"), "w") as f:
            cols = ["seq_key", "insert_len_nt", "cds_cluster_id", "cluster_rep",
                    "dominant_cds", "dominant_cds_source", "dominant_cds_len_aa",
                    "dominant_cds_frac", "n_cds", "coding_frac",
                    "partial_begin", "partial_end", "noncoding_regions"]
            f.write("\t".join(cols) + "\n")
            for key in sorted(table):
                it = info[key]
                rep = clu.get(key, "")
                f.write("\t".join(str(x) for x in [
                    key, L[key], bagid.get(rep, ""), rep, it["dominant_cds"],
                    it["dominant_cds_source"], it["dominant_cds_len_aa"],
                    it["dominant_cds_frac"], it["n_cds"], it["coding_frac"],
                    it["partial_begin"], it["partial_end"], it["noncoding_regions"]]) + "\n")
        res[tag] = dict(insertions_with_dominant_cds=len(table),
                        clusters=len(bagid), joined=joined)
        print(f"[{tag}] {len(table)} dominant CDS -> {len(bagid)} clusters", flush=True)

    res["insertions_total"] = len(L)
    res["insertions_without_any_prodigal_cds"] = len(L) - len(dom)
    res["params"] = dict(min_seq_id=a.min_seq_id, cov=a.cov, cov_mode=0)
    with open(os.path.join(a.out, "cluster.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()

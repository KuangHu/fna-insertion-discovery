#!/usr/bin/env python3
"""Scan both routes' proteins against Pfam-A under one fixed set of rules.

ONE PFAM, ONE HMMER, ONE THRESHOLD, BOTH ROUTES.  If route A and route B were
scanned with different cutoffs, every difference between them would be partly
an artefact of the scan.  So this script takes an arbitrary number of .faa
files and searches them in a SINGLE pass over the HMM library, with identical
options, writing one table whose `route` column is the only thing that differs.

THRESHOLD: each family's own GATHERING THRESHOLD (`--cut_ga`).  Pfam curates
GA per family precisely because one global E-value is wrong for 24,000 models
of very different length and information content.

COVERAGE IS RECORDED, NEVER USED AS A FILTER.  Demanding, say, >=80% HMM
coverage would discard exactly the short domain fragments this comparison
exists to detect -- a CDS cut by an extraction boundary produces a partial
domain, and that partial domain is the evidence.  So `hmm_cov` is a column.
Full vs fragment is a later interpretation, made explicitly in 224, not a
silent gate here.

SENSITIVE RE-PASS.  `--max` turns off HMMER's heuristic filters (F1/F2/F3 and
the bias filter) while KEEPING the gathering threshold, which is the official
way to trade speed for sensitivity without lowering the bar for what counts as
a hit.  It is reported separately (`scan_mode` column) and never merged into
the main result, because a hit found only with filters off is a weaker
statement than one found with them on.

Hits are per DOMAIN, not per protein: one protein may carry several families,
and a family may occur twice.  Each row is one included domain.
"""
import argparse, glob, json, os, time, zlib

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile

def dec(x):
    """pyhmmer returns bytes or str for names depending on build; accept both."""
    return x.decode() if isinstance(x, (bytes, bytearray)) else x


def enc(x):
    return x if isinstance(x, (bytes, bytearray)) else str(x).encode()


COLS = ["protein_id", "route", "seq_key", "pfam_acc", "pfam_name", "scan_mode",
        "bitscore", "i_evalue", "c_evalue", "hmm_len", "hmm_from", "hmm_to",
        "hmm_cov", "ali_from", "ali_to", "env_from", "env_to", "prot_len",
        "domain_index", "n_domains_in_protein"]


# The route is carried by the call_id prefix that 221/225 assign, so it stays
# correct no matter what a file is named; the filename is only a fallback.
PREFIX_ROUTE = {"A": "prodigal", "B": "sixframe", "R": "context_rescued"}


def load_faa(paths, alphabet):
    seqs, meta = [], {}
    for p in paths:
        b = os.path.basename(p)
        fallback = ("prodigal" if "_prodigal" in b else
                    "context_rescued" if "_rescued" in b else "sixframe")
        name, parts = None, []
        with open(p) as f:
            for line in f:
                if line[0] == ">":
                    if name:
                        s = "".join(parts).replace("*", "")
                        seqs.append(TextSequence(name=name.encode(), sequence=s)
                                    .digitize(alphabet))
                        meta[name] = (PREFIX_ROUTE.get(name.split("|", 1)[0], fallback),
                                      len(s))
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
        if name:
            s = "".join(parts).replace("*", "")
            seqs.append(TextSequence(name=name.encode(), sequence=s).digitize(alphabet))
            meta[name] = (PREFIX_ROUTE.get(name.split("|", 1)[0], fallback), len(s))
    return seqs, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", nargs="+", required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=0)
    ap.add_argument("--max", action="store_true",
                    help="disable heuristic filters (hmmsearch --max); GA is kept")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()

    alpha = Alphabet.amino()
    paths = []
    for p in a.faa:
        paths.extend(sorted(glob.glob(p)))
    seqs, meta = load_faa(paths, alpha)
    ntotal = len(seqs)
    if a.nshard > 1:
        # Python's str hash is salted per process, so it cannot partition work;
        # crc32 is stable across processes and runs.
        seqs = [s for s in seqs
                if zlib.crc32(enc(s.name)) % a.nshard == a.shard % a.nshard]
    print(f"proteins {len(seqs)}/{ntotal} (shard {a.shard}/{a.nshard}) "
          f"from {len(paths)} files", flush=True)

    # Z is pinned to the FULL protein count, not this shard's, so `i_evalue` is
    # the value an unsharded run would report. Inclusion is by GA bit score and
    # is unaffected either way, but a shard-dependent E-value column would be a
    # trap for anyone who later filtered on it.
    # domZ is deliberately NOT pinned: the conditional E-value is defined
    # relative to the sequences that passed THIS search, and forcing domZ=Z
    # would collapse c_evalue onto i_evalue and destroy the distinction.
    # c_evalue is therefore search-scoped; compare bit scores, not c_evalue,
    # across shards.
    opts = dict(bit_cutoffs="gathering", Z=float(ntotal))
    mode = "ga"
    if a.max:
        opts.update(F1=1.0, F2=1.0, F3=1.0, bias_filter=False)
        mode = "ga_max"

    nhmm, nhit = 0, 0
    per_prot = {}
    rows = []
    with HMMFile(a.hmm) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=a.cpus, **opts):
            nhmm += 1
            acc = dec(hits.query.accession or hits.query.name)
            nm = dec(hits.query.name)
            M = hits.query.M
            for hit in hits:
                if not hit.included:
                    continue
                pid = dec(hit.name)
                route, plen = meta[pid]
                doms = [d for d in hit.domains if d.included]
                for di, d in enumerate(doms, 1):
                    al = d.alignment
                    cov = (al.hmm_to - al.hmm_from + 1) / M
                    rows.append([pid, route, pid.split("|")[1] if "|" in pid else "",
                                 acc, nm, mode, f"{d.score:.2f}", f"{d.i_evalue:.3g}",
                                 f"{d.c_evalue:.3g}", M, al.hmm_from, al.hmm_to,
                                 f"{cov:.4f}", al.target_from, al.target_to,
                                 d.env_from, d.env_to, plen, di, len(doms)])
                    nhit += 1
                per_prot[pid] = per_prot.get(pid, 0) + len(doms)
            if nhmm % 2000 == 0:
                print(f"  {nhmm} models, {nhit} domain hits, "
                      f"{time.time()-t0:.0f}s", flush=True)

    with open(a.out, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in sorted(rows, key=lambda r: (r[0], int(r[13]))):
            f.write("\t".join(str(x) for x in r) + "\n")

    vf = os.path.join(os.path.dirname(a.hmm), "Pfam.version")
    release = open(vf).read().strip() if os.path.exists(vf) else "UNKNOWN"
    summ = dict(hmm_models=nhmm, proteins=len(seqs), proteins_total=ntotal,
                shard=a.shard, nshard=a.nshard, proteins_with_hit=len(per_prot),
                domain_hits=nhit, scan_mode=mode, hmm_db=a.hmm,
                pfam_release=release,
                pyhmmer=pyhmmer.__version__, seconds=round(time.time() - t0, 1))
    with open(a.out + ".json", "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()

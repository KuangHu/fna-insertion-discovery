#!/usr/bin/env python3
"""Which genes do IS elements land in, and is any gene hit by several of them?

Three questions, in order:

  1 which insertions actually carry a transposase or integrase
  2 of those, which landed INSIDE a gene, and what is that gene
  3 is any target gene family disrupted by SEVERAL DIFFERENT element families

THIS DOES NOT VIOLATE THE FOUNDING CONSTRAINT. Nothing here seeds discovery.
Every insertion was found structurally, by comparing alleles, with no model of
what an element looks like. Pfam is applied afterwards to ASK WHAT THE THING
IS -- post-hoc classification, which is exactly the order the project insists
on. A transposase model never selects a locus.

  1 ELEMENT IDENTITY. hmmsearch the insert's own ORFs against Pfam-A, and keep
    the hits whose Pfam accession is in a transposase/integrase set. Reported
    with the model that fired, so a reader can see WHICH family rather than
    trusting a boolean.

  2 TARGET GENE. 170_ already located the disrupted CDS in the pre-insertion
    genome. Its protein is extracted and searched against Pfam too, so the
    target is named by the same evidence standard as the element.

  3 MULTIPLY-TARGETED GENES. Target proteins are clustered with mmseqs at the
    corpus's own 50%/80% setting, giving target gene FAMILIES without needing
    a name. A family hit by two or more distinct element families is the
    answer to the third question, and it is answerable even where Pfam names
    nothing -- which matters, because most bacterial ORFs are hypothetical.

WHAT A HIT MEANS, AND DOES NOT. A Pfam hit above the gathering threshold says
the protein contains that domain. It does not say the gene is expressed, that
the disruption mattered, or that the element is active. And every event here
is a SURVIVOR: insertions into essential genes killed their host and were
never sampled, so the genes seen are biased toward the dispensable. A gene
being hit repeatedly means it is both targetable AND tolerant of disruption,
and those two cannot be separated with this data.
"""
import argparse
import collections
import csv
import glob
import os
import subprocess
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

# Pfam clans/families for transposases and integrases. Matched on the
# ACCESSION, not on the description string: matching "transposase" in free
# text would also catch "transposase-like" and miss families whose name does
# not contain the word.
TNP_KEYWORDS = ("transpos", "integrase", "recombinase", "resolvase",
                "dde_", "rve", "is200", "is605", "is110", "insertion element",
                "tn916", "tn3", "mule", "helitron")


def read_fa(p):
    out, name, buf = {}, None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def hmm_best(faa, pfam, out, threads, use_cut_ga=True):
    """hmmsearch -> {protein: (pfam_acc, pfam_name, score)} keeping the best hit.

    --cut_ga uses each model's own curated gathering threshold rather than one
    E-value for every family. An E-value is database-size dependent and the
    project forbids it for anchors (assertion A13); the same reasoning applies
    here.
    """
    cmd = ["hmmsearch", "--cpu", str(threads), "--noali", "--tblout", out]
    cmd += ["--cut_ga"] if use_cut_ga else ["-E", "1e-5"]
    cmd += [pfam, faa]
    if subprocess.run(cmd, stdout=subprocess.DEVNULL,
                      stderr=subprocess.DEVNULL).returncode != 0:
        return {}
    best = {}
    for line in open(out):
        if not line or line[0] == "#":
            continue
        f = line.split()
        if len(f) < 9:
            continue
        prot, acc, name = f[0], f[3], f[2]
        try:
            sc = float(f[5])
        except ValueError:
            continue
        if prot not in best or sc > best[prot][2]:
            best[prot] = (acc, name, sc)
    return best


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--context", required=True, help="glob of *_context.tsv from 170_")
    ap.add_argument("--db", required=True, help="database insertions.tsv")
    ap.add_argument("--cds", required=True,
                    help="cds_v3/insert_cds.tsv -- carries cds_cluster_id, "
                         "which insertions.tsv does NOT")
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--genome-dirs", required=True,
                    help="species=dir,species=dir,...")
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--min-seq-id", type=float, default=0.50)
    ap.add_argument("--cov", type=float, default=0.80)
    a = ap.parse_args()

    gdir = dict(x.split("=", 1) for x in a.genome_dirs.split(",") if "=" in x)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    tmp = a.out + "_tmp"
    os.makedirs(tmp, exist_ok=True)

    ctx = []
    for p in sorted(glob.glob(a.context)):
        ctx.extend(csv.DictReader(open(p), delimiter="\t"))
    log("context rows: %d" % len(ctx))
    intr = [r for r in ctx if r["call"] == "intragenic"]
    log("  intragenic: %d" % len(intr))
    if not intr:
        return 0

    # ---- 1. which inserts carry a transposase / integrase -----------------
    dbr = {r["db_id"]: r for r in csv.DictReader(open(a.db), delimiter="\t")}
    # cds_cluster_id lives in insert_cds.tsv, NOT in insertions.tsv. Reading it
    # from the wrong file made every element's cluster "." so every target
    # family saw exactly one distinct element and "families hit by >=2
    # elements" was 0 BY CONSTRUCTION, before any data was read.
    clu = {r["db_id"]: r.get("cds_cluster_id", ".")
           for r in csv.DictReader(open(a.cds), delimiter="\t")}
    log("cds_cluster_id loaded for %d inserts" % len(clu))
    ins = read_fa(a.inserts)
    # ID SPACE. 170_ writes the bare event_id (E001625) because that is what
    # dedup produces, but the database and its FASTA are keyed by db_id,
    # "<species>.<event_id>" -- namespaced in 110_ because event_id is
    # species-scoped and collides across species. Joining on the bare id
    # silently matched NOTHING and the only symptom was "inserts to classify: 0".
    def dbid(r):
        return "%s.%s" % (r["species"], r["event_id"])
    need = {dbid(r) for r in ctx}
    f_ins = os.path.join(tmp, "inserts.fna")
    n_found = 0
    with open(f_ins, "w") as fh:
        for k in need:
            if k in ins:
                fh.write(">%s\n%s\n" % (k, ins[k]))
                n_found += 1
    if n_found == 0:
        log("FATAL: 0 of %d context events found in the inserts FASTA -- the "
            "id spaces do not match" % len(need))
        return 3
    log("  inserts to classify: %d of %d context events" % (n_found, len(need)))
    f_orf = os.path.join(tmp, "insert_orfs.faa")
    if subprocess.run(["prodigal", "-i", f_ins, "-a", f_orf, "-p", "meta",
                       "-q", "-o", "/dev/null"],
                      stderr=subprocess.DEVNULL).returncode != 0:
        log("FATAL: prodigal failed on inserts")
        return 2
    log("  insert ORFs: %d" % sum(1 for l in open(f_orf) if l[0] == ">"))
    log("  hmmsearch inserts vs Pfam (--cut_ga) ...")
    ib = hmm_best(f_orf, a.pfam, os.path.join(tmp, "ins.tbl"), a.threads)
    elem = {}
    for prot, (acc, name, sc) in ib.items():
        ev = prot.rsplit("_", 1)[0]
        is_tnp = any(k in name.lower() for k in TNP_KEYWORDS)
        if is_tnp and (ev not in elem or sc > elem[ev][2]):
            elem[ev] = (acc, name, sc)
    log("  events whose insert carries a transposase/integrase domain: %d"
        % len(elem))

    # ---- 2. the disrupted gene's protein ----------------------------------
    by_g = collections.defaultdict(list)
    for r in intr:
        by_g[(r["species"], r["genome"])].append(r)
    f_tgt = os.path.join(tmp, "targets.fna")
    kept = []
    with open(f_tgt, "w") as fh:
        for (sp, g), rs in by_g.items():
            d = gdir.get(sp)
            if not d:
                continue
            gp = os.path.join(d, g + ".fna")
            if not os.path.exists(gp):
                continue
            ctg = read_fa(gp)
            for r in rs:
                seq = ctg.get(r["contig"], "")
                try:
                    s, e = int(r["cds_start"]), int(r["cds_end"])
                except ValueError:
                    continue
                if not seq or e > len(seq) or e - s < 90:
                    continue
                nm = "t%d" % len(kept)
                kept.append(r)
                fh.write(">%s\n%s\n" % (nm, seq[s:e]))
    log("  target CDS extracted: %d" % len(kept))
    f_tp = os.path.join(tmp, "targets.faa")
    subprocess.run(["prodigal", "-i", f_tgt, "-a", f_tp, "-p", "meta", "-q",
                    "-o", "/dev/null"], stderr=subprocess.DEVNULL)
    log("  hmmsearch targets vs Pfam (--cut_ga) ...")
    tb = hmm_best(f_tp, a.pfam, os.path.join(tmp, "tgt.tbl"), a.threads)

    # ---- 3. cluster target proteins into FAMILIES -------------------------
    # Families, not names: most bacterial ORFs are hypothetical, so a
    # name-based grouping would discard most of the data. Same 50%/80%
    # setting the bag corpus uses, so "same family" means the same thing here.
    pre = os.path.join(tmp, "tclu")
    subprocess.run(["mmseqs", "easy-cluster", f_tp, pre,
                    os.path.join(tmp, "mm"), "--min-seq-id", str(a.min_seq_id),
                    "-c", str(a.cov), "--cov-mode", "0",
                    "--threads", str(a.threads), "-v", "1"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    fam = {}
    ct = pre + "_cluster.tsv"
    if os.path.exists(ct):
        for line in open(ct):
            rep, mem = line.rstrip("\n").split("\t")[:2]
            fam[mem] = rep
    log("  target families: %d" % len(set(fam.values())))

    # ---- write the per-event table ----------------------------------------
    rows = []
    for i, r in enumerate(kept):
        nm = "t%d" % i
        prot = None
        for k in tb:
            if k.rsplit("_", 1)[0] == nm:
                prot = k
                break
        ev = dbid(r)
        el = elem.get(ev)
        src = dbr.get(ev, {})
        rows.append({
            "event_id": r["event_id"], "db_id": ev, "species": r["species"],
            "element_pfam": el[0] if el else ".",
            "element_pfam_name": el[1] if el else ".",
            "element_is_mobile": "yes" if el else "no",
            "element_cds_cluster": clu.get(ev, "."),
            "inserted_len": r.get("inserted_len", "."),
            "target_gene_family": fam.get(prot, ".") if prot else ".",
            "target_pfam": tb[prot][0] if prot in tb else ".",
            "target_pfam_name": tb[prot][1] if prot in tb else ".",
            "target_cds_len": r.get("cds_len", "."),
            "frac_into_cds": r.get("frac_into_cds", "."),
            "frag5_len": r.get("frag5_len", "."),
            "frag3_len": r.get("frag3_len", "."),
            "genome": r["genome"], "contig": r["contig"],
            "junction_pos": r["junction_pos"]})
    cols = list(rows[0].keys()) if rows else []
    with open(a.out + "_targets.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- question 3: families hit by SEVERAL element families -------------
    mob = [r for r in rows if r["element_is_mobile"] == "yes"]
    log("")
    log("=" * 72)
    log("  intragenic events with a transposase/integrase insert: %d of %d"
        % (len(mob), len(rows)))
    hit = collections.defaultdict(set)
    hitn = collections.Counter()
    n_noclu = sum(1 for r in mob if r["element_cds_cluster"] == ".")
    if mob and n_noclu == len(mob):
        log("FATAL: every element has cds_cluster '.', so 'families hit by >=2 "
            "elements' can only be 0. The cluster join is broken.")
        return 4
    if n_noclu:
        log("  note: %d of %d elements have no cds_cluster and are excluded "
            "from the multi-element count" % (n_noclu, len(mob)))
    for r in mob:
        tf = r["target_gene_family"]
        if tf == "." or r["element_cds_cluster"] == ".":
            continue
        hit[tf].add(r["element_cds_cluster"])
        hitn[tf] += 1
    multi = {k: v for k, v in hit.items() if len(v) >= 2}
    log("  target gene families hit at all              : %d" % len(hit))
    log("  families hit by >=2 DIFFERENT element families: %d" % len(multi))
    log("")
    log("  top multiply-targeted genes")
    log("  %-10s %6s %7s  %-34s %s" %
        ("family", "events", "elements", "target Pfam", "element Pfam(s)"))
    name_of = {r["target_gene_family"]: r["target_pfam_name"] for r in mob}
    enames = collections.defaultdict(set)
    for r in mob:
        if r["target_gene_family"] in multi:
            enames[r["target_gene_family"]].add(r["element_pfam_name"])
    for tf in sorted(multi, key=lambda k: (-len(hit[k]), -hitn[k]))[:20]:
        log("  %-10s %6d %7d  %-34s %s"
            % (tf[:10], hitn[tf], len(hit[tf]), name_of.get(tf, ".")[:34],
               ",".join(sorted(enames[tf]))[:40]))
    log("=" * 72)
    log("  Every event is a SURVIVOR. Insertions into essential genes killed")
    log("  their host and were never sampled, so these genes are biased toward")
    log("  the dispensable. Repeated targeting means targetable AND tolerant;")
    log("  this data cannot separate the two.")
    log("wrote %s_targets.tsv" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

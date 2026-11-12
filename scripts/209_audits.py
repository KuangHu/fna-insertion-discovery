#!/usr/bin/env python3
"""The two audits promised earlier: junction anomalies, and delta != 0.

AUDIT A -- JUNCTION ANOMALIES. Re-derives the cases where the insert's own
junction is found in the assembly the catalogue calls EMPTY, and records WHERE
each hit lies rather than only that it exists. The question is whether the
EMPTY LOCUS is empty, not whether the whole assembly lacks that junction: a
mobile element is expected elsewhere in the same genome, so a hit far from the
recorded junction is not an anomaly. Distance to the recorded junction is
therefore the discriminating column.

AUDIT B -- delta != 0. The catalogue-level identity
    (len_derived_allele - len_empty_allele) - inserted_len
equals decompose()'s own ambiguity lcp + lcs - len(short), which is 0 for a
clean decomposition. The events where it is NOT zero are the ones where exact
prefix/suffix matching overran: the shared prefix and suffix together exceed
the short allele, so the two are not independent and the junction cannot be
placed from length alone. Each such event is classified by comparing the allele
pair directly, rather than all being filed as one bucket.
"""
import argparse
import collections
import os
import sqlite3
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    out, k, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if k is not None:
                out[k] = "".join(buf).upper()
            k = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf).upper()
    return out


def load_map(p):
    m = {}
    with open(p) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[f[0]] = f[1]
    return m


def audit_a(a):
    man = load_map(SRC + "/genome_map.tsv")
    ins = read_fasta(SRC + "/insert.fna")
    con = sqlite3.connect("file:%s/locations/insertions.sqlite?mode=ro" % W, uri=True)
    ev = list(con.execute(
        "SELECT db_id, species, coord_frame, coord_contig, junction_L_abs,"
        " junction_R_abs FROM loc WHERE junction_L_abs NOT IN ('.','')"
        " GROUP BY db_id"))
    by = collections.defaultdict(list)
    for r in ev:
        by[r[2]].append(r)
    gs = sorted(by)
    mine = [g for i, g in enumerate(gs) if i % a.nshard == a.shard]
    out = open("%s/junction_anomaly_%03d.tsv" % (a.out, a.shard), "w")
    K = a.probe
    n = 0
    for gi, g in enumerate(mine, 1):
        p = man.get(g)
        if not p:
            continue
        try:
            ctgs = read_fasta(p)
        except (IOError, OSError):
            continue
        whole = {k: v for k, v in ctgs.items()}
        for db_id, sp, cf, ct, jl, jr in by[g]:
            seq = ctgs.get(ct)
            iseq = ins.get(db_id, "")
            if seq is None or not iseq:
                continue
            jl, jr = int(jl), int(jr)
            probe = seq[max(0, jl - K):jl] + iseq[:K]
            if len(probe) < 2 * K:
                continue
            hits = []
            rp = revcomp(probe)
            for c, v in whole.items():
                for pat, st in ((probe, "+"), (rp, "-")):
                    i = v.find(pat)
                    while i >= 0 and len(hits) < 20:
                        hits.append((c, i, st))
                        i = v.find(pat, i + 1)
            if not hits:
                continue
            n += 1
            for c, i, st in hits:
                d = abs(i - jl) if c == ct else -1
                cls = ("at_recorded_junction" if c == ct and d <= 50 else
                       ("same_contig_distant" if c == ct else "other_contig"))
                out.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n"
                          % (db_id, sp, cf, ct, jl, c, i, st, d) )
            if gi % 200 == 0:
                sys.stderr.write("  %d/%d genomes\n" % (gi, len(mine)))
    out.close()
    print("audit A shard %d: %d events with a junction hit in the empty assembly" % (a.shard, n))


def audit_b(a):
    """Classify every event whose length identity is not zero."""
    con = sqlite3.connect("file:%s/locations/insertions.sqlite?mode=ro" % W, uri=True)
    emp, der = {}, {}
    for db_id, role, al in con.execute(
            "SELECT db_id,role,allele_len FROM loc WHERE allele_len NOT IN ('.','')"):
        try:
            al = int(al)
        except ValueError:
            continue
        (emp if role == "carrier_empty" else der).setdefault(db_id, set()).add(al)
    cat = {}
    for f in sorted(os.listdir(W + "/catalogue_v3")):
        if not f.endswith("_events.tsv"):
            continue
        sp = f[:-len("_events.tsv")]
        with open(W + "/catalogue_v3/" + f) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in fh:
                g = line.rstrip("\n").split("\t")
                cat["%s.%s" % (sp, g[ix["event_id"]])] = (
                    int(g[ix["target_len"]]), int(g[ix["inserted_len"]]),
                    int(g[ix["junction_overlap_bp"]]),
                    g[ix["decomposition_method"]],
                    int(g[ix["target_bases_lost"]] or 0),
                    g[ix["placement_status"]])
    cls = collections.Counter()
    rows = []
    for db_id in der:
        if db_id not in emp or db_id not in cat:
            continue
        if len(der[db_id]) != 1 or len(emp[db_id]) != 1:
            cls["multi_length_alleles_not_comparable"] += 1
            continue
        L, S = list(der[db_id])[0], list(emp[db_id])[0]
        tl, il, ov, meth, lost, pst = cat[db_id]
        d = (L - S) - il
        if d == 0:
            cls["zero_identity_holds"] += 1
            continue
        # d != 0 => decompose's own ambiguity is nonzero: lcp+lcs > len(short)
        # ORDER MATTERS. An earlier version tested the decomposition path
        # first, which swallowed every other category and produced a single
        # bucket. The specific, falsifiable conditions are tested first and the
        # path is only the fallback.
        suffix = "_tolerant" if meth == "tolerant_alignment" else "_exact"
        if lost > 0:
            c = "target_bases_lost_at_junction" + suffix
        elif pst == "ambiguous":
            c = "placement_ambiguous" + suffix
        elif d == ov and ov > 0:
            # DESCRIPTIVE ONLY -- THIS IS NOT AN EXPLANATION. Sliding moves the
            # interval without changing its length, so the slide span cannot
            # cause a LENGTH difference. delta == slide span is a coincidence
            # of two numbers, and these rows stay in the unexplained tally.
            c = "delta_coincides_with_slide_span_DESCRIPTIVE" + suffix
        elif d > 0:
            c = "positive_delta_unexplained" + suffix
        else:
            c = "negative_delta_unexplained" + suffix
        cls[c] += 1
        rows.append((db_id, L, S, il, ov, d, meth, lost, pst, c))
    with open(a.out + "/delta_nonzero_classified.tsv", "w") as fh:
        fh.write("db_id\tlen_derived\tlen_empty\tinserted_len\tslide_span\t"
                 "delta\tmethod\ttarget_bases_lost\tplacement\tclass\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    tot = sum(cls.values())
    print("AUDIT B -- length identity delta, n=%d" % tot)
    for k, v in cls.most_common():
        print("   %-46s %7d  %5.1f%%" % (k, v, 100.0 * v / tot))
    nz = sum(v for k, v in cls.items() if k not in
             ("zero_identity_holds", "multi_length_alleles_not_comparable"))
    # Only target_bases_lost and placement_ambiguous are candidate CAUSES of a
    # length difference. A coincidence with the slide span is not one.
    explained = sum(v for k, v in cls.items()
                    if k.startswith(("target_bases_lost", "placement_ambiguous")))
    unexp = nz - explained
    print("   delta != 0 total %d" % nz)
    print("   with a candidate cause (bases lost / placement ambiguous) %d" % explained)
    print("   STILL UNEXPLAINED %d  (includes delta==slide-span coincidences,"
          " which do not explain a length change)" % unexp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("audit", choices=["A", "B"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--probe", type=int, default=40)
    ap.add_argument("--out", default=W + "/audits")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    (audit_a if a.audit == "A" else audit_b)(a)


if __name__ == "__main__":
    main()

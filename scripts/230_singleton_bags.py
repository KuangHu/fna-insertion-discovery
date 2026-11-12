#!/usr/bin/env python3
"""What is in the one-site bags? Specifically: do any of them look like transposases?

A bag seen at exactly ONE independent site is the ambiguous case.  It could be
a genuinely rare or novel element, or a fragment, or a passenger gene that
happens to be the longest ORF in its insert, or an assembly artefact.  The
question asked here is narrower and answerable: how often does the DOMINANT
CDS -- the protein that actually defines the cluster -- carry a mobility
domain, and how does that rate compare with the bags seen many times?

THE MGE LABEL IS A KEYWORD MATCH AND IS TREATED AS ONE.  There is no curated
"transposase" list in Pfam, so families are matched on their ID and
description text from Pfam-A.hmm.dat.  That is fuzzy in both directions: it
can catch a family merely described as transposase-ASSOCIATED, and it will
miss a mobility protein whose name says nothing about mobility.  So the script
prints the full matched family list, with counts, for inspection -- the label
is an index into evidence, never the conclusion.

DOMINANT-CDS DOMAINS ARE SEPARATED FROM ANY-DOMAIN.  A bag is defined by the
cluster of its dominant CDS, so a transposase domain sitting on some OTHER ORF
in the insert says the insert carries a transposase, not that the BAG is one.
Both numbers are reported and never added together.

RATES ARE COMPARED AGAINST THE MULTI-SITE BAGS.  "412 one-site bags look like
transposases" means nothing on its own; what matters is whether one-site bags
are enriched or depleted in mobility domains relative to bags seen often.
"""
import argparse, collections, gzip, json, os, re

MGE_PAT = re.compile(
    r"transpos|integrase|recombinase|resolvase|invertase|\bDDE\b|DDE_|Tnp|"
    r"\brve\b|IS\d{2,4}|InsA|InsB|insertion sequence|insertion element|"
    r"MULE|helitron|Mu DNA|phage integrase|tyrosine recombinase|"
    r"serine recombinase|relaxase|conjugat|mobilis|\bTn\d",
    re.I)


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def load_pfam_dat(p):
    """acc -> (id, description, clan). Reads the plain or gzipped .dat."""
    op = gzip.open if p.endswith(".gz") else open
    out, cur = {}, {}
    with op(p, "rt") as f:
        for line in f:
            if line.startswith("#=GF ID"):
                cur = {"id": line.split(None, 2)[2].strip()}
            elif line.startswith("#=GF AC"):
                cur["acc"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF DE"):
                cur["de"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF CL"):
                cur["cl"] = line.split(None, 2)[2].strip()
            elif line.startswith("//") and cur.get("acc"):
                out[cur["acc"].split(".")[0]] = (
                    cur.get("id", ""), cur.get("de", ""), cur.get("cl", ""))
                cur = {}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds", required=True, help="insertion_cds_strict.tsv")
    ap.add_argument("--bags", required=True, help="strict_bags.tsv")
    ap.add_argument("--hits", required=True, help="all_hits_on_dna.tsv")
    ap.add_argument("--calls", required=True, help="all_calls.tsv")
    ap.add_argument("--occ", required=True, help="strict_occurrence_bags.tsv")
    ap.add_argument("--dat", required=True, help="Pfam-A.hmm.dat[.gz]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    dat = load_pfam_dat(a.dat)
    mge = {acc for acc, (i, d, c) in dat.items() if MGE_PAT.search(i + " " + d)}
    print(f"Pfam families {len(dat)}; keyword-matched as mobility-related {len(mge)}",
          flush=True)

    cds = {r["seq_key"]: r for r in read_tsv(a.cds)}
    bag_of = {k: r["cds_cluster_id"] for k, r in cds.items() if r["cds_cluster_id"]}
    dom_call = {k: r["dominant_cds"] for k, r in cds.items()}
    sites = {r["bag_id"]: int(r["n_independent_sites"]) for r in read_tsv(a.bags)}

    partial = {}
    for r in read_tsv(a.calls):
        if r["route"] == "prodigal":
            partial[r["call_id"]] = (r["partial_begin"], r["partial_end"])

    dom_fams = collections.defaultdict(set)     # bag -> families on dominant CDS
    any_fams = collections.defaultdict(set)     # bag -> families anywhere in the insert
    key_dom = collections.defaultdict(set)
    for h in read_tsv(a.hits):
        b = bag_of.get(h["seq_key"])
        if not b:
            continue
        acc = h["pfam_acc"].split(".")[0]
        any_fams[b].add(acc)
        if h["protein_id"] == dom_call.get(h["seq_key"]):
            dom_fams[b].add(acc)
            key_dom[h["seq_key"]].add(acc)

    def band(n):
        if n == 1:
            return "1 site"
        if n <= 3:
            return "2-3 sites"
        if n <= 9:
            return "4-9 sites"
        if n <= 19:
            return "10-19 sites"
        return ">=20 sites"

    tab = collections.defaultdict(lambda: collections.Counter())
    for b, n in sites.items():
        g = band(n)
        tab[g]["bags"] += 1
        d, y = dom_fams.get(b, set()), any_fams.get(b, set())
        tab[g]["dominant_has_any_pfam"] += bool(d)
        tab[g]["dominant_has_mge"] += bool(d & mge)
        tab[g]["insert_has_any_pfam"] += bool(y)
        tab[g]["insert_has_mge"] += bool(y & mge)
        tab[g]["no_pfam_anywhere"] += (not y)

    order = ["1 site", "2-3 sites", "4-9 sites", "10-19 sites", ">=20 sites"]
    print(f"\n{'band':>12} {'bags':>7} {'dom Pfam':>10} {'dom MGE':>9} "
          f"{'%MGE':>6} {'insert MGE':>11} {'no Pfam':>8}")
    rows = {}
    for g in order:
        t = tab[g]
        pct = 100 * t["dominant_has_mge"] / max(1, t["bags"])
        rows[g] = dict(t)
        rows[g]["pct_dominant_mge"] = round(pct, 1)
        print(f"{g:>12} {t['bags']:>7} {t['dominant_has_any_pfam']:>10} "
              f"{t['dominant_has_mge']:>9} {pct:>5.1f}% {t['insert_has_mge']:>11} "
              f"{t['no_pfam_anywhere']:>8}")

    # the one-site bags whose dominant CDS carries a mobility domain
    occ1 = {}
    for r in read_tsv(a.occ):
        if sites.get(r["cds_cluster_id"]) == 1 and r["cds_cluster_id"] not in occ1:
            occ1[r["cds_cluster_id"]] = r
    singles = [b for b, n in sites.items() if n == 1]
    fam_count = collections.Counter()
    for b in singles:
        for f in dom_fams.get(b, ()):
            fam_count[f] += 1

    with open(a.out + "_singleton_mge.tsv", "w") as f:
        f.write("bag_id\tseq_key\toccurrence_id\tspecies\tinsert_len_nt\t"
                "dominant_cds_len_aa\tdominant_cds_frac\tpartial_begin\tpartial_end\t"
                "n_cds\tdominant_families\tdominant_family_names\n")
        for b in sorted(singles):
            d = dom_fams.get(b, set())
            if not (d & mge):
                continue
            o = occ1.get(b, {})
            k = o.get("seq_key", "")
            c = cds.get(k, {})
            pb, pe = partial.get(c.get("dominant_cds", ""), ("", ""))
            f.write("\t".join([
                b, k, o.get("occurrence_id", ""), o.get("species", ""),
                c.get("insert_len_nt", ""), c.get("dominant_cds_len_aa", ""),
                c.get("dominant_cds_frac", ""), pb, pe, c.get("n_cds", ""),
                ";".join(sorted(d)),
                ";".join(dat.get(x, ("?",))[0] for x in sorted(d))]) + "\n")

    with open(a.out + "_singleton_families.tsv", "w") as f:
        f.write("pfam_acc\tpfam_id\tdescription\tclan\tis_mge_keyword\tn_singleton_bags\n")
        for acc, n in fam_count.most_common():
            i, d, c = dat.get(acc, ("?", "?", ""))
            f.write(f"{acc}\t{i}\t{d}\t{c}\t{int(acc in mge)}\t{n}\n")

    with open(a.out + "_singleton.json", "w") as f:
        json.dump(dict(bands=rows, mge_families_in_pfam=len(mge),
                       singleton_bags=len(singles),
                       singleton_dominant_mge=sum(
                           1 for b in singles if dom_fams.get(b, set()) & mge)), f, indent=2)
    print(f"\none-site bags: {len(singles)}; dominant CDS carries a mobility family in "
          f"{sum(1 for b in singles if dom_fams.get(b, set()) & mge)}")


if __name__ == "__main__":
    main()

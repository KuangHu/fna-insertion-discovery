#!/usr/bin/env python3
"""Four objects, kept separate, at record level. No summary score.

Every previous attempt collapsed four different things into one "background"
count and then argued about the number. They are not the same object and do
not have the same cardinality at one locus:

  physical_site     assembly, contig, interval, orientation. One place in one
                    assembly. Deduplicated on the actual alignment interval,
                    not a position bin.
  element_instance  the sequence actually present at that site, its boundary
                    evidence, and EVERY seed that hit it. When several element
                    families hit the same interval the conflict is RECORDED,
                    not resolved -- family is marked ambiguous and the site
                    still deduplicates.
  host_locus        the same locus across strains, by flank homology. Carries
                    the flank distance that justified the merge so the
                    grouping can be re-cut at a different tolerance.
  pairing_variant   one (element instance, target context) combination. A host
                    locus with three element versions has one host_locus and
                    three pairing_variants, and that is not a contradiction.

WHAT IS DELIBERATELY NOT DONE. No score, no per-core-hour figure, no claim
about which method wins. Those were withdrawn: the counts feeding them mixed
these objects together. This emits records a person can read one at a time.

AMBIGUITY IS A VALUE, NOT A FAILURE. A site hit by six element families gets
family=AMBIGUOUS with all six listed. Previous versions either merged them
(losing the conflict) or counted them separately (inflating ~6x).
"""
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)
W = "/global/scratch/users/kh36969/fna_ins_discovery"
SP = sys.argv[1] if len(sys.argv) > 1 else "spneumoniae"
OUT = sys.argv[2] if len(sys.argv) > 2 else W + "/xspecies/%s_objects" % SP
OVL = 0.90      # interval overlap to call two hits the same physical site


def iv_same(a, b):
    """Two alignment intervals are one physical site if they reciprocally
    overlap >= OVL. Replaces position//50 binning, which split a hit whose
    boundary straddled a bin edge and merged two starts inside one bin."""
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    if e <= s:
        return False
    o = e - s
    return o / (a[1] - a[0]) >= OVL and o / (b[1] - b[0]) >= OVL


rows = list(csv.DictReader(open(W + "/xspecies/%s_hits.tsv" % SP), delimiter="\t"))
print("hit records: %d" % len(rows))

# ---- physical sites: cluster hits per (genome, contig) by interval overlap
bykey = collections.defaultdict(list)
for r in rows:
    bykey[(r["genome"], r["contig"])].append(r)

sites = []
for (g, c), rs in bykey.items():
    rs.sort(key=lambda x: int(x["elem_start"]))
    cur = []
    for r in rs:
        iv = (int(r["elem_start"]), int(r["elem_end"]))
        if cur and iv_same((int(cur[0]["elem_start"]), int(cur[0]["elem_end"])), iv):
            cur.append(r)
        else:
            if cur:
                sites.append(cur)
            cur = [r]
    if cur:
        sites.append(cur)
print("CANDIDATE physical site groups (reciprocal overlap >= %.2f): %d"
      % (OVL, len(sites)))
print("  grouped hits, NOT verified insertion events")

amb = sum(1 for s in sites if len({x["cds_cluster"] for x in s}) > 1)
print("  of which family-AMBIGUOUS (>1 element family hit the same interval): "
      "%d (%.1f%%)" % (amb, 100.0 * amb / max(1, len(sites))))

# ---- host loci: group physical sites by flank context, keeping the distance
def fk(r):
    lf, rf = r.get("left_flank") or "", r.get("right_flank") or ""
    return (lf[-60:].upper(), rf[:60].upper()) if len(lf) >= 60 and len(rf) >= 60 else None


exact = collections.defaultdict(list)
for s in sites:
    # the boundary representative, not an arbitrary first record
    b = max(s, key=lambda x: (float(x["pct_ident"]), float(x["q_cov"])))
    k = fk(b)
    if k:
        exact[k].append(s)
print("EXACT flank-sequence key classes: %d" % len(exact))
print("  NOT a count of homologous host loci: a SNP in the 60bp window splits")
print("  one locus in two, and a repeated host sequence can give two distinct")
print("  loci the same key.")


cols = ["site_id", "genome", "contig", "start", "end", "elem_len",
        "family", "family_ambiguous", "all_seeds", "n_seeds",
        "best_ident", "best_qcov", "flank_complete", "tsd_len",
        "host_locus_key", "n_sites_sharing_host_locus", "boundary_spread_bp",
        "left_flank60", "right_flank60"]
with open(OUT + ".tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for i, s in enumerate(sites):
        fams = sorted({x["cds_cluster"] for x in s})
        # BOUNDARY REPRESENTATIVE, stated explicitly. Hits merged into one
        # site have DIFFERENT boundaries (1826030 vs 1826032, 1827318 vs
        # 1827323), so taking s[0] let flank, TSD and flank_complete depend on
        # record order. The highest-identity hit, tie-broken on query
        # coverage, defines the boundary; boundary_spread_bp reports how much
        # the group disagreed, so an unresolved boundary stays visible rather
        # than being hidden by the choice.
        best = max(s, key=lambda x: (float(x["pct_ident"]), float(x["q_cov"])))
        starts = sorted({int(x["elem_start"]) for x in s})
        ends = sorted({int(x["elem_end"]) for x in s})
        bspread = (starts[-1] - starts[0]) + (ends[-1] - ends[0])
        k = fk(best)
        hk = "%s_%s" % (k[0][:12], k[1][:12]) if k else "."
        fh.write("\t".join(str(v) for v in [
            "S%06d" % i, best["genome"], best["contig"],
            best["elem_start"], best["elem_end"],
            int(best["elem_end"]) - int(best["elem_start"]),
            fams[0] if len(fams) == 1 else "AMBIGUOUS",
            "no" if len(fams) == 1 else "yes",
            ",".join(fams), len(fams),
            best["pct_ident"], best["q_cov"],
            best.get("flank_complete", "."), best.get("tsd_len", "."),
            hk, len(exact.get(k, [])) if k else 0, bspread,
            (best.get("left_flank") or "")[-60:],
            (best.get("right_flank") or "")[:60]]) + "\n")
print("wrote %s.tsv" % OUT)

# ---- the six-family case, both sides, as asked -------------------------
print()
print("THE MULTI-FAMILY CASE, examined from both sides")
cand = [s for s in sites if len({x["cds_cluster"] for x in s}) >= 5]
print("  physical sites hit by >=5 element families: %d" % len(cand))
if cand:
    s = max(cand, key=lambda x: len({y["cds_cluster"] for y in x}))
    b0 = max(s, key=lambda x: (float(x["pct_ident"]), float(x["q_cov"])))
    st0 = sorted({int(x["elem_start"]) for x in s})
    en0 = sorted({int(x["elem_end"]) for x in s})
    print("  boundary spread across the group: start %d-%d, end %d-%d"
          % (st0[0], st0[-1], en0[0], en0[-1]))
    print("  example %s:%s %s-%s"
          % (b0["genome"], b0["contig"], b0["elem_start"], b0["elem_end"]))
    print("  %-12s %-20s %9s %9s %10s %10s"
          % ("family", "source_event", "ident", "q_cov", "start", "end"))
    for x in sorted(s, key=lambda y: -float(y["pct_ident"]))[:8]:
        print("  %-12s %-20s %9s %9s %10s %10s"
              % (x["cds_cluster"], x["source_event"][:20], x["pct_ident"],
                 x["q_cov"], x["elem_start"], x["elem_end"]))
    print()
    print("  identity and coverage decide whether these are redundant seeds of")
    print("  one family or distinct families sharing a homologous region. High")
    print("  identity AND high coverage on all of them means redundant seeds;")
    print("  high identity with PARTIAL coverage means a shared sub-region.")
    fams = sorted({x["cds_cluster"] for x in s})
    rep = collections.Counter()
    for s2 in sites:
        f2 = {x["cds_cluster"] for x in s2}
        if len(f2 & set(fams)) >= 2:
            rep[len(f2 & set(fams))] += 1
    print()
    print("  does this co-hit pattern repeat elsewhere?")
    for k in sorted(rep, reverse=True)[:5]:
        print("    %d of these families co-hit at %d other sites" % (k, rep[k]))

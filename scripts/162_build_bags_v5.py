#!/usr/bin/env python3
"""bags_v5 = observed sites (v4) + sites reconstructed from element COPIES.

v4 holds one site per EVENT: 78,120 of them. Every additional copy of those
same elements -- 783,577 at distinct loci -- contributed nothing, because a
bag site needs the pre-insertion sequence and a homology hit gives only
presence. 160_ excises those copies to recover the site; this merges the two
into one corpus.

THE TWO CLASSES STAY DISTINGUISHABLE, which is the point of building it:

  empty_site_source = observed            the empty allele was seen in a
                                          genome that lacks the insertion
  empty_site_source = tsd_reconstructed   the site was recovered by excising
                                          the element from a carrier

Measured agreement between the two methods, where both are available: 84.7%
of reconstructed flanks match the observed one (152_). That figure is a FLOOR
on accuracy and a CEILING on what can be measured -- it was computed on loci
that HAVE an observed empty allele, i.e. sites demonstrably once empty, while
the copies here have none and may be ancestral. No gate was found that
predicts correctness: three candidates all proved circular.

So the corpus carries both and labels them. Training on observed-only,
reconstructed-only and combined answers the question this analysis could not.

DEDUP IS ACROSS BOTH CLASSES. The site key is (bag_id, flank), the same
independence rule that cut v3's inflated bag counts by 28%. A reconstructed
site landing on a flank an observed site already occupies is NOT a second
observation, and OBSERVED WINS -- a real empty allele is better evidence than
an excision, so the reconstructed duplicate is dropped rather than the
observed one.
"""
import argparse
import collections
import csv
import glob
import json
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI
from lib.util import log                                          # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v4-sites", required=True)
    ap.add_argument("--copysites", required=True, help="glob of *_sites.tsv")
    ap.add_argument("--db", required=True, help="database insertions.tsv")
    ap.add_argument("--cds", required=True, help="insert_cds.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-sites", type=int, default=1)
    a = ap.parse_args()

    # bag + non-coding payload come from the SOURCE EVENT: a copy is the same
    # element, so it inherits the element's cluster and nc regions. No new
    # prodigal or mmseqs run is needed or wanted -- re-predicting ORFs per copy
    # would let the same element land in different bags by chance.
    cds = {r["db_id"]: r for r in csv.DictReader(open(a.cds), delimiter="\t")}
    dbr = {r["db_id"]: r for r in csv.DictReader(open(a.db), delimiter="\t")}
    log("cds records %d, db rows %d" % (len(cds), len(dbr)))

    out, seen = [], {}
    nc_of = {}
    n_obs = 0
    for line in open(a.v4_sites):
        line = line.strip()
        if not line:
            continue
        d = json.loads(line)
        k = (d.get("bag_id"), d.get("flank"))
        if k in seen:
            continue
        seen[k] = "observed"
        d["empty_site_source"] = d.get("empty_site_source", "observed")
        d["n_alleles_at_site"] = d.get("n_alleles_at_site", 1)
        out.append(d)
        n_obs += 1
        # index by source event so a copy can inherit the element's
        # non-coding payload. A copy IS the same element, so re-deriving nc
        # per copy would be wrong as well as wasteful -- and leaving it empty
        # would hand the model a feature that separates the two classes for a
        # reason unrelated to site quality.
        sid = str(d.get("site_id", ""))
        if sid.endswith(".site"):
            nc_of[sid[:-5]] = (d.get("noncoding_regions", []),
                               d.get("nc_sequence_hash", "."),
                               d.get("nc_region_count", 0),
                               d.get("nc_total_len", 0))
    log("observed sites carried over: %d" % n_obs)

    st = collections.Counter()
    # --copysites may carry several space-separated patterns: whole-species
    # files plus per-shard files from the big species.
    files = []
    for pat in a.copysites.split():
        files.extend(sorted(glob.glob(pat)))
    files = sorted(set(files))
    log("copysite files: %d" % len(files))
    for p in files:
        for r in csv.DictReader(open(p), delimiter="\t"):
            st["read"] += 1
            src = r.get("source_event", "")
            c = cds.get(src)
            db = dbr.get(src)
            if not c or not db:
                st["no_source_event"] += 1
                continue
            bag = c.get("cds_cluster_id", "")
            if not bag or bag == ".":
                st["no_cds_cluster"] += 1
                continue
            if src not in nc_of:
                # the source event never produced an observed site (it was
                # dropped for no_noncoding_region or window fit), so there is
                # no nc payload to inherit and the copy would be a different
                # kind of record from every other site in its bag.
                st["source_has_no_site"] += 1
                continue
            fl = r.get("flank", "")
            if len(fl) != 120:
                st["bad_flank_len"] += 1
                continue
            k = (bag, fl)
            if k in seen:
                st["dup_%s" % seen[k]] += 1
                continue
            seen[k] = "reconstructed"
            st["EMITTED"] += 1
            out.append({
                "site_id": r["db_id"],
                "bag_id": bag,
                "cds_cluster_id": bag,
                "corpus": "fna_ins_discovery",
                "species": r.get("species", "."),
                "phylum": db.get("phylum", "."),
                "flank": fl,
                "flank_len": 120,
                "flank_side": "joined",
                "flank_spacer": "",
                "flank_source": "excised_from_carrier",
                "empty_site_source": "tsd_reconstructed",
                "empty_site_n_carrier_genomes": 0,
                "insertion_point_in_flank": 60,
                "orient": r.get("orient", "."),
                "orient_granularity": "per_site",
                "orient_source": "canonical_lexicographic",
                "insert_md5": r.get("insert_md5", "."),
                "inserted_len": r.get("inserted_len", "."),
                "source_event": src,
                "carrier": r.get("carrier", "."),
                "tsd_len": r.get("tsd_len", "."),
                "copy_pct_ident": r.get("pct_ident", "."),
                # nc payload inherited from the source element's own site
                "noncoding_regions": nc_of.get(src, ([], ".", 0, 0))[0],
                "nc_sequence_hash": nc_of.get(src, ([], ".", 0, 0))[1],
                "nc_region_count": nc_of.get(src, ([], ".", 0, 0))[2],
                "nc_total_len": nc_of.get(src, ([], ".", 0, 0))[3],
                "n_alleles_at_site": 1,
                "data_source": "fna_ins_discovery",
                "generator_version_or_commit": "162_build_bags_v5.py",
            })

    # bags below --min-sites are dropped AFTER the merge, so a bag that only
    # reaches the floor thanks to reconstructed sites is kept
    bysize = collections.Counter(s["bag_id"] for s in out)
    out = [s for s in out if bysize[s["bag_id"]] >= a.min_sites]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out + "_sites.jsonl", "w") as fh:
        for s in out:
            fh.write(json.dumps(s) + "\n")

    bags = collections.defaultdict(list)
    for s in out:
        bags[s["bag_id"]].append(s)
    with open(a.out + "_bags.tsv", "w") as fh:
        fh.write("bag_id\tn_sites\tn_observed\tn_reconstructed\tn_species\n")
        for b, ss in sorted(bags.items(), key=lambda kv: -len(kv[1])):
            o = sum(1 for x in ss if x["empty_site_source"] == "observed")
            fh.write("%s\t%d\t%d\t%d\t%d\n"
                     % (b, len(ss), o, len(ss) - o,
                        len({x["species"] for x in ss})))

    nr = sum(1 for s in out if s["empty_site_source"] == "tsd_reconstructed")
    log("=" * 74)
    log("BAGS v5")
    log("  sites total            %8d" % len(out))
    log("    observed             %8d" % (len(out) - nr))
    log("    tsd_reconstructed    %8d" % nr)
    log("  bags                   %8d" % len(bags))
    for k in (5, 10):
        log("  bags with >=%-2d sites   %8d   (observed-only: %d)"
            % (k, sum(1 for ss in bags.values() if len(ss) >= k),
               sum(1 for ss in bags.values()
                   if sum(1 for x in ss
                          if x["empty_site_source"] == "observed") >= k)))
    log("")
    log("  copysite rows read     %8d" % st["read"])
    for k in ("EMITTED", "dup_observed", "dup_reconstructed",
              "no_source_event", "no_cds_cluster", "source_has_no_site",
              "bad_flank_len"):
        if st[k]:
            log("    %-22s %8d" % (k, st[k]))
    log("")
    log("  `dup_observed` are reconstructed sites landing on a flank an")
    log("  observed site already holds. Observed wins -- a real empty allele")
    log("  outranks an excision.")
    log("=" * 74)
    log("wrote %s_sites.jsonl and %s_bags.tsv" % (a.out, a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())

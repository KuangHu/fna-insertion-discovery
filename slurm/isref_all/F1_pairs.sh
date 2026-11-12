#!/bin/bash
#SBATCH --job-name=isall_F1
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
# STAGE F1 -- choose the genome pairs ISCompare will run on.
#
# "Suitable pairs within each species" are not arbitrary pairs and certainly
# not all-vs-all: insertion_db_v1 already holds, for every occurrence, a
# carrier and a reference whose locus was verified empty.  Those ARE the
# suitable pairs, already validated, so this reuses finished work instead of
# re-deriving it.  Greedy set cover then picks the fewest pairs that witness
# every bag in a species, because one good pair witnesses many bags at once
# (the abaum pilot pair carried 54 differentially located ISs).
set -euo pipefail
/global/home/users/kh36969/.conda/envs/claude-env/bin/python3 - <<'PY'
import csv, collections, json
D="/global/scratch/users/kh36969/fna_ins_discovery/insertion_db_v1"
C="/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all"
V="/global/scratch/users/kh36969/fna_ins_discovery/isref_v1"
occbag={}; occsp={}
for r in csv.DictReader(open(C+"/bags/strict_occurrence_bags.tsv"),delimiter="\t"):
    occbag[r["occurrence_id"]]=r["cds_cluster_id"]; occsp[r["occurrence_id"]]=r["species"]
have=set()
for r in csv.DictReader(open(V+"/corpus_v2.tsv"),delimiter="\t"): have.add(r["assembly"])

pairs=collections.defaultdict(set)      # (species, carrier, ref) -> bags witnessed
nocc=skip_noref=skip_missing=0
for r in csv.DictReader(open(D+"/occurrences_browse.tsv"),delimiter="\t"):
    oid=r["occurrence_id"]; b=occbag.get(oid)
    if not b: continue
    nocc+=1
    q=r["carrier_assembly_accession_version"]; ref=r["reference_assembly_accession_version"]
    if not ref or ref=="": skip_noref+=1; continue
    if q not in have or ref not in have: skip_missing+=1; continue
    pairs[(occsp[oid],q,ref)].add(b)
print(json.dumps({"bagged_occurrences_scanned":nocc,
  "skipped_no_reference":skip_noref,"skipped_genome_not_in_corpus":skip_missing,
  "candidate_pairs":len(pairs)},indent=2))

bysp=collections.defaultdict(list)
for (sp,q,ref),bs in pairs.items(): bysp[sp].append((q,ref,bs))
BUDGET=int(__import__("os").environ.get("BUDGET","40"))
chosen=[]; report={}
for sp,cands in sorted(bysp.items()):
    need=set().union(*[c[2] for c in cands])
    cov=set(); picked=[]
    pool=sorted(cands,key=lambda c:(-len(c[2]),c[0],c[1]))
    while pool and len(picked)<BUDGET and cov!=need:
        pool.sort(key=lambda c:-len(c[2]-cov))
        q,ref,bs=pool.pop(0)
        if not (bs-cov): break
        picked.append((q,ref,bs)); cov|=bs
    for q,ref,bs in picked: chosen.append((sp,q,ref,len(bs)))
    report[sp]={"bags_in_species":len(need),"pairs_available":len(cands),
                "pairs_chosen":len(picked),"bags_covered":len(cov),
                "bags_uncovered":len(need-cov)}
with open(W+"/iscompare_pairs.tsv","w") as o:
    o.write("species\tquery_assembly\tref_assembly\tbags_witnessed\n")
    for t in sorted(chosen): o.write("%s\t%s\t%s\t%d\n"%t)
print("\nper species:")
print("%-16s %8s %10s %8s %9s %10s"%("species","bags","pairs_avail","chosen","covered","uncovered"))
for sp,r in sorted(report.items()):
    print("%-16s %8d %10d %8d %9d %10d"%(sp,r["bags_in_species"],r["pairs_available"],
          r["pairs_chosen"],r["bags_covered"],r["bags_uncovered"]))
print("\ntotal pairs chosen:",len(chosen))
PY

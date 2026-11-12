#!/bin/bash
#SBATCH --job-name=isref_loci
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
/global/home/users/kh36969/.conda/envs/claude-env/bin/python3 - <<'PY'
import csv, collections, json
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_v1"
# 181 queries are not 181 different elements: many are rve / DDE_Tnp_1
# variants, so one physical copy on a contig is hit by several queries and
# counted several times by 233.  Collapse to genomic LOCI before quoting any
# number as "how many copies exist".
loci=collections.defaultdict(list)
n=0
for r in csv.DictReader(open(W+"/isref_copies.tsv"),delimiter="\t"):
    if r["copy_kind"]!="full_copy": continue
    n+=1
    loci[(r["assembly"],r["contig"])].append(
        (int(r["s0"]),int(r["s1"]),r["bag_id"],r["catalogue_status"],r["species"]))
merged=[]
for (asm,ctg),v in loci.items():
    v.sort()
    cs,ce,bags,stat,sp = v[0][0],v[0][1],{v[0][2]},{v[0][3]},v[0][4]
    for s0,s1,b,st,s in v[1:]:
        if s0 < ce - 50:                      # same physical copy
            ce=max(ce,s1); bags.add(b); stat.add(st)
        else:
            merged.append((asm,ctg,cs,ce,bags,stat,sp)); cs,ce,bags,stat = s0,s1,{b},{st}
    merged.append((asm,ctg,cs,ce,bags,stat,sp))
cat=sum(1 for m in merged if "catalogued" in m[5])
gen=len({m[0] for m in merged})
spp=collections.Counter(m[6] for m in merged)
multi=sum(1 for m in merged if len(m[4])>1)
print(json.dumps({
 "full_copy_rows_from_233": n,
 "distinct_genomic_loci": len(merged),
 "loci_hit_by_more_than_one_bag_query": multi,
 "loci_overlapping_a_catalogued_occurrence": cat,
 "loci_not_in_the_catalogue": len(merged)-cat,
 "genomes_carrying_at_least_one": gen,
 "corpus_genomes": 40541}, indent=2))
print("\nloci by species:")
for s,c in spp.most_common(): print("  %-16s %7d"%(s,c))
with open(W+"/isref_loci.tsv","w") as o:
    o.write("assembly\tspecies\tcontig\ts0\ts1\tn_bag_queries\tbags\tcatalogued\n")
    for asm,ctg,cs,ce,bags,stat,sp in sorted(merged):
        o.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\n"%(asm,sp,ctg,cs,ce,len(bags),
                ";".join(sorted(bags)),int("catalogued" in stat)))
PY
echo "=== BAGS000173 specifically ==="
awk -F'\t' 'NR>1 && $7 ~ /BAGS000173/ {n++; g[$1]=1; c+=$8} END{
  printf "  loci %d in %d genomes; %d overlap a catalogued occurrence\n", n, length(g), c}' \
  /global/scratch/users/kh36969/fna_ins_discovery/isref_v1/isref_loci.tsv

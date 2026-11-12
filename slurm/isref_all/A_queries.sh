#!/bin/bash
#SBATCH --job-name=isall_A
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
# STAGE A -- build the query set for ALL bags.
#
# Not one sequence per bag.  Bag membership was decided on dominant-CDS PROTEIN
# similarity while the search runs on DNA, and in the 181-reference run 53 of
# 181 references collapsed into 24 DNA-level components.  So the query set is
# built by grouping the insertion DNA itself:
#   * a bag whose members are diverse yields several groups, hence several
#     representatives -- the diversity requirement;
#   * two bags holding the same element land in one group, so the same physical
#     copy is not counted twice -- dedup applied at the query level rather than
#     patched up afterwards.
#
# WHY NOT mmseqs.  `easy-cluster` on this nucleotide set emitted 119 internal
# errors ("deleted from ... now is empty, but not assigned to a cluster") and
# mmseqs 18 has no --search-type switch on easy-cluster to force DNA mode.  A
# clustering that reports internal errors cannot be trusted as the backbone of
# everything downstream, so this uses blastn all-vs-all with the SAME criterion
# already used on the 181 references (>=90% identity over >=80% of the shorter
# sequence, reciprocal) -- fully controlled, and comparable to that result.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
[ -s "$W/queries.fna" ] && { echo "stage A already done"; exit 0; }
mkdir -p "$W/tmp"

if [ ! -s "$W/selfaln.tsv" ]; then
  "$BL/makeblastdb" -in "$C/dedup_insertions.fna" -dbtype nucl -out "$W/tmp/alldb" -title all >/dev/null
  "$BL/blastn" -task megablast -query "$C/dedup_insertions.fna" -db "$W/tmp/alldb" \
     -evalue 1e-50 -perc_identity 90 -max_target_seqs 100000 -num_threads 64 \
     -outfmt '6 qseqid sseqid pident length qlen slen' -out "$W/selfaln.tsv" 2>"$W/tmp/selfaln.log"
fi
echo "self-alignment rows: $(wc -l < "$W/selfaln.tsv")"

"$PY" - <<'PY'
import csv, collections, json
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all"
C="/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"
keybags=collections.defaultdict(set); bags=collections.defaultdict(set)
for r in csv.DictReader(open(C+"/bags/strict_occurrence_bags.tsv"),delimiter="\t"):
    keybags[r["seq_key"]].add(r["cds_cluster_id"]); bags[r["cds_cluster_id"]].add(r["seq_key"])
meta={r["bag_id"]:r for r in csv.DictReader(open(C+"/bags/strict_bags.tsv"),delimiter="\t")}
print("bags %d, bagged seq_keys %d"%(len(bags),len(keybags)))

seqs={}; n=None; buf=[]
def fl():
    if n: seqs[n.split()[0].replace("SEQ|","")]="".join(buf)
for line in open(C+"/dedup_insertions.fna"):
    if line[0]==">": fl(); n=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
fl()
print("insertion sequences:",len(seqs))

# reciprocal >=90% identity over >=80% of the shorter sequence
fwd=collections.defaultdict(set)
for l in open(W+"/selfaln.tsv"):
    q,s,pid,aln,ql,sl=l.rstrip("\n").split("\t")
    q=q.replace("SEQ|",""); s=s.replace("SEQ|","")
    if q==s: continue
    if int(aln)>=0.80*min(int(ql),int(sl)) and float(pid)>=90: fwd[q].add(s)
adj=collections.defaultdict(set)
for q,ss in fwd.items():
    for s in ss:
        if q in fwd.get(s,()): adj[q].add(s); adj[s].add(q)

comp_of={}; comps=[]
for k in seqs:
    if k in comp_of: continue
    stack=[k]; c=set()
    while stack:
        x=stack.pop()
        if x in c: continue
        c.add(x); stack.extend(adj[x]-c)
    i=len(comps); comps.append(c)
    for x in c: comp_of[x]=i
print("DNA groups over all insertions:",len(comps))

# a query = one DNA group holding at least one bagged insertion; its
# representative is the LONGEST member, so a partial copy never speaks for the
# group and shorten the search
cl=collections.defaultdict(set)
for k,b in keybags.items(): cl[comp_of[k]] |= b
rows=[]
for i,bs in sorted(cl.items()):
    members=[m for m in comps[i] if m in seqs]
    rep=max(members,key=lambda m:(len(seqs[m]),m))
    sites=sum(int(meta[b]["n_independent_sites"]) for b in bs if b in meta)
    rows.append((rep,sorted(bs),len(seqs[rep]),sites,len(members)))
with open(W+"/queries.fna","w") as fa, open(W+"/queries.tsv","w") as t:
    t.write("query_id\tn_bags\tbags\tlen_nt\tcatalogue_sites\tgroup_members\n")
    for rep,bs,L,sites,nm in rows:
        fa.write(">Q_%s\n%s\n"%(rep,seqs[rep]))
        t.write("Q_%s\t%d\t%s\t%d\t%d\t%d\n"%(rep,len(bs),";".join(bs),L,sites,nm))
perbag=collections.Counter()
for rep,bs,L,s,nm in rows:
    for b in bs: perbag[b]+=1
uncov=[b for b in bags if b not in perbag]
print(json.dumps({
 "queries":len(rows),
 "total_query_bp":sum(r[2] for r in rows),
 "bags_total":len(bags),
 "bags_covered":len(perbag),
 "bags_NOT_covered":len(uncov),
 "bags_with_multiple_representatives":sum(1 for b,n in perbag.items() if n>1),
 "max_representatives_for_one_bag":max(perbag.values()) if perbag else 0,
 "queries_shared_by_more_than_one_bag":sum(1 for r in rows if len(r[1])>1)},indent=2))
if uncov: print("UNCOVERED BAGS (first 10):",uncov[:10])
PY

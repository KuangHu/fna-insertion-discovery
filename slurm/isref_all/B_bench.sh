#!/bin/bash
#SBATCH --job-name=isall_B
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.out
# STAGE B -- measure before committing node-hours.
#
# The 181-reference sweep cost ~2.3 node-hours for 274 kb of query against the
# 174 GB corpus.  The all-bag query set is far larger, and naive scaling would
# put the full sweep in the hundreds of node-hours.  So: time the real query
# set against ONE shard with both dc-megablast (what the pilot used) and
# megablast (faster, less sensitive), and measure what megablast actually
# LOSES in recovered loci rather than assuming it is fine.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
mkdir -p "$W/bench"
SH=${SH:-0}

# rebuild shard SH exactly as the v1 sweep sharded it, so the benchmark is on
# a real slice of the corpus and the numbers extrapolate honestly
sf="$W/bench/shard.fna"
if [ ! -s "$sf" ]; then
  : > "$sf"
  awk -v n=64 -v i="$SH" 'NR>1 && (NR-2)%n==i {print $1"\t"$3}' "$V/corpus_v2.tsv" |
  while IFS=$'\t' read -r acc path; do
    awk -v a="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"a"|"p[1]; next} {print}' "$path" >> "$sf"
  done
fi
echo "shard bytes: $(stat -c%s "$sf")   genomes: $(awk -v n=64 -v i="$SH" 'NR>1 && (NR-2)%n==i' "$V/corpus_v2.tsv" | wc -l)"
echo "query bp: $(grep -v '^>' "$W/queries.fna" | tr -d '\n' | wc -c)   queries: $(grep -c '^>' "$W/queries.fna")"
"$BL/makeblastdb" -in "$sf" -dbtype nucl -out "$W/bench/db" -title bench >/dev/null 2>&1

for task in dc-megablast megablast; do
  s=$(date +%s)
  "$BL/blastn" -task $task -query "$W/queries.fna" -db "$W/bench/db" \
     -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads 64 \
     -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
     -out "$W/bench/$task.raw" 2>"$W/bench/$task.log"
  e=$(date +%s)
  awk -F'\t' '$4>=200' "$W/bench/$task.raw" > "$W/bench/$task.tsv"
  echo "$task: $((e-s)) s, $(wc -l < "$W/bench/$task.tsv") hsps"
  echo "$task $((e-s))" >> "$W/bench/timing.txt"
done

echo "=== what megablast loses, at locus level ==="
/global/home/users/kh36969/.conda/envs/claude-env/bin/python3 - <<'PY'
import collections
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all/bench"
def loci(p):
    byc=collections.defaultdict(list)
    for l in open(p):
        f=l.split("\t")
        if int(f[3])<200: continue
        ss,se=int(f[8]),int(f[9]); 
        if se<ss: ss,se=se,ss
        byc[f[1]].append((ss,se))
    out=set()
    for c,v in byc.items():
        v.sort(); cs,ce=v[0]
        for s,e in v[1:]:
            if s<ce-50: ce=max(ce,e)
            else: out.add((c,cs//500)); cs,ce=s,e
        out.add((c,cs//500))
    return out
d=loci(W+"/dc-megablast.tsv"); m=loci(W+"/megablast.tsv")
print("dc-megablast loci: %d"%len(d))
print("megablast loci:    %d"%len(m))
print("found by dc only:  %d (%.2f%% of dc)"%(len(d-m),100*len(d-m)/max(1,len(d))))
print("found by mb only:  %d"%len(m-d))
t={}
for l in open(W+"/timing.txt"):
    k,v=l.split(); t[k]=int(v)
if len(t)==2:
    print("\nspeedup megablast vs dc-megablast: %.1fx"%(t["dc-megablast"]/max(1,t["megablast"])))
    for k,v in t.items():
        print("  full 64-shard sweep with %-14s ~%.1f node-hours"%(k, v*64/3600))
PY

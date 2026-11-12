#!/bin/bash
#SBATCH --job-name=isall_B2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
# STAGE B2 -- what does the all-sequence query set actually cost, and what
# would the cheaper choices lose?
#
# Naive scaling from the completed 181-query sweep says 782 node-hours for
# 63,095 queries.  That number is almost certainly wrong, because a large part
# of a blast run is building and scanning the genome database and that does
# not scale with query size at all.  Wrong in which direction and by how much
# is the whole question, so this measures it instead of arguing about it.
#
# Four runs over ONE stratified genome sample, same database, 64 threads each:
#   all x dc-megablast   all x megablast
#   rep x dc-megablast   rep x megablast
# Timing gives the projection; the locus comparison gives the cost of being
# cheap.  Nothing is swept until this reports.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
B=$W/bench; mkdir -p $B
export PATH="$BL:$PATH"

NG=64   # genomes in the sample
if [ ! -s "$B/sample.fna" ]; then
  # stratified: proportional across the 19 species, deterministic
  awk -F'\t' -v n=$NG 'NR>1{c[$2]++; rows[$2"|"c[$2]]=$0; tot++}
    END{ for(s in c){ k=int(n*c[s]/tot+0.5); if(k<1)k=1;
           for(i=1;i<=k && i<=c[s];i++){ step=int(c[s]/k); if(step<1)step=1;
             print rows[s"|"((i-1)*step+1)] } } }' "$V/corpus_v2.tsv" > "$B/sample.tsv"
  : > "$B/sample.fna"
  while IFS=$'\t' read -r acc sp path; do
    awk -v a="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"a"|"p[1]; next}{print}' "$path" >> "$B/sample.fna"
  done < "$B/sample.tsv"
fi
echo "sample genomes: $(wc -l < "$B/sample.tsv")  bp: $(grep -v '^>' "$B/sample.fna" | tr -d '\n' | wc -c)"
[ -s "$B/db.nsq" ] || makeblastdb -in "$B/sample.fna" -dbtype nucl -out "$B/db" -title s >/dev/null

run() {   # name querfile task
  local nm=$1 q=$2 tk=$3
  [ -s "$B/$nm.tsv" ] && { echo "$nm already done"; return 0; }
  local t0=$(date +%s)
  blastn -task "$tk" -query "$q" -db "$B/db" -evalue 1e-30 -perc_identity 80 \
     -max_target_seqs 1000000 -num_threads 64 \
     -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
     -out "$B/$nm.raw" 2>"$B/$nm.err"
  local t1=$(date +%s)
  awk -F'\t' '$4>=200' "$B/$nm.raw" > "$B/$nm.tsv"
  echo "$nm seconds=$((t1-t0)) hsps_raw=$(wc -l < "$B/$nm.raw") hsps_200=$(wc -l < "$B/$nm.tsv")" \
      | tee -a "$B/timings.txt"
}

run all_dcmega "$W/all_queries.fna" dc-megablast
run all_mega   "$W/all_queries.fna" megablast
run rep_dcmega "$W/queries.fna"     dc-megablast
run rep_mega   "$W/queries.fna"     megablast

echo; echo "=== TIMINGS ==="; cat "$B/timings.txt"
echo; echo "=== LOCUS COMPARISON ==="
for n in all_dcmega all_mega rep_dcmega rep_mega; do
  awk -F'\t' '{split($2,p,"|"); s=$7<$8?$7:$8; e=$7<$8?$8:$7;
               printf "%s\t%s\t%d\n", p[1], p[2], int(s/500)}' "$B/$n.tsv" | sort -u > "$B/$n.loci"
  echo "$n loci=$(wc -l < "$B/$n.loci")"
done
echo
echo "lost by each cheaper option, relative to all x dc-megablast:"
for n in all_mega rep_dcmega rep_mega; do
  lost=$(comm -23 "$B/all_dcmega.loci" "$B/$n.loci" | wc -l)
  gain=$(comm -13 "$B/all_dcmega.loci" "$B/$n.loci" | wc -l)
  tot=$(wc -l < "$B/all_dcmega.loci")
  awk -v n="$n" -v l="$lost" -v g="$gain" -v t="$tot" 'BEGIN{
    printf "  %-12s missing %7d of %7d (%.2f%%), extra %d\n", n, l, t, 100*l/t, g}'
done

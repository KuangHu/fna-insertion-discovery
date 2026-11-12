#!/bin/bash
#SBATCH --job-name=isall_C
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=14:00:00
#SBATCH --array=0-15
#SBATCH --output=%x_%A_%a.out
# STAGE C -- full corpus sweep, ALL 63,095 insertion DNA queries, dc-megablast.
#
# Measured by B2 (job 27721664) on a 66-genome stratified sample, not guessed:
# 822 s on 64 cores -> ~140 node-hours for 40,541 genomes.  megablast would cut
# that to 75 but miss 4.35% of loci, which are exactly the diverged copies
# dc-megablast exists to find, so it is not used.  Every bag is included,
# questionable ones among them.
#
# TWO FAILURES FROM THE EARLIER SWEEP ARE GUARDED HERE.
# 1. A scratch I/O fault made shards short while still writing a hits file, so
#    a shard that cannot read every genome it was assigned writes a FAILED
#    marker and NO hits file, rather than a quietly incomplete one.
# 2. B2 showed the raw HSPs come to ~281 GB for the corpus.  Chaining to copies
#    happens per shard, inside the job, and the raw HSPs are gzipped; the full
#    281 GB is never staged as an intermediate.
#
# Resumable at shard grain: a shard is skipped only when its copies file AND
# its .done marker exist.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
D=/global/scratch/users/kh36969/fna_ins_discovery/insertion_db_v1
S=/global/home/users/kh36969/fna_based_mgefinder_project/scripts
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
NSHARD=${NSHARD:-64}
THR=64
head -c 100 "$V/corpus_v2.tsv" >/dev/null || { echo "FATAL: scratch unreadable on $(hostname)"; exit 1; }

run_shard() {
  local i=$1
  [ -s "$W/copies/shard_$i.tsv" ] && [ -s "$W/shardlog/done_$i" ] && { echo "shard $i already done"; return 0; }
  local t0=$(date +%s)
  local sf="${TMPDIR:-/tmp}/shard_$i.fna"
  local want got acc path
  want=$(awk -v n="$NSHARD" -v k="$i" 'NR>1 && (NR-2)%n==k' "$V/corpus_v2.tsv" | wc -l)
  : > "$sf"; got=0
  while IFS=$'\t' read -r acc path; do
    if ! awk -v a="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"a"|"p[1]; next}{print}' \
            "$path" >> "$sf" 2>/dev/null; then
      echo "UNREADABLE $acc $path" >> "$W/shardlog/FAILED_$i"; continue
    fi
    got=$((got+1))
  done < <(awk -v n="$NSHARD" -v k="$i" 'NR>1 && (NR-2)%n==k {print $1"\t"$3}' "$V/corpus_v2.tsv")
  if [ "$got" -ne "$want" ]; then
    echo "shard $i SHORT: $got/$want genomes -- refusing to write a hits file" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf"; return 1
  fi

  local db="${TMPDIR:-/tmp}/db_$i"
  "$BL/makeblastdb" -in "$sf" -dbtype nucl -out "$db" -title "s$i" >/dev/null 2>&1 || {
    echo "shard $i makeblastdb failed" | tee -a "$W/shardlog/FAILED_$i"; rm -f "$sf"; return 1; }
  if ! "$BL/blastn" -task dc-megablast -query "$W/all_queries.fna" -db "$db" \
        -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads $THR \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -out "$W/hits/shard_$i.raw" 2>"$W/shardlog/blast_$i.err"; then
    echo "shard $i blastn failed" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$W/hits/shard_$i.raw" "$sf" "$db".*; return 1
  fi
  awk -F'\t' '$4>=200' "$W/hits/shard_$i.raw" > "$W/hits/shard_$i.tsv"

  # chain to copies immediately, then drop the raw HSPs to gzip
  if ! "$PY" -I "$S/239_chain_copies.py" --hits "$W/hits/shard_$i.tsv" \
        --queries "$W/all_queries.tsv" --corpus "$V/corpus_v2.tsv" \
        --occurrences "$D/occurrences_browse.tsv" \
        --occbags "$C/bags/strict_occurrence_bags.tsv" \
        --out "$W/copies/shard_$i" > "$W/shardlog/chain_$i.log" 2>&1; then
    echo "shard $i chaining failed" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf" "$db".*; return 1
  fi
  gzip -f "$W/hits/shard_$i.tsv"; rm -f "$W/hits/shard_$i.raw" "$sf" "$db".*
  echo "$(($(date +%s)-t0))" > "$W/shardlog/done_$i"
  echo "shard $i done in $(cat "$W/shardlog/done_$i")s"
}
export -f run_shard; export W V C D S PY BL NSHARD THR

T=${SLURM_ARRAY_TASK_ID:-0}; NT=${SLURM_ARRAY_TASK_COUNT:-1}
for i in $(seq 0 $((NSHARD-1))); do
  [ $((i % NT)) -eq "$T" ] || continue
  run_shard "$i" || echo "shard $i returned non-zero (recorded, continuing)"
done
echo "task $T finished; shards done: $(ls "$W/shardlog"/done_* 2>/dev/null | wc -l)/$NSHARD; failed markers: $(ls "$W/shardlog"/FAILED_* 2>/dev/null | wc -l)"

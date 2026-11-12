#!/bin/bash
#SBATCH --job-name=isref_sweep2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --array=0-3
#SBATCH --exclude=n0006.cf1,n0008.cf1
#SBATCH --output=%x_%A_%a.out
# Re-run of the shards lost when scratch dropped out under array tasks 1 and 2
# of job 27693571 ("Cannot send after transport endpoint shutdown", 504 and 257
# times).  The first version treated an unreadable genome as `continue`, so an
# I/O fault became a silently short shard that still wrote a hits file -- two
# of those shards reported 0 hsps and two reported none at all.  A short shard
# is worse than a failed one, so this version:
#   * counts the genomes it actually appended and refuses to emit a hits file
#     unless that equals the number of genomes the shard was assigned;
#   * checks the exit status of makeblastdb and blastn instead of trusting the
#     pipeline;
#   * leaves a FAILED_<i> marker so a short shard is visible, not invisible.
# Shards already produced by the clean tasks 0 and 3 are skipped on sight.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
NSHARD=${NSHARD:-64}
PAR=${PAR:-8}
THR=${THR:-8}
mkdir -p "$W/hits" "$W/shardlog" "$W/tmp"

# fail fast if this node cannot see scratch at all
head -c 100 "$W/corpus.tsv" > /dev/null || { echo "FATAL: scratch unreadable on $(hostname)"; exit 1; }

run_shard() {
  local i=$1
  local sf="$W/tmp/s2_$i.fna" db="$W/tmp/s2_$i" out="$W/hits/shard_$i.tsv"
  [ -s "$out" ] && { echo "shard $i already done"; return 0; }
  rm -f "$W/shardlog/FAILED_$i"
  local want got
  want=$(awk -v n="$NSHARD" -v i="$i" 'NR>1 && (NR-2)%n==i' "$W/corpus.tsv" | wc -l)
  : > "$sf"; got=0
  while IFS=$'\t' read -r acc path; do
    if ! awk -v a="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"a"|"p[1]; next} {print}' \
            "$path" >> "$sf" 2>/dev/null; then
      echo "UNREADABLE $acc $path" >> "$W/shardlog/FAILED_$i"; continue
    fi
    got=$((got+1))
  done < <(awk -v n="$NSHARD" -v i="$i" 'NR>1 && (NR-2)%n==i {print $1"\t"$3}' "$W/corpus.tsv")
  if [ "$got" -ne "$want" ]; then
    echo "shard $i SHORT: $got/$want genomes -- refusing to write a hits file" \
      | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf"; return 1
  fi
  if ! "$BL/makeblastdb" -in "$sf" -dbtype nucl -out "$db" -title "s2_$i" \
        > "$W/shardlog/mkdb2_$i.log" 2>&1; then
    echo "shard $i makeblastdb FAILED" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf" "$db".n*; return 1
  fi
  if ! "$BL/blastn" -task dc-megablast -query "$W/is_like_queries.fna" -db "$db" \
        -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads "$THR" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -out "$W/tmp/s2_$i.raw" 2> "$W/shardlog/blast2_$i.log"; then
    echo "shard $i blastn FAILED" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf" "$db".n* "$W/tmp/s2_$i.raw"; return 1
  fi
  awk -F'\t' '$4>=200' "$W/tmp/s2_$i.raw" > "$out.part"
  mv "$out.part" "$out"
  rm -f "$sf" "$db".n* "$W/tmp/s2_$i.raw"
  echo "shard $i done: $want genomes, $(wc -l < "$out") hsps"
}
export -f run_shard; export W BL NSHARD THR

T=${SLURM_ARRAY_TASK_ID:-0}; NT=${SLURM_ARRAY_TASK_COUNT:-1}
# only the shards that are still missing
missing=$(for i in $(seq 0 $((NSHARD-1))); do [ -s "$W/hits/shard_$i.tsv" ] || echo "$i"; done)
echo "$missing" | awk -v n="$NT" -v t="$T" 'NF && (NR-1)%n==t' \
  | xargs -r -P "$PAR" -I{} bash -c 'run_shard {} || true'
echo "array task $T finished"
nfail=$(ls "$W"/shardlog/FAILED_* 2>/dev/null | wc -l)
ndone=$(ls "$W"/hits/shard_*.tsv 2>/dev/null | wc -l)
echo "shards complete: $ndone/$NSHARD, failed markers: $nfail"

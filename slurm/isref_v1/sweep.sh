#!/bin/bash
#SBATCH --job-name=isref_sweep
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out
# Reference-guided copy sweep: blastn the 181 IS-like bag representatives
# against every assembly in the corpus.
#
# This is ISseeker step 1 (blast the IS against a contig collection) done at
# corpus scale with BLAST+ instead of the legacy blastall/formatdb that
# ISseeker needs and that no longer ships.  Flank extraction and the
# insertion-event test are NOT done here -- they are done afterwards against
# the validated insertion_db_v1 junction geometry, which is a better reference
# than ISseeker's single-replicon assumption allows.
#
# Shard fasta and shard BLAST db are deleted as soon as each shard's blastn
# finishes, so peak disk is ~(concurrency x 3.5 GB), not 220 GB.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
NSHARD=${NSHARD:-64}
PAR=${PAR:-16}
THR=${THR:-4}
mkdir -p "$W/hits" "$W/shardlog" "$W/tmp"

# one shard: build a renamed multi-fasta, make the db, blast, drop the db
run_shard() {
  local i=$1
  local sf="$W/tmp/shard_$i.fna" db="$W/tmp/shard_$i" out="$W/hits/shard_$i.tsv"
  [ -s "$out" ] && { echo "shard $i already done"; return 0; }
  : > "$sf"
  awk -v n="$NSHARD" -v i="$i" 'NR>1 && (NR-2)%n==i {print $1"\t"$3}' "$W/corpus.tsv" |
  while IFS=$'\t' read -r acc path; do
    [ -s "$path" ] || { echo "MISSING $acc $path" >> "$W/shardlog/missing_$i.txt"; continue; }
    # prefix every contig id with the assembly so a hit is traceable; keep the
    # original id after the first '|' so coordinates can be replayed.
    awk -v a="$acc" '/^>/{sub(/^>/,""); split($0,p," "); print ">"a"|"p[1]; next} {print}' "$path" >> "$sf"
  done
  [ -s "$sf" ] || { echo "shard $i empty"; rm -f "$sf"; return 0; }
  "$BL/makeblastdb" -in "$sf" -dbtype nucl -out "$db" -title "shard_$i" > "$W/shardlog/mkdb_$i.log" 2>&1
  "$BL/blastn" -task dc-megablast -query "$W/is_like_queries.fna" -db "$db" \
      -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads "$THR" \
      -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
      2> "$W/shardlog/blast_$i.log" \
    | awk -F'\t' '$4>=200' > "$out.part"
  mv "$out.part" "$out"
  rm -f "$sf" "$db".n*
  echo "shard $i done: $(wc -l < "$out") hsps"
}
export -f run_shard; export W BL NSHARD THR

T=${SLURM_ARRAY_TASK_ID:-0}; NT=${SLURM_ARRAY_TASK_COUNT:-1}
seq 0 $((NSHARD-1)) | awk -v n="$NT" -v t="$T" '(NR-1)%n==t' \
  | xargs -P "$PAR" -I{} bash -c 'run_shard {}'
echo "array task $T finished its shards"

# only the last task standing concatenates
done_n=$(ls "$W"/hits/shard_*.tsv 2>/dev/null | wc -l)
if [ "$done_n" -lt "$NSHARD" ]; then
  echo "shards complete: $done_n/$NSHARD -- not concatenating yet"; exit 0
fi

cat "$W"/hits/shard_*.tsv > "$W/is_like_hits.tsv"
echo "TOTAL HSPs: $(wc -l < "$W/is_like_hits.tsv")"
cut -f1 "$W/is_like_hits.tsv" | sort | uniq -c | sort -rn | head -20
( cat "$W"/shardlog/missing_*.txt 2>/dev/null || true ) | wc -l | sed "s/^/missing genome files: /"

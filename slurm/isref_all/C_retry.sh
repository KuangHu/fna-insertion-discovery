#!/bin/bash
#SBATCH --job-name=isall_Cr
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --array=0-8
#SBATCH --output=%x_%A_%a.out
# STAGE C RETRY -- redo the shards job 27725308 got wrong, two different ways.
#
# WHAT WENT WRONG, AND WHY THE ORIGINAL GUARDS DID NOT CATCH IT.
# Shards 6, 9, 11, 14 wrote a `done` marker over partial data.  blastn read
# only ~34% of all_queries.fna (21,486 of 63,095 distinct query ids reached
# the output of shard 9, against 61,772 for a healthy shard), treated the
# short read as end-of-input, and EXITED 0 WITH EMPTY STDERR.  The original
# `if ! blastn` guard tests the exit status, so it saw nothing wrong.  The
# database was intact in every case -- 634 distinct assemblies in good and bad
# shards alike -- so this was the query side, not the subject side.
#
# Cause: 16 nodes re-read the same 95 MB query FASTA off Lustre for every one
# of their four shards, and a contended read came back short.  Two changes:
#   1. the query file is staged to node-local disk ONCE per task and verified
#      by byte count and record count before any blast runs, so the shared
#      filesystem is out of the hot path.
#   2. a post-blast guard counts distinct query ids in the output.  Exit
#      status cannot detect a silent short read; only the content can.
#
# The unguarded `awk '$4>=200'` filter is also checked now.  set -e does not
# apply inside run_shard (it is called in a `||` context), so a failing filter
# would previously have let chaining proceed over a truncated file.
#
# TASK MAP
#   tasks 0-7 : full re-blast of shards 6 9 11 14 33 42 43 49
#               6 9 11 14 -- silent short query read, done marker overwritten
#               33 49     -- transient scratch fault reading one genome; both
#                            genomes verified readable again before this ran
#               42 43     -- Errno 5 during chaining, raw truncated mid-line
#   task 8    : shard 47 CHAIN ONLY.  Its blast finished and its filtered hits
#               are intact -- chaining read all 55,310,096 HSPs and then died
#               writing output (BrokenPipeError 108).  Re-blasting it would
#               burn 3.3 node-hours to reproduce a file that is already on
#               disk, so it is not re-blasted.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
D=/global/scratch/users/kh36969/fna_ins_discovery/insertion_db_v1
S=/global/home/users/kh36969/fna_based_mgefinder_project/scripts
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
NSHARD=64
THR=64
MINQ=55000          # distinct query ids a healthy shard must reach; good shards hit ~61,772
REDO=(6 9 11 14 33 42 43 49)
T=${SLURM_ARRAY_TASK_ID:-0}

# ---- stage the query file node-locally and verify it ------------------------
TD=${TMPDIR:-/tmp}
QSRC="$W/all_queries.fna"
QL="$TD/all_queries.fna"
want_b=$(stat -c%s "$QSRC")
want_n=$(grep -c '^>' "$QSRC")
cp "$QSRC" "$QL"
got_b=$(stat -c%s "$QL")
got_n=$(grep -c '^>' "$QL")
if [ "$got_b" -ne "$want_b" ] || [ "$got_n" -ne "$want_n" ]; then
  echo "FATAL: query stage mismatch on $(hostname): bytes $got_b/$want_b records $got_n/$want_n"; exit 1
fi
echo "query file staged on $(hostname): $got_n records, $got_b bytes"

chain_shard() {
  local i=$1
  "$PY" -I "$S/239_chain_copies.py" --hits "$W/hits/shard_$i.tsv" \
      --queries "$W/all_queries.tsv" --corpus "$V/corpus_v2.tsv" \
      --occurrences "$D/occurrences_browse.tsv" \
      --occbags "$C/bags/strict_occurrence_bags.tsv" \
      --out "$W/copies/shard_$i" > "$W/shardlog/chain_$i.log" 2>&1
}

run_shard() {
  local i=$1
  local t0=$(date +%s)
  # clear every artifact of the bad attempt so a crash here cannot leave a
  # half-old half-new shard behind
  rm -f "$W/shardlog/done_$i" "$W/shardlog/FAILED_$i" "$W/copies/shard_$i.tsv" \
        "$W/hits/shard_$i.tsv" "$W/hits/shard_$i.tsv.gz" "$W/hits/shard_$i.raw"

  local sf="$TD/shard_$i.fna" want got acc path
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

  local db="$TD/db_$i"
  "$BL/makeblastdb" -in "$sf" -dbtype nucl -out "$db" -title "s$i" >/dev/null 2>&1 || {
    echo "shard $i makeblastdb failed" | tee -a "$W/shardlog/FAILED_$i"; rm -f "$sf"; return 1; }
  if ! "$BL/blastn" -task dc-megablast -query "$QL" -db "$db" \
        -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads $THR \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -out "$W/hits/shard_$i.raw" 2>"$W/shardlog/blast_$i.err"; then
    echo "shard $i blastn failed" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$W/hits/shard_$i.raw" "$sf" "$db".*; return 1
  fi

  # filter and count distinct query ids in ONE pass.  blast emits a query's
  # hits contiguously, so a change of field 1 is a new query.  This is the
  # check that exit status cannot do.
  local nq
  if ! nq=$(awk -F'\t' -v out="$W/hits/shard_$i.tsv" \
        '{if($1!=p){nq++;p=$1} if($4>=200) print >> out} END{print nq+0}' \
        "$W/hits/shard_$i.raw"); then
    echo "shard $i filter failed" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$W/hits/shard_$i.raw" "$W/hits/shard_$i.tsv" "$sf" "$db".*; return 1
  fi
  if [ "$nq" -lt "$MINQ" ]; then
    echo "shard $i SHORT QUERY READ: only $nq distinct query ids reached the output (floor $MINQ of $want_n) -- blastn exited 0 but did not read the whole query file" \
      | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$W/hits/shard_$i.raw" "$W/hits/shard_$i.tsv" "$sf" "$db".*; return 1
  fi
  echo "shard $i distinct query ids: $nq / $want_n"

  if ! chain_shard "$i"; then
    echo "shard $i chaining failed" | tee -a "$W/shardlog/FAILED_$i"
    rm -f "$sf" "$db".*; return 1
  fi
  gzip -f "$W/hits/shard_$i.tsv"; rm -f "$W/hits/shard_$i.raw" "$sf" "$db".*
  echo "$(($(date +%s)-t0))" > "$W/shardlog/done_$i"
  echo "shard $i done in $(cat "$W/shardlog/done_$i")s"
}

if [ "$T" -eq 8 ]; then
  # chain-only recovery for shard 47
  i=47
  if [ ! -s "$W/hits/shard_$i.tsv" ]; then
    echo "FATAL: shard $i filtered hits absent; it must be re-blasted, not re-chained"; exit 1
  fi
  nq=$(cut -f1 "$W/hits/shard_$i.tsv" | uniq | wc -l)
  echo "shard $i filtered hits carry $nq distinct query ids"
  if [ "$nq" -lt "$MINQ" ]; then
    echo "shard $i hits are themselves short ($nq); re-chaining would launder bad data" | tee -a "$W/shardlog/FAILED_$i"; exit 1
  fi
  t0=$(date +%s)
  rm -f "$W/shardlog/FAILED_$i" "$W/copies/shard_$i.tsv"
  if ! chain_shard "$i"; then
    echo "shard $i chaining failed again" | tee -a "$W/shardlog/FAILED_$i"; exit 1
  fi
  gzip -f "$W/hits/shard_$i.tsv"; rm -f "$W/hits/shard_$i.raw"
  echo "$(($(date +%s)-t0))" > "$W/shardlog/done_$i"
  echo "shard $i re-chained in $(cat "$W/shardlog/done_$i")s"
else
  run_shard "${REDO[$T]}" || echo "shard ${REDO[$T]} returned non-zero (recorded)"
fi
echo "task $T finished; shards done: $(ls "$W/shardlog"/done_* 2>/dev/null | wc -l)/$NSHARD; failed markers: $(ls "$W/shardlog"/FAILED_* 2>/dev/null | wc -l)"

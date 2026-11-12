#!/bin/bash
#SBATCH --job-name=isall_A3c
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
# A3c -- finish BAGS005506 (7,014 members) with threads instead of one core.
# A3 gave every bag one core so 60 bags ran at once; that was right while 6,143
# bags were outstanding and wrong once this was the only one left, with 59
# cores idle and ~10 h projected against a 12 h wall.  blastn parallelises over
# queries and there are 7,014 of them, so this is the same computation, not a
# different one -- identical task, evalue, identity, max_target_seqs and
# outfmt.  Whichever run finishes first writes the same file.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
b=BAGS005506
d=$(mktemp -d "${TMPDIR:-/tmp}/bagdb.XXXXXX")
"$BL/makeblastdb" -in "$W/bagfa/$b.fna" -dbtype nucl -out "$d/db" -title "$b" >/dev/null
"$BL/blastn" -task megablast -query "$W/bagfa/$b.fna" -db "$d/db" \
  -evalue 1e-50 -perc_identity 90 -max_target_seqs 100000 -num_threads 48 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
  -out "$W/bagaln/$b.tsv.fast"
rm -rf "$d"
rows=$(wc -l < "$W/bagaln/$b.tsv.fast")
echo "$b rows: $rows"
[ "$rows" -gt 0 ] || { echo "FATAL: empty"; exit 1; }
# publish only if A3's own slower copy has not already landed
if [ ! -s "$W/bagaln/$b.done" ]; then
  mv "$W/bagaln/$b.tsv.fast" "$W/bagaln/$b.tsv"
  grep -c '^>' "$W/bagfa/$b.fna" > "$W/bagaln/$b.done"
  echo "published $b"
else
  echo "A3 finished it first; leaving its file in place"
  rm -f "$W/bagaln/$b.tsv.fast"
fi

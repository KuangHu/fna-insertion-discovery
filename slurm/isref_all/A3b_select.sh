#!/bin/bash
#SBATCH --job-name=isall_A3b
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.out
# STAGE A3b -- representative selection, split out of A3.
#
# A3's own selection step would have died the same way A2 did.  Measured pair
# density inside a bag is 0.89-1.00 of n^2, so the 6,144 bags hold ~1e8
# aligned pairs and BAGS005506 alone ~4.9e7; a dict keyed by pair is hundreds
# of GB.  236 streams instead, one member's alignment block at a time: 18 MB
# peak on the densest bag.  A3's alignment work is kept -- it is checkpointed
# per bag and is what this consumes.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
S=/global/home/users/kh36969/fna_based_mgefinder_project/scripts
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin

nd=$(find "$W/bagaln" -name '*.done' -printf . | wc -c)
nw=$(find "$W/bagfa" -name '*.fna' -printf . | wc -c)
echo "bags aligned $nd / $nw"
[ "$nd" -eq "$nw" ] || { echo "FATAL: alignment incomplete -- not selecting on a partial set"; exit 1; }
[ -s "$W/tmp/A3_failed.txt" ] && { echo "FATAL: alignment failures present"; cat "$W/tmp/A3_failed.txt"; exit 1; }

"$PY" -I "$S/236_select_queries.py" \
  --bagaln "$W/bagaln" --bagfa "$W/bagfa" \
  --cds "$C/bags/insertion_cds_strict.tsv" \
  --occbags "$C/bags/strict_occurrence_bags.tsv" \
  --fna "$C/dedup_insertions.fna" --out "$W"

echo "queries: $(grep -c '^>' "$W/queries.fna")"

# cross-bag equivalence on the REPRESENTATIVES ONLY -- the set it is used for
"$BL/makeblastdb" -in "$W/queries.fna" -dbtype nucl -out "$W/tmp/qdb" -title q >/dev/null
"$BL/blastn" -task megablast -query "$W/queries.fna" -db "$W/tmp/qdb" \
  -evalue 1e-50 -perc_identity 90 -max_target_seqs 100000 -num_threads 60 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
  -out "$W/queries_selfaln.tsv"
echo "query self-alignment rows: $(wc -l < "$W/queries_selfaln.tsv")"

"$PY" -I "$S/237_query_equivalence.py" --aln "$W/queries_selfaln.tsv" \
  --queries "$W/queries.tsv" --out "$W/query_equivalence.tsv"

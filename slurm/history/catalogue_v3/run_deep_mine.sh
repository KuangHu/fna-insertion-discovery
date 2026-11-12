#!/bin/bash
#SBATCH --job-name=deep_mine
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=72:00:00
# EXHAUSTIVE MINING: one panel per eligible seed, k=12.
#
# The rarefaction said E. coli was still at 64% of its 100->200 marginal yield
# at 400 panels -- not saturated. So seed every genome that can anchor a k=12
# clique instead of sampling 200-400 of them.
#
# k=12 rather than 24: measured 95.9% of k=24's USABLE events per genome from
# 27% fewer genomes, so at fixed compute it reaches more of the pool.
#
# Arm A writes into the EXISTING full_<sp>/armA so already-processed genomes are
# skipped; only the newly recruited ones are computed.
#
# Old and new Layer 2 output are deduplicated TOGETHER at the end: the event key
# is (canonical insert, canonical +/-100bp context), so a panel that re-finds a
# known event merges rather than double-counting.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
for t in minimap2 minigraph gfatools samtools python3; do
  command -v "$t" >/dev/null || { echo "FATAL: $t missing"; exit 4; }
done

O=$W/deep_$SP            # new panels/armB/l2
A=$W/full_$SP/armA       # SHARED Arm A, reused across the old and deep runs
mkdir -p "$O/armB" "$O/l2" "$A"

echo "### 1. panels: one per eligible seed, k=12"
if [ ! -s "$O/panels_panels.tsv" ]; then
  python3 $P/scripts/91_seed_panel_sampler.py --edges "$EDGES" --k 12 \
      --n-seeds "$NSEEDS" --genome-dir "$GDIR" --out "$O/panels" || exit 5
fi
cat $O/panels_manifests/*.manifest | sort -u > $O/armA_all.txt
NM=$(ls $O/panels_manifests/*.manifest | wc -l)
echo "   $NM panels, $(wc -l < $O/armA_all.txt) genomes recruited"

echo "### 2. Arm A on the newly recruited genomes (existing are skipped)"
armA_one(){ b=$(basename "$1" .fna); [ -s "$2/${b}_families.tsv" ] && return 0
  python3 "$3/scripts/20_armA_multicopy.py" "$1" --sample-id "$b" \
      --out "$2/$b" --threads 8 >/dev/null 2>&1; }
export -f armA_one
xargs -P 8 -I{} bash -c 'armA_one "$@"' _ {} "$A" "$P" < $O/armA_all.txt
echo "   armA present: $(ls $A/*_families.tsv 2>/dev/null | wc -l)"

echo "### 3. Arm B"
armB_one(){ CID=$(basename "$1" .manifest); OUT=$2/armB/$CID
  [ -s "$OUT/armB_events.tsv" ] && return 0; mkdir -p "$OUT"
  python3 "$3/scripts/30_armB_graph.py" --manifest "$1" --outdir "$OUT" \
      --cluster-id "$CID" --threads 16 >/dev/null 2>&1
  rc=$?; [ $rc -eq 3 ] && echo "$CID" >> "$2/armB/zero_bubble_panels.txt"
  [ $rc -ne 0 ] && [ $rc -ne 3 ] && echo "$CID rc=$rc" >> "$2/armB/failed_panels.txt"
  return 0; }
export -f armB_one
ls $O/panels_manifests/*.manifest | xargs -P 4 -I{} bash -c 'armB_one "$@"' _ {} "$O" "$P"
NB=$(ls $O/armB/*/armB_events.tsv 2>/dev/null | wc -l)
echo "   armB $NB/$NM"

echo "### 4. Layer 2"
l2_one(){ CID=$(basename "$1" .manifest); AB=$2/armB/$CID; OUT=$2/l2/$CID
  [ -s "$OUT/loci.tsv" ] && return 0; [ -s "$AB/armB_events.tsv" ] || return 0
  mkdir -p "$OUT"
  awk -F'\t' -v OFS='\t' 'NR==1{for(j=1;j<=NF;j++)h[$j]=j; print "locus_id","seed_genome","contig","start","end"; next}
   {s=$h["junction_L"]-200; e=$h["junction_R"]+200; if(s<0)s=0;
    if($h["empty_ref_genome"]!="" && $h["empty_ref_contig"]!="")
       print $h["event_id"],$h["empty_ref_genome"],$h["empty_ref_contig"],s,e}' \
   "$AB/armB_events.tsv" > "$OUT/loci_seed.tsv"
  # never leave an empty dir behind: downstream globs on l2/*/ and an empty one
  # crashed the M join twice after dedup had already succeeded
  [ $(($(wc -l < "$OUT/loci_seed.tsv")-1)) -le 0 ] && { rm -rf "$OUT"; return 0; }
  python3 "$3/scripts/70_allele_reconstructor.py" --manifest "$1" \
      --loci "$OUT/loci_seed.tsv" --outdir "$OUT" --min-insert 500 >/dev/null 2>&1; }
export -f l2_one
ls $O/panels_manifests/*.manifest | xargs -P 4 -I{} bash -c 'l2_one "$@"' _ {} "$O" "$P"
NL=$(ls $O/l2/*/loci.tsv 2>/dev/null | wc -l)
echo "   l2 $NL panels"
[ "$NL" -lt 1 ] && { echo "FATAL: no Layer 2 output"; exit 8; }

echo "### 5. dedup OLD + DEEP together, then M join and catalogue"
M=$W/merged_$SP; mkdir -p $M
python3 $P/scripts/100_event_dedup.py --recon "$W/full_$SP/l2/*/" "$O/l2/*/" \
    --out $M/dedup || exit 9
python3 $P/scripts/99_mobility_join.py --recon "$W/full_$SP/l2/*/" "$O/l2/*/" \
    --arma "$A" --events $M/dedup_loci_to_event.tsv --out $M/mobility || exit 10
python3 $P/scripts/86_catalogue.py --recon "$W/full_$SP/l2/*/" "$O/l2/*/" \
    --events $M/dedup_loci_to_event.tsv --mobility $M/mobility_mobility.tsv \
    --species "$SP" --out $W/catalogue_deep/$SP || exit 11
echo "### $SP DEEP COMPLETE"

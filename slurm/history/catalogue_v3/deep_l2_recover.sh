#!/bin/bash
#SBATCH --job-name=deep_l2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=48:00:00
# STEP 1 -- redeem the deep Arm B that is already computed.
#
# The four exhaustive jobs died on the 72h wall INSIDE Arm B, so 3,176 panels
# have a finished armB_events.tsv and no Layer 2 at all. Arm B is the expensive
# stage; running Layer 2 on what exists costs a fraction of re-mining and is
# independent of whether the exhaustive run is ever resumed.
#
# Panels are taken one at a time -- a species does not have to finish Arm B
# before its completed panels are converted.
#
# All four species run inside ONE allocation: cf1 is EXCLUSIVE and bills 64
# cores per job regardless of the request.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
command -v minimap2 >/dev/null || { echo "FATAL: minimap2 missing"; exit 4; }

l2_one(){
  AB=$1; O=$2; P=$3; M=$4
  CID=$(basename "$AB")
  OUT=$O/l2/$CID
  [ -s "$OUT/loci.tsv" ] && return 0
  [ -s "$AB/armB_events.tsv" ] || return 0
  MAN=$M/$CID.manifest
  [ -s "$MAN" ] || return 0
  mkdir -p "$OUT"
  awk -F'\t' -v OFS='\t' 'NR==1{for(j=1;j<=NF;j++)h[$j]=j; print "locus_id","seed_genome","contig","start","end"; next}
   {s=$h["junction_L"]-200; e=$h["junction_R"]+200; if(s<0)s=0;
    if($h["empty_ref_genome"]!="" && $h["empty_ref_contig"]!="")
       print $h["event_id"],$h["empty_ref_genome"],$h["empty_ref_contig"],s,e}' \
   "$AB/armB_events.tsv" > "$OUT/loci_seed.tsv"
  # an empty dir left behind crashes the downstream globs
  [ $(($(wc -l < "$OUT/loci_seed.tsv")-1)) -le 0 ] && { rm -rf "$OUT"; return 0; }
  python3 "$P/scripts/70_allele_reconstructor.py" --manifest "$MAN" \
      --loci "$OUT/loci_seed.tsv" --outdir "$OUT" --min-insert 500 >/dev/null 2>&1
  [ -s "$OUT/loci.tsv" ] || rm -rf "$OUT"
  return 0
}
export -f l2_one

for SP in sa_k24 kpneu ecoli efaecium; do
  O=$W/deep_$SP
  N=$(ls -d $O/armB/*/ 2>/dev/null | wc -l)
  echo "### $SP : $N panels with Arm B output, $(ls $O/l2/*/loci.tsv 2>/dev/null|wc -l) already have Layer 2"
  ls -d $O/armB/*/ 2>/dev/null | xargs -P 16 -I{} bash -c 'l2_one "$@"' _ {} "$O" "$P" "$O/panels_manifests"
  echo "### $SP : Layer 2 now on $(ls $O/l2/*/loci.tsv 2>/dev/null|wc -l) panels"
done
echo "### deep Layer 2 recovery complete"

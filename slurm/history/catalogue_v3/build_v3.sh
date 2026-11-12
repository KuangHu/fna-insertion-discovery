#!/bin/bash
#SBATCH --job-name=v3
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=24:00:00
# v3 = v2 (corrected first round) + the redeemed deep panels, for every species.
#
# A species with a deep run is keyed across BOTH runs in one dedup pass, which
# is the only correct way to combine them: the event key is (canonical insert,
# canonical +/-100bp context) and run-scoped locus ids keep the two runs'
# panels apart. Never add two catalogues.
#
# The 14 species with no deep run are rebuilt too, unchanged in content, so
# every species in v3 carries the same locus-id format and is comparable.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
command -v minimap2 >/dev/null || { echo "FATAL: minimap2 missing"; exit 4; }
O=$W/catalogue_v3; mkdir -p $O $W/review/v3logs

one(){
  SP=$1; W=$2; P=$3; O=$4
  B=$W/full_$SP; V=$B/v3; mkdir -p $V
  L=$W/review/v3logs/$SP.log
  # Arm A is SHARED between the two runs; sa_k12 has none of its own and draws
  # on sa_k24's, the two being the k-arms of one panel-size experiment.
  A=$B/armA; [ "$SP" = "sa_k12" ] && A=$W/full_sa_k24/armA
  RECON=("$B/l2/"*/)
  if ls $W/deep_$SP/l2/*/loci.tsv >/dev/null 2>&1; then
      RECON+=("$W/deep_$SP/l2/"*/); fi
  { echo "recon dirs: ${#RECON[@]}"
    python3 $P/scripts/100_event_dedup.py --recon "${RECON[@]}" --out $V/dedup || exit 2
    python3 $P/scripts/99_mobility_join.py --recon "${RECON[@]}" --arma $A \
        --events $V/dedup_loci_to_event.tsv --out $V/mobility --threads 4 || exit 3
    python3 $P/scripts/86_catalogue.py --recon "${RECON[@]}" \
        --events $V/dedup_loci_to_event.tsv --mobility $V/mobility_mobility.tsv \
        --species $SP --out $O/$SP || exit 4
  } > $L 2>&1
  echo "$SP rc=$?"
}
export -f one
SPECIES=$(ls $W/catalogue_v2/*_loci.tsv | xargs -n1 basename | sed 's/_loci.tsv//')
echo "$SPECIES" | tr ' ' '\n' | xargs -P 6 -I{} bash -c 'one "$@"' _ {} "$W" "$P" "$O"

echo; echo "### v2 -> v3 per species"
printf "%-16s %9s %9s %9s\n" species v2 v3 delta
t2=0; t3=0
for f in $W/catalogue_v2/*_events.tsv; do
  sp=$(basename $f _events.tsv); g=$O/${sp}_events.tsv
  a=$(( $(wc -l < $f) - 1 )); b=0; [ -s "$g" ] && b=$(( $(wc -l < $g) - 1 ))
  t2=$((t2+a)); t3=$((t3+b))
  printf "%-16s %9d %9d %+9d\n" "$sp" "$a" "$b" "$((b-a))"
done
printf "%-16s %9d %9d %+9d\n" TOTAL $t2 $t3 $((t3-t2))

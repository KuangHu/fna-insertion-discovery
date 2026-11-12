#!/bin/bash
#SBATCH --job-name=xrun4
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=24:00:00
# Cross-run merge for the four species whose deep Arm B was redeemed.
#
# Same shape as the S. pneumoniae merge: first round alone, deep alone, and
# merged, ALL with the same fixed code and the same key, so the reportable
# number is (merged - first round) and never (old catalogue vs new catalogue)
# and never a sum of two catalogues.
#
# The deep run here is PARTIAL -- these are the panels whose Arm B finished
# before the 72h wall, not the full exhaustive design -- so the gain is a
# sample of what the design would give, not its total.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"

one(){
  SP=$1; W=$2; P=$3
  X=$W/xrun_$SP; mkdir -p $X
  L=$W/review/v2logs/xrun_$SP.log
  { echo "### $SP"
    NP=$(ls $W/deep_$SP/l2/*/loci.tsv 2>/dev/null | wc -l)
    echo "deep panels with Layer 2: $NP"
    [ "$NP" -lt 1 ] && { echo "no deep Layer 2, skipping"; exit 0; }
    python3 $P/scripts/100_event_dedup.py --recon "$W/full_$SP/l2/"*/ --out $X/r1   || exit 2
    python3 $P/scripts/100_event_dedup.py --recon "$W/deep_$SP/l2/"*/ --out $X/deep || exit 3
    python3 $P/scripts/100_event_dedup.py \
        --recon "$W/full_$SP/l2/"*/ "$W/deep_$SP/l2/"*/ --out $X/merged || exit 4
    python3 $P/scripts/126_newevent_support.py --run-a $W/full_$SP --run-b $W/deep_$SP \
        --merged $X/merged_loci_to_event.tsv || exit 5
  } > $L 2>&1
  echo "$SP rc=$?"
}
export -f one
echo sa_k24 kpneu ecoli efaecium | tr ' ' '\n' | xargs -P 4 -I{} bash -c 'one "$@"' _ {} "$W" "$P"

echo; echo "################ SUMMARY ################"
printf "%-10s %8s %8s %8s %9s %8s\n" species r1 deep merged gain pct
for SP in ecoli kpneu sa_k24 efaecium; do
  X=$W/xrun_$SP
  g(){ grep -m1 "DISTINCT EVENTS" $1 2>/dev/null | awk '{print $NF}'; }
  r=$(grep -A99 "out $X/r1" /dev/null 2>/dev/null); :
  a=$(awk '/DISTINCT EVENTS/{print $NF}' $W/review/v2logs/xrun_$SP.log 2>/dev/null | sed -n 1p)
  b=$(awk '/DISTINCT EVENTS/{print $NF}' $W/review/v2logs/xrun_$SP.log 2>/dev/null | sed -n 2p)
  c=$(awk '/DISTINCT EVENTS/{print $NF}' $W/review/v2logs/xrun_$SP.log 2>/dev/null | sed -n 3p)
  [ -z "$c" ] && continue
  printf "%-10s %8s %8s %8s %+9d %7.1f%%\n" "$SP" "$a" "$b" "$c" "$((c-a))" \
     "$(awk -v x=$c -v y=$a 'BEGIN{print 100*(x-y)/y}')"
done
echo; echo "### support source per species"
for SP in ecoli kpneu sa_k24 efaecium; do
  echo "--- $SP"; sed -n '/SUPPORT SOURCE/,/bounds ABOVE/p' $W/review/v2logs/xrun_$SP.log 2>/dev/null | head -6
done

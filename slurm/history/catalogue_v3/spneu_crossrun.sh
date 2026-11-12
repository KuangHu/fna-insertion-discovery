#!/bin/bash
#SBATCH --job-name=spneu_xrun
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=8:00:00
# STEP 1b -- the first cross-run merge, on the one species where the deep run
# finished. Also the test of the run-scoped id fix: spneumoniae shares 143 of
# its 200 first-round panel ids with the deep run, so before the fix this exact
# merge silently dropped loci.
#
# The number that answers "what did deep mining buy" is
#     merged distinct events  -  first-round distinct events
# BOTH sides computed with the SAME (fixed) code on the SAME key. The old
# catalogue_deep figure of 2,153 is not comparable: v1 code, and built from a
# colliding merge.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
X=$W/xrun_spneumoniae; mkdir -p $X

echo "##### A. first round alone, fixed code (the baseline to subtract)"
python3 $P/scripts/100_event_dedup.py --recon "$W/full_spneumoniae/l2/"*/ \
    --out $X/r1 || exit 2

echo; echo "##### B. deep alone, fixed code"
python3 $P/scripts/100_event_dedup.py --recon "$W/deep_spneumoniae/l2/"*/ \
    --out $X/deep || exit 3

echo; echo "##### C. MERGED -- both runs keyed together"
python3 $P/scripts/100_event_dedup.py \
    --recon "$W/full_spneumoniae/l2/"*/ "$W/deep_spneumoniae/l2/"*/ \
    --out $X/merged || exit 4

echo; echo "##### D. provenance of the merged events"
python3 - <<'PY'
import csv, collections, sys
csv.field_size_limit(sys.maxsize)
X = "/global/scratch/users/kh36969/fna_ins_discovery/xrun_spneumoniae"
def load(p):
    ev = collections.defaultdict(set)
    for r in csv.DictReader(open(p), delimiter="\t"):
        ev[r["event_id"]].add(r["locus_id"].split("/", 1)[0])
    return ev
r1 = load(X + "/r1_loci_to_event.tsv")
dp = load(X + "/deep_loci_to_event.tsv")
mg = load(X + "/merged_loci_to_event.tsv")
src = collections.Counter()
for e, runs in mg.items():
    has1 = any(x.startswith("full_") for x in runs)
    has2 = any(x.startswith("deep_") for x in runs)
    src["both" if (has1 and has2) else ("first_round_only" if has1 else "deep_only")] += 1
print()
print("  first round alone        %6d events" % len(r1))
print("  deep alone               %6d events" % len(dp))
print("  MERGED                   %6d events" % len(mg))
print()
print("  of the merged events:")
for k in ("both", "first_round_only", "deep_only"):
    print("    %-20s %6d" % (k, src[k]))
print()
print("  NET GAIN from deep mining = merged - first round = %+d  (%.1f%%)"
      % (len(mg) - len(r1), 100.0 * (len(mg) - len(r1)) / max(1, len(r1))))
print()
print("  naive addition would have claimed %d (+%d) -- wrong, the runs overlap"
      % (len(r1) + len(dp), len(dp)))
PY

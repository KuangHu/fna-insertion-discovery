#!/bin/bash
#SBATCH --job-name=v3bags
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=24:00:00
# Database and bag corpus on v3 (first round + redeemed deep panels, 109,709 events).
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
for t in prodigal mmseqs; do command -v $t >/dev/null || { echo "FATAL: $t"; exit 4; }; done
D=$W/database_v3; C=$W/cds_v3; B=$W/bags_v3; mkdir -p $D $C $B
python3 $P/scripts/110_build_database.py --catalogue $W/catalogue_v3 --out $D/insertions || exit 2
python3 $P/scripts/111_cds_cluster.py --inserts $D/insertions_inserts.fna \
    --db $D/insertions.tsv --out $C --threads 60 || exit 3
python3 $P/scripts/112_build_bags.py --db $D/insertions.tsv --cds $C/insert_cds.tsv \
    --inserts $D/insertions_inserts.fna --orf-gff $C/orfs.gff --out $B/insertions || exit 4
echo "### bags_v3 complete"

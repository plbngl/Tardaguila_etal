#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=multiome_test1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=job_name_%j.log
#SBATCH --mem=52G
#SBATCH --time=36:00:00


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

export PATH=/home/paola.benaglio/cellranger-arc-2.0.2:$PATH

cd /group/soranzo/paola.benaglio/k562_multiome
cellranger-arc count --id=MCO_1281 \
                      --reference=/group/soranzo/paola.benaglio/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --libraries=/group/soranzo/paola.benaglio/k562_multiome/library_1281.csv \
                      --localcores=16 \
                      --localmem=47 \
                      --jobmode=local

echo "========================"
echo "Completed: $(date)"


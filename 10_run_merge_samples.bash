#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=SeuratMerge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --output=SeuratMerge_%j.log
#SBATCH --mem=160G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


eval "$(conda shell.bash hook)"
conda activate renv_multiome


Rscript merge_samples.R  >> merge_log

conda deactivate

echo "========================"
echo "Completed: $(date)"

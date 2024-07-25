#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=SeuratSecondPass
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=SeuratSecondPass_%j.log
#SBATCH --mem=108G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


eval "$(conda shell.bash hook)"
conda activate renv_multiome

samp=$NAME

Rscript /group/soranzo/paola.benaglio/k562_multiome/scripts/seurat_secondpass_rep.R $samp


conda deactivate

echo "========================"
echo "Completed: $(date)"

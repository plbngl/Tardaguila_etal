#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=SeuratFirstPass
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --output=SeuratFirstPass_%j.log
#SBATCH --mem=108G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


eval "$(conda shell.bash hook)"
conda activate renv_multiome

samp=$NAME
rna_min_feat=500
atac_min_frag=1000
MITO_max=10

Rscript /group/soranzo/paola.benaglio/k562_multiome/scripts/seurat_firstpass.R $samp $rna_min_feat $atac_min_frag $MITO_max >> seurat_log.$samp

conda deactivate

echo "========================"
echo "Completed: $(date)"

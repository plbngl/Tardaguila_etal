#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=Cellbender
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=gpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --output=Cellbender_%j.log
#SBATCH --mem=48G
#SBATCH --gres=gpu:2



echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

eval "$(conda shell.bash hook)" 
conda activate cellbender

cd /group/soranzo/paola.benaglio/k562_multiome

samp=$NAME

mkdir -p processing_outputs/$samp


cellbender remove-background \
      --cuda \
      --input $samp/outs/raw_feature_bc_matrix.h5 \
      --output processing_outputs/$samp/cellbender_gex.h5 \
      --exclude-feature-types Peaks

conda deactivate

conda activate cellbender_0.2.2

cd processing_outputs/$samp

ptrepack --complevel 5 cellbender_gex.h5:/matrix cellbender_gex_seurat.h5:/matrix

conda deactivate

echo "========================"
echo "Completed: $(date)"


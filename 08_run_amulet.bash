#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=amulet
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=amulet_%j.log
#SBATCH --mem=64G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"



eval "$(conda shell.bash hook)"
conda activate renv_multiome

samp=$NAME
BARCODES=/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/$NAME/intermediate/keep_barcodes_step1.txt

#mkdir -p /group/soranzo/paola.benaglio/k562_multiome/processing_outputs/$samp

cd /group/soranzo/paola.benaglio/k562_multiome/processing_outputs/$samp

Rscript /group/soranzo/paola.benaglio/k562_multiome/scripts/amulet_2.R /group/soranzo/paola.benaglio/k562_multiome/$samp/outs/atac_fragments.tsv.gz $BARCODES >> amulet_log.$samp

conda deactivate

echo "========================"
echo "Completed: $(date)"

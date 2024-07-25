#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=snATAC
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=snATAC_%j.log
#SBATCH --mem=160G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


module load samtools
module load bedtools
module load R/4.3.1

eval "$(conda shell.bash hook)"
conda activate python_snatac

BAM=/group/soranzo/paola.benaglio/k562_multiome/$NAME/outs/atac_possorted_bam.bam
OUTDIR=/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/$NAME/snATAC_matrices
BARCODES=/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/$NAME/intermediate/keep_barcodes_step1.txt


mkdir -p $OUTDIR
cd $OUTDIR

python3 /group/soranzo/paola.benaglio/k562_multiome/scripts/snATAC_10x_matrices_pipeline.py -m 6 -b $BAM -o $OUTDIR -n $NAME --keep-bc $BARCODES 
#--skip-convert




echo "========================"
echo "Completed: $(date)"

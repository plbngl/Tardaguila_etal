#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=merge_frag
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=merge_frag_%j.log
#SBATCH --mem=96G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


module load htslib-tools
cd /group/soranzo/paola.benaglio/k562_multiome/

for SAMPLE in MCO_1278 MCO_1279 MCO_1280 MCO_1281; do zcat ${SAMPLE}/outs/atac_fragments.tsv.gz | awk -v SAMPLE=$SAMPLE \ 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,SAMPLE"_"$4,$5}\'; done | sort -k1,1 -k2,2n -S 72G | bgzip -c -@ 16 > processing_outputs/merged.atac_fragments.tsv.gz 
tabix -p bed processing_outputs/merged.atac_fragments.tsv.gz



echo "========================"
echo "Completed: $(date)"


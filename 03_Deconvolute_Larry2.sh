#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=larryCounts2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paola.benaglio@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=larryCounts2_%j.log
#SBATCH --mem=64G



echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

module load samtools
module load bwa-mem2


newref='/group/soranzo/paola.benaglio/references/modified_site/GFP_transgenev4.fa'
mkdir -p /group/soranzo/paola.benaglio/k562_multiome/deconvolute_LARRY/
cd /group/soranzo/paola.benaglio/k562_multiome/deconvolute_LARRY/

for samp in  MCO_1278 MCO_1279 MCO_1280 MCO_1281;

do

echo Process sample $samp

dir=/group/soranzo/paola.benaglio/k562_multiome/$samp/outs

samtools view -b -f 4 $dir/gex_possorted_bam.bam -@ 32 > $samp.unmapped.bam
samtools sort -@ 32 -n -o $samp.unmapped.sorted.bam $samp.unmapped.bam


samtools fastq -@ 32 $samp.unmapped.sorted.bam -T CB,UB,xf > $samp.unmapped.fq
sed -i 's/\t/_/g' $samp.unmapped.fq


### Remap to modified reference with deletions 
bwa-mem2 mem -M -t 16 $newref $samp.unmapped.fq | samtools sort -o $samp.realigned.bam - 
samtools index $samp.realigned.bam

### Remove multimappers to keep only reads uniquely mapping to each GFP barcode

samtools view -q 5 -b $samp.realigned.bam > $samp.realigned.unique.bam 
samtools index $samp.realigned.unique.bam

while read -r geno;
do 
samtools view $samp.realigned.unique.bam $geno \
 | cut -f1 | sed 's/_/ /g' | awk '{print $3}' | sort | uniq > $samp.$geno.unique.bcs
done < ../list_gfp


while read -r geno;
do 
samtools view $samp.realigned.unique.bam $geno \
 | cut -f1 | sed 's/_/ /g' | awk -v var="$geno" '{print $3,$4,var}' | sort | uniq > $samp.$geno.unique_bcs_umi
done < ../list_gfp


while read -r geno;
do 
samtools view $samp.realigned.unique.bam $geno \
 | cut -f1 | sed 's/_/ /g' | awk -v var="$geno" '{print $2,$3,$4,var}' | sort  > $geno.$samp.bcs_umi
done < ../list_gfp

cat *.$samp.bcs_umi | sort  > $samp.all_geno_bc_umi


done
echo "========================"
echo "Completed: $(date)"


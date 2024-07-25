#!/bin/bash

dir='/group/soranzo/paola.benaglio/k562_multiome'
newref='/group/soranzo/paola.benaglio/references/modified_site/chr7_mod.fa'

module load samtools   ###  /ssu/gassu/software/modulefiles/samtools (v1.14)
module load bwa-mem2   ###  /ssu/gassu/software/modulefiles/bwa-mem2 (v2.2.1) 
module load R/4.3.1    ###  /share/apps/modulefiles

mkdir -p $dir/deconvolute_ATAC

for samp in MCO_1278 MCO_1279 MCO_1280 MCO_1281;

do

echo "~~~~ Process sample $samp ~~~~~~~~~"
### Extract reads mapping to region of interest
samtools view  -b $dir/$samp/outs/atac_possorted_bam.bam chr7:101855762-101857539 -o $dir/deconvolute_ATAC/$samp.region.bam
samtools index $dir/deconvolute_ATAC/$samp.region.bam 
samtools sort -n -o $dir/deconvolute_ATAC/$samp.region.qsort.bam $dir/deconvolute_ATAC/$samp.region.bam

samtools fastq -@ 8 $dir/deconvolute_ATAC/$samp.region.qsort.bam \
    -1 $dir/deconvolute_ATAC/$samp.region.end1.fq \
    -2 $dir/deconvolute_ATAC/$samp.region.end2.fq \
    -T CB \
    -0 /dev/null -s /dev/null -n

### append barcode to read name
sed -i 's/\t/_/' $dir/deconvolute_ATAC/$samp.region.end1.fq
sed -i 's/\t/_/' $dir/deconvolute_ATAC/$samp.region.end2.fq

### Remap to modified reference with deletions 
bwa-mem2 mem -M -t 8 $newref $dir/deconvolute_ATAC/$samp.region.end1.fq $dir/deconvolute_ATAC/$samp.region.end2.fq | samtools sort -o $dir/deconvolute_ATAC/$samp.region.realigned.bam - 
samtools index $dir/deconvolute_ATAC/$samp.region.realigned.bam

### Remove multimappers to keep only reads uniquely mapping to each genotype (wt, del16bp or del80bp)

samtools view -q 1 -b $dir/deconvolute_ATAC/$samp.region.realigned.bam > $dir/deconvolute_ATAC/$samp.region.unique.bam 
samtools index $dir/deconvolute_ATAC/$samp.region.unique.bam 



### Extract read names / barcodes
samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam Unmodified | cut -f1 | sort | uniq | sed 's/_/ /' | awk '{print $2}' > $dir/deconvolute_ATAC/$samp.reads_unmod
samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam 16bp_del | cut -f1 | sort | uniq | sed 's/_/ /' | awk '{print $2}' > $dir/deconvolute_ATAC/$samp.reads_16bp_del
samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam 80bp_del | cut -f1 | sort | uniq | sed 's/_/ /' | awk '{print $2}' > $dir/deconvolute_ATAC/$samp.reads_80bp_del

samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam Unmodified | cut -f1 | sed 's/_/ /' | awk '{print $2}' | sort | uniq > $dir/deconvolute_ATAC/$samp.bcs_unmod
samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam 16bp_del | cut -f1 | sed 's/_/ /' | awk '{print $2}' | sort | uniq > $dir/deconvolute_ATAC/$samp.bcs_16bp_del
samtools view $dir/deconvolute_ATAC/$samp.region.unique.bam 80bp_del | cut -f1 | sed 's/_/ /' | awk '{print $2}' | sort | uniq > $dir/deconvolute_ATAC/$samp.bcs_80bp_del

#### Genotype the base edit
samtools mpileup -f $newref  -A -q 1 -B -d 0 -Q 0 --output-QNAME --no-output-ins --no-output-del --no-output-ends \
-r Unmodified $dir/deconvolute_ATAC/$samp.region.unique.bam > $dir/deconvolute_ATAC/$samp.mpileup.txt


#### Genotype the base edit using original cellranger bam

samtools mpileup -f /group/soranzo/paola.benaglio/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
  -q 1 -B -d 0 -Q 0 --output-extra QNAME,CB --no-output-ins --no-output-del --no-output-ends \
-r chr7:101856500-101856800 $dir/$samp/outs/atac_possorted_bam.bam > $dir/deconvolute_ATAC/$samp.mpileup_300bp.txt

Rscript $dir/scripts/Barcode_counts.R  $samp $dir/deconvolute_ATAC

done

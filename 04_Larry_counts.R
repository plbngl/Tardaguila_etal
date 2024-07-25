.libPaths(new ="/home/paola.benaglio/R/x86_64-pc-linux-gnu-library/4.3")
library(pheatmap)
library(stringr)
library(tidyr)
library(RColorBrewer)

setwd("/group/soranzo/paola.benaglio/k562_multiome/deconvolute_LARRY/")




for (samp in c('MCO_1278','MCO_1279','MCO_1280','MCO_1281')) {

bcs = read.table(paste0(samp, ".all_geno_bc_umi"), sep="\ ")
bcs = bcs[bcs$V2!="",-1] ## xf:i:0 for all

colnames(bcs)= c("CellBC", "MolBC", "GFPbc")
bcs$UMI_GFP = paste(bcs$MolBC, bcs$GFPbc, sep="-")
bcs$Cell_GFP = paste(bcs$CellBC, bcs$GFPbc, sep="-")
bcs$UMI_GFP_Cell = paste(bcs$MolBC, bcs$GFPbc,bcs$CellBC, sep="-")
length(unique(bcs$UMI_GFP_Cell))   ## OUTPUT

#number of reads supporting a given UMI/GFP-BC pair in a particular cell
n_reads = table(bcs$UMI_GFP_Cell)  
summary(c(n_reads)) ## OUTPUT
bcs = bcs[!duplicated(bcs$UMI_GFP_Cell),]
#number of number of UMIs supporting a given cell/GFP-BC pair
n_umis = aggregate(MolBC~Cell_GFP,bcs, length)
colnames(n_umis)[2] = "No_UMIs"
summary(n_umis$No_UMIs)  ## OUTPUT
filt2 = n_umis[n_umis$No_UMIs>=3,] ### filter for min 3 
bcs_filt = bcs[bcs$Cell_GFP %in% filt2$Cell_GFP, ]
bcs_filt_un = bcs_filt[!duplicated(bcs_filt$Cell_GFP),]
##### There should be only one GFP (GT) per cell:
ag = aggregate(GFPbc ~ CellBC,bcs_filt_un, length )
colnames(ag)[2] = "No_GFPbcs"
uniGT = ag$CellBC[ag$No_GFPbcs==1]



######## Plot Barcode Statistics ######## 
######################################### 
pdf(paste0(samp, ".barcode_statistics.pdf"), width = 11,height = 8)
par(mfrow=c(2,3))

bp2=table(bcs_filt_un$GFPbc)
cb2=spread(data.frame(bp2, gts = c('16bp_del','80bp_del','80bp_del','80bp_del',
                                   "HET", 'MUT','MUT','MUT', 'WT','WT','WT')), Var1, Freq, fill=0)
rownames(cb2)=cb2$gts
cb2 = cb2[,-1]

bar1 = barplot(bp2, col =rainbow(11), names.arg = NA, legend=rownames(bp2), main = "#Cell BCs per experiment", 
        ylim=c(0, max(bp2)+(max(bp2)*0.2)))
text(bar1, bp2,bp2,  pos=3)
y=rowSums(cb2)
bar2 = barplot(t(cb2), col =rainbow(11), 
               #legend=colnames(cb2), 
               main = "#Cell BCs per genotype",
               ylim=c(0, max(y)+(max(y)*0.2)))
text(bar2, y,y,  pos=3)


bp=table(ag$No_GFPbcs) 
bb = barplot(bp,ylim=c(0, max(bp)+(max(bp)*0.2)), main="Cell barcodes with 1+ experiment BC", xlab="# GFP BCs ", ylab="# Cell BCs")
text(bb, bp,bp, pos=3)

mm = merge(bcs_filt_un, ag, by=1)

tb = table(mm$GFPbc, mm$No_GFPbcs)
barplot(tb, col =rainbow(11), beside = T, legend=rownames(tb), main = "Cell BCs with 1+ experiment BC",
        xlab="# GFP BCs ", ylab="# Cell BCs")

#### look at how many umi in each category
mm2 = merge(mm,n_umis ,by='Cell_GFP')
boxplot(No_UMIs~No_GFPbcs, mm2, main = "#UMIs per cell with 1+ assigned GFP bcs")

mms = subset(mm2,No_GFPbcs==2 )
mms = mms[order(mms$No_UMIs),]
agg = aggregate(No_UMIs~CellBC, mms, paste, collapse=",")
vect = str_split_fixed(agg[,2], "\\,",2)
hist(as.numeric(vect[,2])/as.numeric(vect[,1]), breaks=20, xaxt='n', xlab="UMIs Cell1 / UMI Cell2", 
     main = "UMI ratio for Cell BC assigned to 2 experiments")
axis(side=1, at=1:20)

dev.off()

##############################################
##############################################

write.csv(mm[,c("CellBC", "GFPbc" ,"Cell_GFP",'No_GFPbcs')], paste0(samp,".larry_barcodes_assignemnts.csv"))

####### Compare with ATAC genotypes ##########
##############################################

gts = c('80bp_del','16bp_del','A_gt','G_gt','AG_gt')

atacbc=data.frame()
for (gg in gts){
  
  bc = readLines(paste0('/group/soranzo/paola.benaglio/k562_multiome/deconvolute_ATAC/',samp ,'.bcs_',gg))
  bc = data.frame(bc = bc, gt = gg)
  atacbc = rbind(atacbc,bc)
}

both = merge(atacbc, bcs_filt_un, by=1)
tb =table(both$GFPbc, both$gt)
cols = brewer.pal(9, "Blues")
pheatmap(tb, cluster_rows = F,cluster_cols = F, 
         file =paste0( samp, ".compare_withATAC.pdf"),
         cellwidth = 25, cellheight = 20,
         display_numbers = tb, col=cols ,fontsize_number = 12 , number_color = "gray60")


bcs_filt_unique = bcs_filt_un[bcs_filt_un$CellBC %in% uniGT,]

both = merge(atacbc, bcs_filt_unique, by=1)
tb =table(both$GFPbc, both$gt)
cols = brewer.pal(9, "Blues")
pheatmap(tb, cluster_rows = F,cluster_cols = F, 
         file =paste0( samp, ".unique.compare_withATAC.pdf"),
         cellwidth = 25, cellheight = 20,
         display_numbers = tb, col=cols ,fontsize_number = 12 , number_color = "gray60")


gts = setNames(nm = sort(unique (bcs_filt_unique$GFPbc)), 
               object = c('16bp_del','80bp_del','80bp_del','80bp_del',
                          "HET", 'MUT','MUT','MUT', 'WT','WT','WT'))
bcs_filt_unique$GT = gts[bcs_filt_unique$GFPbc]

df = cbind(ATAC=table(atacbc$gt), GFP=table(bcs_filt_unique$GT), Shared=0, Union=0)
df = rbind(df, c(0,0, nrow(both), length(unique(c(atacbc$bc, bcs_filt_unique$CellBC)))))
rownames(df)[6]="all"
df = df[,c(4,3,2,1)]

pdf(paste0(samp, ".assigned_bcs_summary.pdf"), height = 6, width = 6)
pp = barplot(df, col=rainbow(6), legend=rownames(df),
             ylim=c(0, max(df)+(max(df)*0.2)), main ="Assigned cell BCs unique", las=1)
text(pp,colSums(df), colSums((df)), pos=3)
dev.off()



####### Compare with ATAC- v2 - more stringent ############
####### ####### ####### ####### ####### ####### ####### 

atacbc2 = read.csv(paste0('/group/soranzo/paola.benaglio/k562_multiome/deconvolute_ATAC_v2/',
                          samp , '.genotype_assignment.csv'))
cols = brewer.pal(9, "Oranges")
both = merge(atacbc2, bcs_filt_un, by=1)
tb =table(both$GFPbc, both$geno)
pheatmap(tb, cluster_rows = F,cluster_cols = F, 
         file =paste0( samp, ".compare_withATAC_v2.pdf"),
         cellwidth = 25, cellheight = 20,
         display_numbers = tb, col=cols ,fontsize_number = 12 , number_color = "gray60")

both = merge(atacbc2, bcs_filt_unique, by=1)
tb =table(both$GFPbc, both$geno)
pheatmap(tb, cluster_rows = F,cluster_cols = F, 
         file =paste0( samp, ".unique.compare_withATAC_v2.pdf"),
         cellwidth = 25, cellheight = 20,
         display_numbers = tb, col=cols ,fontsize_number = 12 , number_color = "gray60")

df = cbind(ATAC=table(atacbc2$geno), GFP=table(bcs_filt_unique$GT), Shared=0, Union=0)
df = rbind(df, c(0,0, nrow(both), length(unique(c(atacbc2$geno, bcs_filt_unique$CellBC)))))
rownames(df)[6]="all"
df = df[,c(4,3,2,1)]

pdf(paste0(samp, ".assigned_bcs_summary_v2.pdf"), height = 6, width = 6)
pp = barplot(df, col=rainbow(6), legend=rownames(df),
             ylim=c(0, max(df)+(max(df)*0.2)), main ="Assigned cell BCs unique", las=1)
text(pp,colSums(df), colSums((df)), pos=3)
dev.off()
}

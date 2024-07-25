.libPaths(new = c("/home/paola.benaglio/conda_envs/renv_multiome/lib/R/library",
          "/group/soranzo/paola.benaglio/r_libraries"))
### requested memory: 72 G
Sys.setenv(RETICULATE_PYTHON="/home/paola.benaglio/conda_envs/renv_multiome/bin/python")
library(reticulate)
reticulate::use_python("/home/paola.benaglio/conda_envs/renv_multiome/bin/python")
reticulate::use_condaenv("/home/paola.benaglio/conda_envs/renv_multiome")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg') 
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat)) 
suppressMessages(library(Signac)) 
suppressMessages(library(EnsDb.Hsapiens.v86)) 
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
warnLevel <- getOption('warn')
options(warn = -1)
library(ggplot2)
suppressMessages(library(scDblFinder))
suppressMessages(library(glmGamPoi))


args = commandArgs(trailingOnly=TRUE)
sample_name        =  args[1]
sample_dir         =  sprintf("/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/%s/", sample_name)
output_dir         =  sprintf("/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/%s/pre_merge/", sample_name)
crange_dir         =  sprintf("/group/soranzo/paola.benaglio/k562_multiome/%s/outs/", sample_name)

if (FALSE) {
dir.create(output_dir)

### New Seurat object
#### load previous filtered seurat object
adata <- readRDS(file = file.path(sample_dir,'intermediate','preliminary_filtered.rds'))

#### cellbender counts -- filter them for previous bcs -- uses it as the main RNA modality
cb = Read10X_h5(file.path(sample_dir, 'cellbender_gex_seurat.h5'))
cb_counts   <- cb$'Gene Expression'
cb_counts   <- cb_counts[,colnames(adata)]

print(dim(cb_counts))

adata2 = CreateSeuratObject(counts = cb_counts)
adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')

#add in previous raw RNA data as another assay (RNA_raw) for comparison
DefaultAssay(adata) <- 'RNA'
raw_rna <-  GetAssayData(object = adata, slot = "counts")
raw_rna_assay <- CreateAssayObject(counts = raw_rna)
adata2[['RNA_raw']] <- raw_rna_assay

### Add ATAC modality using previously generated 5kb windows matrix
atac_lfmtx = read.table(paste0(sample_dir, 'snATAC_matrices/', sample_name,'.long_fmt_mtx.txt.gz'))

FALSE %in% (colnames(adata2) %in%  atac_lfmtx$V1)
FALSE %in% ( atac_lfmtx$V1 %in% colnames(adata2))

atac_lfmtx$V1 <- factor(atac_lfmtx$V1, levels=colnames(adata2))
reordered_lfm <- atac_lfmtx[order(atac_lfmtx$V1),]

atac_sm <- with(reordered_lfm,
                sparseMatrix(i=as.numeric(as.factor(V2)), j=as.numeric(V1), 
                             x=V3, dimnames=list(levels(as.factor(V2)), levels(V1))))


#create the new chromatin assay object and add to Seurat object
atac_sm       <- atac_sm[,colnames(adata2)]
grange.counts <- StringToGRanges(rownames(atac_sm), sep = c(':', '-'))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_sm       <- atac_sm[as.vector(grange.use), ]
suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations)  <- 'UCSC'
genome(annotations)          <- 'hg38'

frag.file <- file.path(crange_dir, 'atac_fragments.tsv.gz')
suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_sm, sep=c(':', '-'), 
                                                     genome='hg38', fragments=frag.file, 
                                                     min.cells=-1, min.features=-1, 
                                                     annotation=annotations))
adata2[['ATAC']] <- chrom_assay

invisible(gc())

### Compute new  metrics for cellbender rna and windows atac
qc <- read.table(file.path(crange_dir, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
qc <- as.data.frame(qc)
rownames(qc) <- qc$gex_barcode
qc <- qc[Cells(adata2), 6:length(colnames(qc))]
adata2 <- AddMetaData(adata2, qc)

DefaultAssay(adata2) <- 'ATAC'
adata2 <- TSSEnrichment(adata2)


#add previous metadata
old_meta = adata@meta.data

colkeep = c('scDblFinder.class','scDblFinder.score','scDblFinder.weighted','scDblFinder.cxds_score',
'scDblFinder.class_atac','scDblFinder.score_atac','scDblFinder.weighted_atac','scDblFinder.cxds_score_atac',
'No_assigned_GFPbc','Assigned_GFPbc','Assigned_GFPgenotype','DBL_comb')

adata2@meta.data = cbind(adata2@meta.data,old_meta[,colkeep])


amures           = read.table(file.path(sample_dir, "Amulet_selected_bc.tsv"))
colnames(amures) = paste0("amulet_",colnames(amures))
amures           = amures[rownames(adata2@meta.data), ]
adata2@meta.data = cbind (adata2@meta.data, amures)
adata2@meta.data$doublet_amulet = adata2@meta.data$amulet_q.value <0.05

saveRDS(adata2, file = file.path(output_dir,'pre_merged.rds'))
}



adata2 = readRDS(file.path(output_dir,'pre_merged.rds'))
adata2 <- adata2[, unname(which( colSums(GetAssayData(adata2, slot = "counts", assay = "RNA"))!=0))]

### Analyze and cluster
# RNA analysis
DefaultAssay(adata2) <- 'RNA'
adata2 <- SCTransform(adata2, verbose = FALSE) 
adata2 <- RunPCA(adata2) 
adata2 <- RunUMAP(adata2, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(adata2) <- 'ATAC'
adata2 <- RunTFIDF(adata2)
adata2 <- FindTopFeatures(adata2, min.cutoff='q0')
adata2 <- RunSVD(adata2)
adata2 <- RunUMAP(adata2, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')

# Multimodal analysis
adata2 <- FindMultiModalNeighbors(adata2, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
adata2 <- RunUMAP(adata2, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = .2, verbose=FALSE)


################## Plots ####################################
#############################################################

options(repr.plot.width = 12, repr.plot.height = 5)
png(file.path(output_dir,'Premerge_UMAPs_clusters.png'), width =1000, height = 350)
p1 <- DimPlot(adata2, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(adata2, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(adata2, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(file.path(output_dir,'Premerge_doublets.png'), width =1000, height = 300)
p1<- DimPlot(adata2, group.by = c("doublet_amulet"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "blue"))
p2<- DimPlot(adata2, group.by = c("scDblFinder.class_atac"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "red4"))
p3<- DimPlot(adata2, group.by = c("scDblFinder.class"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "green4"))


p1 + p2 + p3
dev.off()

options(repr.plot.width = 14, repr.plot.height = 8)
png(file.path(output_dir,'Premerge_UMAPs_UMAPs_qc.png'), width =1000, height = 600)
DefaultAssay(adata2) <- 'SCT'
p5 <- FeaturePlot(adata2, features = c("CUX1",'FYB1', 'PECAM1', 'ITGA2B'),
                  reduction = 'umap.wnn', 
                  cols = c("lightgrey","darkgreen"), ncol = 4)
p6 <- FeaturePlot(adata2, features = c("nCount_SCT", "nCount_RNA", "nCount_ATAC",'TSS.enrichment'), ncol = 4,
                  reduction = 'umap.wnn')
p7 <- FeaturePlot(adata2, features = c("nFeature_SCT", "nFeature_RNA", "nFeature_ATAC",'percent.mt'), ncol = 4,
                  reduction = 'umap.wnn')
p6 / p7 / p5
dev.off()



deconv_dir = '/group/soranzo/paola.benaglio/k562_multiome/deconvolute_LARRY/'

#### add deconvolution results
################################

larry     = read.csv(paste0(deconv_dir, sample_name, '.larry_barcodes_assignemnts.csv'), row.names=1)
larry_uni = larry[larry$No_GFPbcs==1,]
gfp1      = setNames(larry$No_GFPbcs, stringr::str_split_fixed(larry$CellBC, "\\:",3)[,3])
adata2@meta.data$No_assigned_GFPbc = gfp1[rownames(adata2@meta.data)]

gfp2      = setNames(larry_uni$GFPbc, stringr::str_split_fixed(larry_uni$CellBC, "\\:",3)[,3])
adata2@meta.data$Assigned_GFPbc = gfp2[rownames(adata2@meta.data)]

gts = setNames(nm = sort(unique (larry_uni$GFPbc)), 
               object = c('16bp_del','80bp_del','80bp_del','80bp_del',
                          "HET", 'MUT','MUT','MUT', 'WT','WT','WT'))

adata2@meta.data$Assigned_GFPgenotype = gts[adata2@meta.data$Assigned_GFPbc]

adata_assigned = subset(adata2, subset = No_assigned_GFPbc ==1 )



options(repr.plot.width = 12, repr.plot.height = 6)
png(file.path(output_dir,'Premerge_UMAPs_unique_GFP_bc.png'), width =1000, height = 600)
p1<- DimPlot(adata_assigned, reduction = "umap.rna", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("RNA") 
p2<- DimPlot(adata_assigned, reduction = "umap.atac", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPgenotype", label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 +p3 +p4 
dev.off()

##########################################################################################

saveRDS(adata2, file = file.path(output_dir,'pre_merged_clustered.rds'))

##################################################################################
#### perform higher resolution clustering and plot them in relations to QC metrics



adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = 1, verbose=FALSE)

p1 <- VlnPlot(adata2, features='nCount_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
p2 <- VlnPlot(adata2, features='nFeature_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
p3 <- VlnPlot(adata2, features='nCount_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_ATAC), linetype='dashed')
p4 <- VlnPlot(adata2, features='nFeature_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_ATAC), linetype='dashed')
p5 <- VlnPlot(adata2, features='nCount_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
p6 <- VlnPlot(adata2, features='nFeature_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')


#options(repr.plot.width=24, repr.plot.height=8)
#ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2)


p9 <- VlnPlot(adata2, features='TSS.enrichment', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$TSS.enrichment), linetype='dashed')
p10 <- VlnPlot(adata2, features='percent.mt', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$percent.mt), linetype='dashed')
p11 <- VlnPlot(adata2, features='amulet_nFrags', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$amulet_nFrags), linetype='dashed')
p12 <- VlnPlot(adata2, features='scDblFinder.score', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score), linetype='dashed')
p13 <- VlnPlot(adata2, features='scDblFinder.score_atac', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
p14 <- VlnPlot(adata2, features='amulet_q.value', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')


#options(repr.plot.width=24, repr.plot.height=8)
#ggarrange(p9, p10, p11, p12, p13, p14, ncol = 3, nrow = 2)

print(
png(file.path(output_dir,'Violin_plots_QC_byCluster.png'), width =1200, height = 1200))
ggarrange(p1& NoLegend() , p3& NoLegend(), p5& NoLegend(), p2& NoLegend(), p4& NoLegend(), p6& NoLegend(),
          p9& NoLegend(), p10& NoLegend(), p11& NoLegend(), p12& NoLegend(), p13& NoLegend(), p14 & NoLegend(), ncol = 3, nrow = 4)
dev.off()

    print(
png(file.path(output_dir,'Premerge_UMAPs_clusters_res1.png'), width =1000, height = 350))
p1 <- DimPlot(adata2, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(adata2, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(adata2, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Any sample-specific clusters?

p1 <- VlnPlot(adata2, features='nCount_SCT', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
p2 <- VlnPlot(adata2, features='nFeature_SCT', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
p3 <- VlnPlot(adata2, features='nCount_ATAC', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_ATAC), linetype='dashed')
p4 <- VlnPlot(adata2, features='nFeature_ATAC', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_ATAC), linetype='dashed')
p5 <- VlnPlot(adata2, features='nCount_RNA', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
p6 <- VlnPlot(adata2, features='nFeature_RNA', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')


#options(repr.plot.width=24, repr.plot.height=8)
#ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2)

p9 <- VlnPlot(adata2, features='TSS.enrichment', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$TSS.enrichment), linetype='dashed')
p10 <- VlnPlot(adata2, features='percent.mt', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$percent.mt), linetype='dashed')
p11 <- VlnPlot(adata2, features='amulet_nFrags', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$amulet_nFrags), linetype='dashed')
p12 <- VlnPlot(adata2, features='scDblFinder.score', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score), linetype='dashed')
p13 <- VlnPlot(adata2, features='scDblFinder.score_atac', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
p14 <- VlnPlot(adata2, features='amulet_q.value', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')


#options(repr.plot.width=24, repr.plot.height=8)
#ggarrange(p9, p10, p11, p12, p13, p14, ncol = 3, nrow = 2)


qc_clus = data.frame(barcodes = rownames(adata2@meta.data), clusters = adata2@meta.data[,"seurat_clusters"])
write.table(qc_clus, file.path(output_dir,'qc_clusters.tsv'), quote=F, row.names=FALSE)
print(png(file.path(output_dir,'Violin_plots_QC_bysample.png'), width =1200, height = 1200))
ggarrange(p1& NoLegend() , p3& NoLegend(), p5& NoLegend(), p2& NoLegend(), p4& NoLegend(), p6& NoLegend(),
          p9& NoLegend(), p10& NoLegend(), p11& NoLegend(), p12& NoLegend(), p13& NoLegend(), p14 & NoLegend(), ncol = 3, nrow = 4)
dev.off()
    
    

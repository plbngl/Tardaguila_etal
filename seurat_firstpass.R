.libPaths(new = c("/home/paola.benaglio/conda_envs/renv_multiome/lib/R/library",
          "/group/soranzo/paola.benaglio/r_libraries"))

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
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
warnLevel <- getOption('warn')
options(warn = -1)


args = commandArgs(trailingOnly=TRUE)
sample_name        =  args[1]
sample_dir         =  sprintf("/group/soranzo/paola.benaglio/k562_multiome/%s/outs/", sample_name)
output_dir         =  sprintf("/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/%s/intermediate/", sample_name)
deconv_dir         = '/group/soranzo/paola.benaglio/k562_multiome/deconvolute_LARRY/'
rna_min_features   = args[2]
atac_min_fragments = args[3]
MITO_max           = args[4]



dir.create(output_dir)
sample_color = as.numeric(substring(sample_name, 7))
ncells = read.csv(file.path(sample_dir, "summary.csv"))$Estimated.number.of.cells

# Barcode rank plots based on cellranger cell estimate

png(file.path(output_dir,'Barcodes_rank_plots.png'), width = 1000, height = 500)
par(mfrow=c(1,2), mar=c(6,6,6,6))
bc_metrics <- read.csv(file.path(sample_dir, 'per_barcode_metrics.csv'),  , stringsAsFactors=1)
plot( sort(bc_metrics$gex_umis_count,decreasing = T), log='xy', type='l', main=paste(sample_name,"RNA"),
     xlab="Barcode Rank", ylab="Number RNA UMIs", col="gray" ,lwd=3)
lines( sort(bc_metrics$gex_umis_count,decreasing = T)[1:ncells],  col="blue4" ,lwd=3)
grid()
plot(sort(bc_metrics$atac_fragments, decreasing = T), log='xy', type='l', main=paste(sample_name,"ATAC"),
     xlab="Barcode Rank", ylab="Number ATAC fragment", col='gray', lwd=3)
lines( sort(bc_metrics$atac_fragments,decreasing = T)[1:ncells],  col="dodgerblue" ,lwd=3)
grid()
dev.off()
invisible(gc())



stored_filters = c()

## Read raw data
inputdata.10x         <- Read10X_h5(file.path(sample_dir, 'raw_feature_bc_matrix.h5'))
rna_counts            <- inputdata.10x$'Gene Expression'
atac_counts           <- inputdata.10x$'Peaks'
adata                 <- CreateSeuratObject(counts=rna_counts)
adata[['percent.mt']] <- PercentageFeatureSet(adata, pattern = '^MT-')

stored_filters['total_bc'] = length(colnames(adata[["RNA"]]))
adata_sub <- subset(x = adata,subset = nFeature_RNA >= as.numeric(rna_min_features))
stored_filters['after_rna_minfeat'] = length(colnames(adata_sub[["RNA"]]))

#add cellranger QC
qc <- read.table(file.path(sample_dir, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
qc           <- as.data.frame(qc)
rownames(qc) <- qc$gex_barcode
qc           <- qc[Cells(adata_sub), ]
adata_sub    <- AddMetaData(adata_sub, qc)

### Add in ATAC data for the barcodes that passed RNA filters

atac_counts   <- atac_counts[,colnames(adata_sub)]
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts   <- atac_counts[as.vector(grange.use), ]
suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations)  <- 'UCSC'
genome(annotations)          <- 'hg38'

frag.file <- file.path(sample_dir, 'atac_fragments.tsv.gz')
suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_counts, sep=c(':', '-'), 
                                                     genome='hg38', fragments=frag.file, 
                                                     min.cells=-1, min.features=-1, 
                                                     annotation=annotations))
adata_sub[['ATAC']] <- chrom_assay
invisible(gc())

### Fist doublet detection RNA with scDblFinder to be performed on lightly filtered data

sce <- scDblFinder(GetAssayData(object = adata_sub, slot = "counts"))
sce_results = data.frame(SummarizedExperiment::colData(sce))
adata_sub@meta.data = cbind(adata_sub@meta.data,sce_results)


#filter by min ATAC number fragments 
adata_sub_atac <- subset( x = adata_sub, subset = atac_fragments >= as.numeric(atac_min_fragments))
stored_filters['after_atac_minfrag'] = length(colnames(adata_sub_atac[["RNA"]]))

### Fist doublet detection ATAC with scDblFinder to be performed on lightly filtered data

DefaultAssay(adata_sub_atac) <- 'ATAC'
sce_atac <- scDblFinder(GetAssayData(object = adata_sub_atac, slot = "counts"), 
                        artificialDoublets=1, aggregateFeatures=TRUE, 
                        nfeatures=25, processing="normFeatures")


sce_results_atac = data.frame(SummarizedExperiment::colData(sce_atac))
colnames(sce_results_atac) = paste(colnames(sce_results_atac), "atac", sep="_")
adata_sub_atac@meta.data = cbind(adata_sub_atac@meta.data,sce_results_atac)


#remove multiplets and other cellranger excluded cells 
adata_sub_multiplet <- subset(
  x = adata_sub_atac,
  subset = excluded_reason != 1 
)
stored_filters['after_cr_multiplets'] = length(colnames(adata_sub_atac[["RNA"]]))



#remove high mito content cells
adata_sub_mito = subset(adata_sub_multiplet, subset = percent.mt < as.numeric(MITO_max))
stored_filters['after_mito'] = length(colnames(adata_sub_mito[["RNA"]]))


### Compute intermediate bulk metrics of filtered object
adata <- adata_sub_mito
DefaultAssay(adata) <- 'ATAC'
adata <- TSSEnrichment(adata)
stored_filters['median_TSSe'] = median(adata[[]][,'TSS.enrichment'])
stored_filters['median_genesperCells_RNA'] = median(adata[[]][,'nFeature_RNA'])
stored_filters['median_hq_atac_fragm'] = median(adata[[]][,'atac_fragments'])
invisible(gc())

adata@meta.data$orig.ident = sample_name

options(repr.plot.width = 10, repr.plot.height = 5)
png(file.path(output_dir,'Intermediate_metrics.png'), width =800, height = 400)
Idents(adata) <- "orig.ident"
VlnPlot(adata, features = c("nCount_ATAC", "nCount_RNA", "percent.mt",'TSS.enrichment'),  
        ncol = 4, cols=sample_color,
  log = TRUE, pt.size = 0)+ NoLegend()
dev.off()


invisible(gc())
#RNA analysis
DefaultAssay(adata) <- 'RNA'
suppressMessages(adata <- SCTransform(adata, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims=1:50, 
                                            reduction.name='umap.rna', reduction.key='rnaUMAP_'))


#ATAC analysis
#We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(adata) <- 'ATAC'
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff='q0')
adata <- RunSVD(adata)
adata <- RunUMAP(adata, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')


#Multimodal analysis
adata <- FindMultiModalNeighbors(adata, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
adata <- RunUMAP(adata, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
adata <- FindClusters(adata, graph.name='wsnn', algorithm=4, resolution = .5, verbose=FALSE)


######### Plots and stats ############
####################################

options(repr.plot.width = 12, repr.plot.height = 5)
png(file.path(output_dir,'Intermediate_UMAPs_clusters.png'), width =1000, height = 400)
p1 <- DimPlot(adata, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(adata, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(adata, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

options(repr.plot.width = 14, repr.plot.height = 8)
png(file.path(output_dir,'Intermediate_UMAPs_qc.png'), width =1000, height = 600)
p5 <- FeaturePlot(adata, features = c("CUX1",'FYB1', 'PECAM1', 'ITGA2B'),
                  reduction = 'umap.wnn', 
                  cols = c("lightgrey","darkgreen"), ncol = 4)
p6 <- FeaturePlot(adata, features = c("nCount_SCT", "nCount_RNA", "nCount_ATAC",'TSS.enrichment'), ncol = 4,
                  reduction = 'umap.wnn')
p7 <- FeaturePlot(adata, features = c("nFeature_SCT", "nFeature_RNA", "nFeature_ATAC",'percent.mt'), ncol = 4,
                  reduction = 'umap.wnn')
p6 / p7 / p5
dev.off()

adata@meta.data$DBL_comb = paste("R",adata@meta.data$scDblFinder.class, 
                                 "A",adata@meta.data$scDblFinder.class_atac, sep=":")  

png(file.path(output_dir,'Intermediate_UMAPs_doublets.png'), width =1000, height = 300)

p1<- DimPlot(adata, reduction = "umap.rna", group.by = "scDblFinder.class", label.size = 2.5, repel = TRUE) + ggtitle("RNA")+ scale_color_manual(values = c("gray", "black"))
p2<- DimPlot(adata, reduction = "umap.atac", group.by = "scDblFinder.class_atac", label.size = 2.5, repel = TRUE) + ggtitle("ATAC") + scale_color_manual(values = c("gray", "black"))
p3<-DimPlot(adata, reduction = "umap.wnn", group.by = "DBL_comb", label.size = 2.5, repel = TRUE) + ggtitle("WNN") + scale_color_manual(values = c("black", "orange", "gold","gray"))

p1 + p2 +p3
dev.off()

#### add deconvolution results
################################

larry     = read.csv(paste0(deconv_dir, sample_name, '.larry_barcodes_assignemnts.csv'), row.names=1)
larry_uni = larry[larry$No_GFPbcs==1,]
gfp1      = setNames(larry$No_GFPbcs, stringr::str_split_fixed(larry$CellBC, "\\:",3)[,3])
adata@meta.data$No_assigned_GFPbc = gfp1[rownames(adata@meta.data)]
#adata@meta.data$No_assigned_GFPbc[is.na(adata@meta.data$No_assigned_GFPbc)] <- 0
gfp2      = setNames(larry_uni$GFPbc, stringr::str_split_fixed(larry_uni$CellBC, "\\:",3)[,3])
adata@meta.data$Assigned_GFPbc = gfp2[rownames(adata@meta.data)]

gts = setNames(nm = sort(unique (larry_uni$GFPbc)), 
               object = c('16bp_del','80bp_del','80bp_del','80bp_del',
                          "HET", 'MUT','MUT','MUT', 'WT','WT','WT'))

adata@meta.data$Assigned_GFPgenotype = gts[adata@meta.data$Assigned_GFPbc]

adata_assigned = subset(adata, subset = No_assigned_GFPbc ==1 )

png(file.path(output_dir,'Intermediate_UMAPs_unique_GFP_bc.png'), width =1000, height = 600)
p1<- DimPlot(adata_assigned, reduction = "umap.rna", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("RNA") 
p2<- DimPlot(adata_assigned, reduction = "umap.atac", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPgenotype", label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 +p3 +p4 
dev.off()



png(file.path(output_dir,'CUX1_peak_coverage_plot.png'), width =600, height = 400)
options(repr.plot.width = 10, repr.plot.height = 6)
ranges.show <- StringToGRanges('chr7-101855500-101858000')
snp =  StringToGRanges('chr7-101856644-101856655')
snp$color <- "black"
CoveragePlot(adata_assigned, group.by='Assigned_GFPgenotype', region = ranges.show, region.highlight = snp,
              assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)
dev.off()

#### save object and selected barcodes

saveRDS(adata, file = file.path(output_dir,'preliminary_filtered.rds'))
filtered_bcs <- colnames(adata[["RNA"]])
write(filtered_bcs, file=(file.path(output_dir,'keep_barcodes_step1.txt')),sep='\n')


#### output metrics

stored_filters["scDBL_RNA"] <- sum(adata@meta.data$scDblFinder.class == 'doublet')
stored_filters["scDBL_ATAC"]<- sum(adata@meta.data$scDblFinder.class_atac == 'doublet')
stored_filters["scDBL_both"]<- sum(adata@meta.data$DBL_comb == 'R:doublet:A:doublet')
stored_filters["genotyped"] <- length(colnames(adata_assigned[["RNA"]]))

write.table(data.frame(stored_filters), (file.path(output_dir,'barcodes_stats.tsv')),
            sep="\t", quote=F, col.names=FALSE)










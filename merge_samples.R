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
warnLevel <- getOption('warn')
options(warn = -1)
library(ggplot2)
suppressMessages(library(scDblFinder))



samples        = c("MCO_1278","MCO_1279", "MCO_1280","MCO_1281")

adatas <- list()
for (samp in samples){
    adata <- readRDS(sprintf("/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/%s/pre_merge/pre_merged.rds",samp))
    adata@meta.data$orig.ident = samp
    adatas[[samp]] <- adata
    DefaultAssay(adatas[[samp]]) <- "RNA"
    adatas[[samp]][['ATAC']] <- NULL
    }
adatas


merged = merge(x =adatas[[1]], y=adatas[2:4], add.cell.ids = samples )
merged

merged[["RNA"]] <-JoinLayers(merged[["RNA"]])
merged

rm(adatas)
gc()

for (samp in samples){
sample_dir = sprintf("/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/%s/", samp)
lfmat      = read.table(paste0(sample_dir, "snATAC_matrices/", samp, ".long_fmt_mtx.txt.gz"))
lfmat$V1   = paste(samp,lfmat$V1, sep="_" )    
if(samp==samples[1]){
    LFM = lfmat } else {
}
    LFM = rbind(LFM, lfmat)   
}

atac_sm <- with(LFM,
                sparseMatrix(i=as.numeric(as.factor(V2)), j=as.numeric(as.factor(V1)), 
                             x=V3, dimnames=list(levels(as.factor(V2)), levels(as.factor(V1)))))

############################################################
#create the new chromatin assay object and add to Seurat object
############################################################

atac_sm       <- atac_sm[,colnames(merged)]
grange.counts <- StringToGRanges(rownames(atac_sm), sep = c(':', '-'))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_sm       <- atac_sm[as.vector(grange.use), ]
suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations)  <- 'UCSC'
genome(annotations)          <- 'hg38'

frag.file <- "/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/merged.atac_fragments.tsv.gz"
suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_sm, sep=c(':', '-'), 
                                                     genome='hg38', fragments=frag.file, 
                                                     min.cells=-1, min.features=-1, 
                                                     annotation=annotations))



merged[['ATAC']] <- chrom_assay

merged

output_dir = "/group/soranzo/paola.benaglio/k562_multiome/processing_outputs/"
saveRDS(merged, file = file.path(output_dir,'merged_unprocessed.rds'))


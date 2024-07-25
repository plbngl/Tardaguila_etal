.libPaths(new = c("/home/paola.benaglio/conda_envs/renv_multiome/lib/R/library",
          "/group/soranzo/paola.benaglio/r_libraries"))

suppressMessages(library(scDblFinder))
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
frag.file <- args[1]
barcode.file<- args[2]

repeats     <- import('/group/soranzo/paola.benaglio/references/blacklist_repeats_segdups_rmsk_hg38.bed')
otherChroms <- GRanges(c("chrM","chrX","chrY","MT"),IRanges(1L,width=10^8)) 
toExclude   <- suppressWarnings(c(repeats, otherChroms))
bcs 	    <- readLines(barcode.file)
res         <- amulet(frag.file, regionsToExclude=toExclude, barcodes=bcs)

write.table(res, "Amulet_selected_bc.tsv", sep="\t", quote=F)

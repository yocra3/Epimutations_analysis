#'#################################################################################
#'#################################################################################
#' Prepare GSE97362 dataset
#'#################################################################################
#'#################################################################################

# GSE97362 ####
library(GEOquery)
library(HDF5Array)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE97362", GSEMatrix = TRUE)[[1]]

## Convert to GenomicRatioSet
geo.full$disease <- gsub("disease state: ", "", geo.full$characteristics_ch1.3)
gse97362 <- makeGenomicRatioSetFromMatrix(exprs(geo.full), 
                                          pData = pData(geo.full))
save(gse97362, file = "results/preprocess/GSE97362/GSE97362.normalizedRaw.Rdata")


grAnnot <- readRDS("data/HM450.hg19.manifest.rds")
gse97362 <- dropMethylationLoci(gse97362)
save(gse97362, file = "results/preprocess/GSE97362/GSE97362.allCpGs.GenomicRatioSet.Rdata")

### Remove crosshibridizing and probes with SNPs
gse97362 <- gse97362[!grAnnot[rownames(gse97362)]$MASK_general, ]
save(gse97362, file = "results/preprocess/GSE97362/GSE97362.filterAnnotatedProbes.GenomicRatioSet.Rdata")

### Remove probes in sexual chromosomes
gse97362 <- gse97362[!seqnames(rowRanges(gse97362)) %in% c("chrX", "chrY"), ]
save(gse97362, file = "results/preprocess/GSE97362/GSE97362.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")


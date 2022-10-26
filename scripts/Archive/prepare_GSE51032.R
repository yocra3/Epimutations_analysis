#'#################################################################################
#'#################################################################################
#' Prepare GSE51032 dataset
#'#################################################################################
#'#################################################################################

# GSE51032 ####
library(GEOquery)
library(minfi)
library(meffil)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE51032", GSEMatrix = TRUE)[[1]]
## Download file in another server
load("data/GSE51032.Rdata")

outPrefix <- "results/preprocess/GSE51032/GSE51032"

gset <- makeGenomicRatioSetFromMatrix(exprs(gse51032), pData = pData(gse51032),
                                      array = "IlluminaHumanMethylation450k",
                                      annotation = "ilmn12.hg19")
save(gset, file = paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata"))



## Load dataset ####
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

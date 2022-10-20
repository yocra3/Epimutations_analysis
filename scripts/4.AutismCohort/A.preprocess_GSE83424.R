#'#################################################################################
#'#################################################################################
#' Prepare GSE83424
#'#################################################################################
#'#################################################################################

## Download data from GEO (bash)
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83424/matrix/GSE83424_series_matrix.txt.gz' -p data

## Load libraries
library(GEOquery)
library(minfi)
library(tidyverse)

## Load data 
geo.full <- parseGEO("data/GSE83424_series_matrix.txt.gz")

## Convert matrix to GenomicRatioSet
gse83424 <- makeGenomicRatioSetFromMatrix(exprs(geo.full), pData = pData(geo.full))

## Set clearer variable names
gse83424$sex <- gse83424$`gender:ch1`
gse83424$status <- gse83424$`status:ch1`

outPrefix <- "results/preprocess/GSE83424/GSE83424"

save(gse83424, file = paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata"))

## Remove bad probes
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gse83424)
save(gset, file = paste0(outPrefix, ".normalizedComBat.allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".normalizedComBat.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))


## Load detection p-vals from idat files ####
library(meffil)
options(mc.cores = 5)

samplesheet <- meffil.create.samplesheet("data/GSE83424", recursive=TRUE)
qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

dp <- meffil.load.detection.pvalues(qc.objects)
save(dp, file = paste0(outPrefix, ".detection_pvalues.Rdata"))

## Set to NA measurements with detection p-values > 2e-16
dp.f <- dp[rownames(gset), colnames(gset)]
beta <- assay(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

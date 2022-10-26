#'#################################################################################
#'#################################################################################
# Download GSE84727
#'#################################################################################
#'#################################################################################

# Load data ####
library(GEOquery)
library(HDF5Array)
library(minfi)
library(tidyverse)
library(meffil)

options(timeout = 10000)

## Download from GEO
geo.full <- getGEO("GSE84727", GSEMatrix = TRUE)[[1]]

beta_raw <- read_delim("data/GSE84727/GSE84727_normalisedBetas.csv.gz", delim = ",")
beta_mat <- data.matrix(beta_raw[, -1])
rownames(beta_mat) <- beta_raw$X1


pheno <- pData(geo.full)
rownames(pheno) <- pheno$`sentrixids:ch1`
colnames(beta_mat) <- pheno[colnames(beta_mat), "geo_accession"]

pheno <- pData(geo.full)
pheno <- pheno[colnames(beta_mat), ]
gset <- makeGenomicRatioSetFromMatrix(beta_mat, pData = pheno)
save(gset, file = "data/GSE84727/GSE84727_gset_raw.Rdata")

## Process gset
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = "data/GSE84727/GSE84727.allCpGs.withNA.GenomicRatioSet.Rdata")

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file = "data/GSE84727/GSE84727.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = "data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

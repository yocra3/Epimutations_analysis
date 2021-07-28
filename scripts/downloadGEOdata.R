#'#################################################################################
#'#################################################################################
# Download data
#'#################################################################################
#'#################################################################################


# GSE82273 ####
library(GEOquery)
library(minfi)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE82273", GSEMatrix = TRUE)[[1]]

geo.full <- loadSum
gse82273 <-  makeGenomicRatioSetFromMatrix(exprs(geo.full), 
                                  colData = pData(geo.full))
save(gse82273, file = "data/gse82273_raw.Rdata")

## Load from HDF5 
se <- loadHDF5SummarizedExperiment(dir = "data/GSE82273/", prefix = "GSE82273_raw")         

gset <- makeGenomicRatioSetFromMatrix(data.matrix(assay(se)), pData = colData(se))
save(gset, file = "data/GSE82273/GSE82273_gset_raw.Rdata")

## Process gset
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = "data/GSE82273/GSE82273.allCpGs.withNA.GenomicRatioSet.Rdata")

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file = "data/GSE82273/GSE82273.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = "data/GSE82273/GSE82273.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

# GSE84727 ####
library(GEOquery)
library(HDF5Array)
library(minfi)
library(tidyverse)
library(meffil)

options(timeout = 10000)

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

beta.pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
## No se ven outlier grandes de metilación
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


# GSE112611 ####
library(GEOquery)
library(minfi)
library(SummarizedExperiment)
library(tidyverse)

options(timeout = 10000)

## Get phenotypes
geo.full <- getGEO("GSE112611", GSEMatrix = TRUE)[[1]]

## Prepare GenomicRatioSet ####
raw_df <- read_delim("data/GSE112611_beta_values.txt.gz", delim = "\t")

beta_mat <- select(raw_df, -starts_with(c("Detection", "ID"))) %>%
  data.matrix()
rownames(beta_mat) <- raw_df$ID_REF

pMap <- data.frame(geo = colnames(geo.full), methy = geo.full$description)
rownames(pMap) <- pMap$methy
colnames(beta_mat) <- pMap[colnames(beta_mat), "geo"] 

gset <- makeGenomicRatioSetFromMatrix(beta_mat, pData = pData(geo.full)[colnames(beta_mat), ],
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")

outPrefix <- "data/GSE112611/GSE112611"

## Make matrix of detection p-values ####
detP <- select(raw_df, starts_with("Detection")) %>%
  data.matrix()
rownames(detP) <- raw_df$ID_REF
colnames(detP) <- colnames(beta_mat)
detP <- detP[rownames(gset), colnames(gset)]
save(gset, file = paste0(outPrefix, ".raw.GenomiRatioSet.Rdata"))
save(detP, file = paste0(outPrefix, ".detectionPvalues.Rdata"))

## Filter CpGs with call rate < 95%
detP <- detP[rowMeans(detP < 2e-16) > 0.95, ]
detP <- detP[, colMeans(detP < 2e-16) > 0.95]

gset <- gset[rownames(detP), colnames(detP) ]
save(gset, file = paste0(outPrefix, ".normalizedBeta.GenomiRatioSet.Rdata"))
save(detP, file = paste0(outPrefix, ".detectionPvaluesfilt.Rdata"))

## Process gset
grAnnot <- readRDS("data/EPIC.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file =  paste0(outPrefix,".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file = paste0(outPrefix,".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix,".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

detP <- detP[rownames(gset), colnames(gset)] 

assay(gset)[detP > 2e-16] <- NA
save(gset, file = paste0(outPrefix,".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))


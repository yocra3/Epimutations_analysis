#! Rscript
#'#################################################################################
#'#################################################################################
#' Run methylation QC data with wateRmelon to include noob+BMIQ 
#' This QC is based on a previous QC performed with meffil. Therefore, the same 
#' probes and samples will be removed.
#'#################################################################################
#'#################################################################################

# Capture arguments ####
# idatsFold <- "data/IDATs_Esteller"
phenoPath <- "results/preprocess/phenotypes/INMA_phenotypes.Rdata"
badSamplesPath <- "Esteller.removed.samples.txt"
badProbesPath <- "Esteller.removed.probes.txt"
outPrefix <- "Esteller"

# Load libraries and data ####
library(minfi)
library(dplyr)
library(wateRmelon)

load(phenoPath)
badSamples <- read.table(badSamplesPath, sep = "\t", header = TRUE)
badSamples <- unique(badSamples$sample.name)
badProbes <- read.table(badProbesPath, sep = "\t", header = TRUE)
badProbes <- unique(badProbes$name)

# # Load IDATs ####
# targets <- read.metharray.sheet(idatsFold)
# 
# ## Add phenos
# targets$SampleID <- as.character(as.numeric(substring(targets$Sample_Name, 7, 10)))
# combSheet <- left_join(targets, pheno, by = "SampleID")
# 
# RGset <- read.metharray.exp(targets = combSheet, verbose = TRUE)

## Load RGset from other server (due to problems in ISGlobal data management)
load("data/RGset.Rdata")

# Remove bad samples
RGset <- RGset[, !RGset$Sample_Name %in% badSamples]
colnames(RGset) <- RGset$Sample_Name

# Normalization: Noob + BMIQ
mset <- preprocessNoob(RGset)
bmiq <- BMIQ(mset)
gset <- makeGenomicRatioSetFromMatrix(bmiq[!rownames(bmiq) %in% badProbes, !colnames(bmiq) %in% badSamples] )
gset$SampleID <- as.character(as.numeric(substring(colnames(gset), 7, 10)))

save(gset, file = paste0(outPrefix, ".BMIQNormalization.normalizedRaw.GenomicRatioSet.Rdata"))

### Add NAs -> Run create GenomicRatioSet and then add NAs from another dataset
load(paste0(outPrefix, ".minfiRawNormalization.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
nas <- gset

com_cpgs <- intersect(rownames(nas), rownames(ori))
gset <- ori[com_cpgs, ]

assay(gset, 1)[is.na(getBeta(nas[com_cpgs, ]))] <- NA
save(gset, file = paste0(outPrefix, ".minfiRawNormalization.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))




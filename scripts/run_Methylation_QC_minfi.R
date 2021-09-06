#! Rscript
#'#################################################################################
#'#################################################################################
#' Run methylation QC data with minfi
#' This script performs different normalization algorithms included in minfi
#' This QC is based on a previous QC performed with minfi. Therefore, the same 
#' probes and samples will be removed.
#'#################################################################################
#'#################################################################################

# Capture arguments ####
idatsFold <- "data/IDATs_Esteller"
phenoPath <- "results/preprocess/phenotypes/INMA_phenotypes.Rdata"
badSamplesPath <- "Esteller.removed.samples.txt"
badProbesPath <- "Esteller.removed.probes.txt"
outPrefix <- "Esteller"

# Load libraries and data ####
library(minfi)
library(dplyr)

load(phenoPath)
badSamples <- read.table(badSamplesPath, sep = "\t", header = TRUE)
badSamples <- unique(badSamples$sample.name)
badProbes <- read.table(badProbesPath, sep = "\t", header = TRUE)
badProbes <- unique(badProbes$name)

# Load IDATs ####
targets <- read.metharray.sheet(idatsFold)

## Add phenos
targets$SampleID <- as.character(as.numeric(substring(targets$Sample_Name, 7, 10)))
combSheet <- left_join(targets, pheno, by = "SampleID")

RGset <- read.metharray.exp(targets = combSheet, verbose = TRUE)

# Remove bad samples
RGset <- RGset[, !RGset$Sample_Name %in% badSamples]
colnames(RGset) <- RGset$Sample_Name

# Normalization ####
## Raw preprocessing 
mset <- preprocessRaw(RGset)
gset <- ratioConvert(mapToGenome(mset))
save(gset, file = paste0(outPrefix, ".minfiRawNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
rm(mset, gset)

## Illumina preprocessing, as performed by Genome Studio
mset <- preprocessIllumina(RGset)
gset <- ratioConvert(mapToGenome(mset))
save(gset, file = paste0(outPrefix, ".minfiIlluminaNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
rm(mset, gset)

## SWAN normalization
mset <- preprocessSWAN(RGset)
gset <- ratioConvert(mapToGenome(mset))
save(gset, file = paste0(outPrefix, ".minfiSWANNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
rm(mset, gset)

## Quantile normalization
gset <- preprocessQuantile(RGset, fixOutliers = FALSE)
save(gset, file = paste0(outPrefix, ".minfiQuantileNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
rm(mset, gset)

## Noob
mset <- preprocessNoob(RGset)
gset <- ratioConvert(mapToGenome(mset))
save(gset, file = paste0(outPrefix, ".minfiNoobNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
rm(mset, gset)

##  Functional normalization
gset <- preprocessFunnorm(RGset)
save(gset, file = paste0(outPrefix, ".minfiFunctionalNormalization.normalizedRaw.GenomicRatioSet.Rdata"))
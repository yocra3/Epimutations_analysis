#! Rscript
#'#################################################################################
#'#################################################################################
#' Create GenomicRatioSets
#' This code prepares a GenomicRatioSet for analysis, by removing different types of
#' probes:
#' - Probes not measuring methylation (SNPs, CH probes)
#' - Crosshibridizing probes
#' - Probes with SNPs
#' - Probes in sexual chromosomes
#' 
#' For crosshibridizing probes and probes with SNPs, we used annotation from PMID: 27924034
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetfile <- args[1]
qcfile <- args[2]
manifest <- args[3]
outPrefix <- args[4]

## Load libraries ####
library(minfi)
library(meffil)

## Load dataset ####
load(gsetfile)
grAnnot <- readRDS(manifest)
load(qcfile)
ori <- gset


### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

## Create Initial and final dataset with missings
final <- gset

dp <- meffil.load.detection.pvalues(qc.objects)
gset <- ori
dp.f <- dp[rownames(gset), colnames(gset)]

beta <- assay(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".allCpGs.withNA.GenomicRatioSet.Rdata"))

gset <- gset[rownames(final), colnames(final)]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

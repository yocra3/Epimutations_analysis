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
manifest <- args[2]
outPrefix <- args[3]

## Load libraries ####
library(minfi)

## Load dataset ####
load(gsetfile)
grAnnot <- readRDS(manifest)

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".normalizedComBat.allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".normalizedComBat.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))
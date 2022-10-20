#! Rscript
#'#################################################################################
#'#################################################################################
#' Create list of common control samples
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)

load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
gset.1 <- gset

load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
gset.2 <- gset


gset.1.con <- gset.1[, gset.1$Batch == "Esteller" & !gset.1$dup]
gset.2.con <- gset.2[, gset.2$Batch == "Esteller" & !gset.2$dup]

samps <- intersect(colnames(gset.1.con), colnames(gset.2.con))
save(samps, file = "INMA.commonControlSamples.Rdata")

#! Rscript
#'#################################################################################
#'#################################################################################
#' Create list of common control samples
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetPath1 <- args[1]
gsetPath2 <- args[2]

## Load libraries
library(minfi)

load(gsetPath1)
gset.1 <- gset

load(gsetPath2)
gset.2 <- gset


gset.1.con <- gset.1[, gset.1$Batch == "Esteller" & !gset.1$dup]
gset.2.con <- gset.2[, gset.2$Batch == "Esteller" & !gset.2$dup]

samps <- intersect(colnames(gset.1.con), colnames(gset.2.con))
save(samps, file = "INMA.commonControlSamples.Rdata")

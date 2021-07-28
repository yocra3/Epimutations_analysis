#! Rscript
#'#################################################################################
#'#################################################################################
#' Adapt INMA comb datasets 
#' - select age 0 samples
#' - From MeDALL: select SAB0 samples
#' - Mark duplicates
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)

## Load GRS
load("INMA_comb.normalizedComBat.GenomicRatioSet.Rdata")

## Select Age 0 samples
gset.0 <- gset[, gset$Sample_Group == "Age_0"]

## Select INMA SAB samples
gset.sab <- gset.0[, grep("^04|SAB", gset.0$Sample_Name)]


## mark duplciated samples
dups <- gset.sab$idnum[duplicated(gset.sab$idnum)]
gset.sab$dup <- gset.sab$idnum %in% dups 
gset <- gset.sab
save(gset, file = "INMA_comb.normalizedComBat.GenomicRatioSet.corrected.Rdata")
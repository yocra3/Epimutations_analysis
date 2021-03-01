#! Rscript
#'#################################################################################
#'#################################################################################
#' Adapt INMA comb datasets 
#' - select age 0 samples
#' - From MeDALL: select SAB0 samples
#' - Mark duplicates
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gset <- args[1]

## Load libraries ####
library(minfi)

## Load GRS
load(gset)

## Select Age 0 samples
gset.0 <- gset[, gset$Sample_Group == "Age_0"]

## Select INMA SAB samples
gset.sab <- gset.0[, grep("^04|SAB", gset.0$Sample_Name)]


## mark duplciated samples
dups <- gset.sab$idnum[duplicated(gset.sab$idnum)]
gset.sab$dup <- gset.sab$idnum %in% dups 

## Remove technical duplicates
gset <- gset.sab[, !gset.sab$Sample_Name %in% c("04_561_0_duplicate", "04_287_0_duplicate", "SAB_C_0636_Rep1")]

save(gset, file = "INMA_comb.normalizedComBat.GenomicRatioSet.corrected.Rdata")
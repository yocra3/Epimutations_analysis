#! Rscript
#'#################################################################################
#'#################################################################################
#' Merge INMA datasets and select age 0 samples
#' - From MeDALL: select SAB0 samples
#' - Add 
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetMeDALL <- args[1]
gsetEsteller <- args[2]

## Load libraries ####
library(minfi)

## Load MeDALL
load(gsetMeDALL)
gsetMeDALL <- gset

## Select Sabadell Samples
gset.SAB <- gsetMeDALL[, grep("^04_", gsetMeDALL$Sample_Name)]
## Select Sabadell Samples
gset.SAB0 <- gset.SAB[, grep("_0$", gset.SAB$Sample_Name)]
gset.SAB0$batch <- "MeDALL"

## Load Estellet
load(gsetEsteller)
gsetEsteller <- gset
gsetEsteller$batch <- "Esteller"

gset <- combineArrays(gset.SAB0, gsetEsteller)

## Remove techincal duplicates
gset <- gset[, !gset$Sample_Name %in% c("SAB_C_0636_Rep1", "SAB_C_0423_Rep1")]

## mark duplciated samples
dups <- gset$idnum[duplicated(gset$idnum)]
gset$dup <- gset$idnum %in% dups 
save(gset, file = "INMA0combined.normalizedComBat.GenomicRatioSet.Rdata")
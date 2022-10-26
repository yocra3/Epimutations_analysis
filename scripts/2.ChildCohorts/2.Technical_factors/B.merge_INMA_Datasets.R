#! Rscript
#'#################################################################################
#'#################################################################################
#' Merge INMA datasets and select age 0 samples
#' - From MeDALL: select SAB0 samples
#' - Add batch variable
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)

## Load MeDALL
load("MeDALL_all.normalizedComBat.GenomicRatioSet.Rdata")
gsetMeDALL <- gset

## Select Sabadell Samples
gset.SAB <- gsetMeDALL[, grep("^04_", gsetMeDALL$Sample_Name)]
## Select Sabadell Samples
gset.SAB0 <- gset.SAB[, grep("_0$", gset.SAB$Sample_Name)]
gset.SAB0$Batch <- "MeDALL"

## Load Esteller
load("Esteller.normalizedRaw.GenomicRatioSet.Rdata")
gsetEsteller <- gset
gsetEsteller$Batch <- "Esteller"

gset <- combineArrays(gset.SAB0, gsetEsteller)

## mark duplciated samples
dups <- gset$idnum[duplicated(gset$idnum)]
gset$dup <- gset$idnum %in% dups 

save(gset, file = "INMA0combined.normalizedComBat.GenomicRatioSet.Rdata")

## Raw dataset 
## Load MeDALL
load("MeDALL_all.normalizedRaw.GenomicRatioSet.Rdata")
gsetMeDALL <- gset

## Select Sabadell Samples
gset.SAB <- gsetMeDALL[, grep("^04_", gsetMeDALL$Sample_Name)]
## Select Sabadell Samples
gset.SAB0 <- gset.SAB[, grep("_0$", gset.SAB$Sample_Name)]
gset.SAB0$Batch <- "MeDALL"

## Load Esteller
load("Esteller.normalizedRaw.GenomicRatioSet.Rdata")
gsetEsteller <- gset
gsetEsteller$Batch <- "Esteller"

gset <- combineArrays(gset.SAB0, gsetEsteller)

## mark duplciated samples
dups <- gset$idnum[duplicated(gset$idnum)]
gset$dup <- gset$idnum %in% dups 

save(gset, file = "INMA0combined.normalizedRaw.GenomicRatioSet.Rdata")
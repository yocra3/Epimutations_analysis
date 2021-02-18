#'#################################################################################
#'#################################################################################
#' Configuration for applying ComBat to all samples of MeDALL
#'#################################################################################
#'#################################################################################

## Output name
outPrefix <- "MeDALL_all"
GRSfile <- paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata")

## Select variable fo adjust
batch.var <- "TEM"
baseModel <- formula(~ pred.sex + Sample_Group)
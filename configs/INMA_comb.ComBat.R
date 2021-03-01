#'#################################################################################
#'#################################################################################
#' Configuration for applying ComBat to all samples of MeDALL
#'#################################################################################
#'#################################################################################

## Output name
outPrefix <- "INMA_comb"
GRSfile <- paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata")

## Select variable fo adjust
batch.var <- "TechVar"
baseModel <- formula(~ pred.sex + Sample_Group)
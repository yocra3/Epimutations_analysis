#! Rscript
#'#################################################################################
#'#################################################################################
#' Create report after ComBat
#'#################################################################################
#'#################################################################################

args <- commandArgs(trailingOnly=TRUE)
configFile <- args[1]
outPrefix <- args[2]

## Load libraries
library(meffil)
library(minfi)

load(paste0(outPrefix, ".normalizedComBat.GenomicRatioSet.Rdata"))
load(paste0(outPrefix, ".norm.obj.pc.Rdata"))

beta.pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(pcs),
  batch.pcs = seq_len(pcs),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.postComBat.html"))

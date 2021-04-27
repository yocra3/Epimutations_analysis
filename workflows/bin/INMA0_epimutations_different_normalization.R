#! Rscript
#'#################################################################################
#'#################################################################################
#' Run epimutations analysis on INMA age 0.
#' Input: 
#' - path to normalized GenomicRatioSet with NAs for failing probes
#' - Name of the dataset to name the output files
#' - Common samples between normalizations
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetPath <- args[1]
outPrefix <- args[2]
samps <- args[3]

## Load libraries
library(minfi)
library(epimutacions)

load(gsetPath)
load(samps)

## Select CpGs in candidate regions
gset.sel <- subsetByOverlaps(gset, candRegsGR)

### Make cases and controls datasets
controls <- gset.sel[, samps]
cases <- gset.sel[, !colnames(gset.sel) %in% samps]


methods <- c("manova", "mlm", "isoforest", "mahdistmcd", "barbosa", "beta")
names(methods) <- methods
res <- lapply(methods, epimutations, case_samples = cases, control_panel = controls)
save(res, file = paste0(outPrefix, ".epimutations.INMA0.duplicates.Rdata"))


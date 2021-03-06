#! Rscript
#'#################################################################################
#'#################################################################################
#' Run epimutations analysis on INMA age 0.
#' This code performs quality control to methylation data. 
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

## Load libraries
library(minfi)
library(epimutacions)

load(gsetPath)

if (length(args) > 2){
  samps <- args[3]
  load(samps)
} else {
  samps <- NULL
}

## Select samples from Esteller and duplicates
gset.sel <- gset[, gset$Batch == "Esteller" | gset$dup]

## Select CpGs in candidate regions
gset.sel <- subsetByOverlaps(gset.sel, candRegsGR)

### Make cases and controls datasets
cases <- gset.sel[, gset.sel$dup]
controls <- gset.sel[, !gset.sel$dup]

if(!is.null(samps)){
  controls <- controls[, samps]
}

methods <- c("manova", "mlm", "isoforest", "mahdistmcd", "barbosa", "beta")
names(methods) <- methods
res <- lapply(methods, epimutations, case_samples = cases, control_panel = controls)
save(res, file = paste0(outPrefix, ".epimutations.INMA0.duplicates.Rdata"))


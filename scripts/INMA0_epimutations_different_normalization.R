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

## Load libraries
library(minfi)
library(epimutacions)
library(meffil)
library(parallel)

load("INMA.commonControlSamples.Rdata")
methods <- names(epi_parameters())
names(methods) <- methods

runEpimutations <- function(gsetPath){
  load(gsetPath)
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset, candRegsGR)
  
  ### Make cases and controls datasets
  controls <- gset.sel[, samps]
  cases <- gset.sel[, !colnames(gset.sel) %in% samps]
  
  
  res <- mclapply(methods, epimutations, case_samples = cases, control_panel = controls, 
                  mc.cores = 6)
}
paths <- dir(pattern = "Esteller.*Normalization")
names(paths) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))

INMA_norm <- lapply(paths, runEpimutations)
names(INMA_norm) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))
save(INMA_norm, file = "results/epimutations/Esteller.epimutations.normalization.Rdata")

getResiduals <- function(gset){
  pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
  m <- getM(gset)
  res <- residuals(lmFit(m, pcs[, 1:10]), m)
  beta <- ilogit2(res)
  assay(gset) <- beta
  gset
}

runEpimutationsSet <- function(gset){
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset, candRegsGR)
  
  ### Make cases and controls datasets
  controls <- gset.sel[, samps]
  cases <- gset.sel[, !colnames(gset.sel) %in% samps]
  
  
  res <- mclapply(methods, epimutations, case_samples = cases, control_panel = controls, mc.cores = 6)
}

gset.residuals <- lapply(paths, function(x){
  load(x)
  getResiduals(gset)
})
names(gset.residuals) <- names(paths)
save(gset.residuals, file = "Esteller.allminfiNormalizations.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

INMA_norm_resid <- lapply(gset.residuals, runEpimutationsSet)

## Rerun execution with error
filterSet <- function(gset){
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset, candRegsGR)
  
  gset.sel <- gset.sel[rowMeans(is.na(getBeta(gset.sel))) < 0.05, ]
  
  ### Make cases and controls datasets
  controls <- gset.sel[, samps]
  cases <- gset.sel[, !colnames(gset.sel) %in% samps]
  list(controls = controls, cases = cases)
  } 
illumina.set <- filterSet(gset.residuals$IlluminaNormalization)
INMA_norm_resid$IlluminaNormalization$manova <- epimutations(case_samples = illumina.set$cases, 
                                                      control_panel = illumina.set$controls,
                                                      method = "manova")
save(INMA_norm_resid, file = "results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")

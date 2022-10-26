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

runEpimutations <- function(gsetPath, params = epi_parameters()){
  load(gsetPath)
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset, candRegsGR)
  
  ### Make cases and controls datasets
  controls <- gset.sel[, samps]
  cases <- gset.sel[, !colnames(gset.sel) %in% samps]
  
  
  res <- mclapply(methods, epimutations, case_samples = cases, control_panel = controls, 
                  epi_params = params, mc.cores = 6)
}
paths <- dir(pattern = "Esteller.*Normalization.*autosomic.*withNA.*")
names(paths) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))

INMA_norm <- lapply(paths, runEpimutations)
names(INMA_norm) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))
save(INMA_norm, file = "results/epimutations/Esteller.epimutations.normalization.Rdata")

params2 <- epi_parameters()
params2$manova$pvalue_cutoff <- 0.4
params2$mlm$pvalue_cutoff <- 0.4
params2$isoforest$outlier_score_cutoff <- 0.65
params2$quantile$offset_abs <- 0.1
params2$beta$pvalue_cutoff <- 1e-5
INMA_norm_lib <- lapply(paths, runEpimutations, params = params2)
save(INMA_norm_lib, file = "results/epimutations/Esteller.epimutations.normalization.liberal.Rdata")


getResiduals <- function(gset){
  beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)
  ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
  pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
  m <- getM(gset)
  res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
  beta <- ilogit2(res)
  assay(gset, "Beta") <- beta
  gset
}

runEpimutationsSet <- function(gset, params = epi_parameters()){
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset, candRegsGR)
  
  ### Make cases and controls datasets
  controls <- gset.sel[, samps]
  cases <- gset.sel[, !colnames(gset.sel) %in% samps]
  
  
  res <- mclapply(methods, epimutations, case_samples = cases, control_panel = controls, 
                  epi_params = params, mc.cores = 6)
}

gset.residuals <- lapply(paths, function(x){
  load(x)
  getResiduals(gset)
})
names(gset.residuals) <- names(paths)
save(gset.residuals, file = "Esteller.allminfiNormalizations.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

INMA_norm_resid <- lapply(gset.residuals, runEpimutationsSet)
save(INMA_norm_resid, file = "results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")

INMA_norm_resid_lib <- lapply(gset.residuals, runEpimutationsSet, params2)
save(INMA_norm_resid_lib, file = "results/epimutations/Esteller.epimutations.normalization.residuals.liberal.Rdata")


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

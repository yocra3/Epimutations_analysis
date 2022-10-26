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


## Load libraries ####
library(minfi)
library(epimutacions)
library(parallel)
library(meffil)

load("INMA.commonControlSamples.Rdata")
load("INMA_comb.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gsetcomb <- gset
load("INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gset0 <- gset

# Crude model ####
methods <- names(epi_parameters())
names(methods) <- methods
runEpimutations <- function(gset, samps, params = epi_parameters()){
  ## Select samples from Esteller and duplicates
  gset.sel <- gset[, gset$Batch == "Esteller" | gset$dup]
  
  ## Select CpGs in candidate regions
  gset.sel <- subsetByOverlaps(gset.sel, candRegsGR)
  
  ### Make cases and controls datasets
  cases <- gset.sel[, gset.sel$dup]
  controls <- gset.sel[, samps]
  
  

  res <- mclapply(methods, epimutations, case_samples = cases, control_panel = controls, 
                  epi_params = params, mc.cores = 6)
}
res_indep <- runEpimutations(gset0, samps)
save(res_indep, file = "results/epimutations/INMA0combined.raw.epimutations.INMA0.duplicates.Rdata")


# Residuals ####
getResiduals <- function(gset){
  beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)
  ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
  pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
  m <- getM(gset)
  res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
  beta <- ilogit2(res)
  assay(gset) <- beta
  gset
}
gsetcomb_res <- getResiduals(gsetcomb)
save(gsetcomb_res, file = "INMA_comb.raw.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

gset0_res <- getResiduals(gset0)
save(gset0_res, file = "INMA0combined.raw.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
res_indep_resid <- runEpimutations(gset0_res, samps)
save(res_indep_resid, file = "results/epimutations/INMA0combined.raw.epimutations.INMA0.duplicates.residuals.Rdata")
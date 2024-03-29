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
load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gsetcomb <- gset
load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
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
res_joint <- runEpimutations(gsetcomb, samps)
save(res_joint, file = "results/epimutations/INMA_comb.epimutations.INMA0.duplicates.Rdata")

res_indep <- runEpimutations(gset0, samps)
save(res_indep, file = "results/epimutations/INMA0combined.epimutations.INMA0.duplicates.Rdata")

# Less restrictive parameters ####
params2 <- epi_parameters()
params2$manova$pvalue_cutoff <- 0.4
params2$mlm$pvalue_cutoff <- 0.4
params2$isoforest$outlier_score_cutoff <- 0.65
params2$quantile$offset_abs <- 0.1
params2$beta$pvalue_cutoff <- 1e-5
res_joint_lib <- runEpimutations(gsetcomb, samps, params2)
save(res_joint_lib, file = "results/epimutations/INMA_comb.epimutations.INMA0.duplicates.liberal.Rdata")

res_indep_lib <- runEpimutations(gset0, samps, params2)
save(res_indep_lib, file = "results/epimutations/INMA0combined.epimutations.INMA0.duplicates.liberal.Rdata")


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
save(gsetcomb_res, file = "INMA_comb.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
res_joint_resid <- runEpimutations(gsetcomb_res, samps)
save(res_joint_resid, file = "results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.Rdata")

gset0_res <- getResiduals(gset0)
save(gset0_res, file = "INMA0combined.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
res_indep_resid <- runEpimutations(gset0_res, samps)
save(res_indep_resid, file = "results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.Rdata")

# Residuals and Less restrictive parameters ####
res_joint_resid_lib <- runEpimutations(gsetcomb_res, samps, params2)
save(res_joint_resid_lib, file = "results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.liberal.Rdata")

res_indep_resid_lib <- runEpimutations(gset0_res, samps, params2)
save(res_indep_resid_lib, file = "results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.liberal.Rdata")

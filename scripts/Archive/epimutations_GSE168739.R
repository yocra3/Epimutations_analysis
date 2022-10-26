#'#################################################################################
#'#################################################################################
#' Detect epimutations in COVID dataset
#'#################################################################################
#'#################################################################################



## Load libraries
library(minfi)
library(meffil)
library(epimutacions)
library(tidyverse)
library(cowplot)

## Load data
load("results/preprocess/GSE168739/GSE168739.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Compute residuals of pcs
beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
m <- getM(gset)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(gset) <- beta
save(gset, file = "results/preprocess/GSE168739/GSE168739.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")


covid.residuals.loo <- epimutations_one_leave_out(methy = gset, method = "quantile",
  BPPARAM = BiocParallel::MulticoreParam(2))
save(res.helix.residuals.list, file = "results/epimutations/HELIX.epimutations.allSamples.residuals.Rdata")

## Sex ####
helix.boys <- helix[, helix$e3_sex == "male"]
res.helix.boys.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = helix.boys,
                              BPPARAM = MulticoreParam(10))
save(res.helix.boys.residuals.list, file = "results/epimutations/HELIX.epimutations.boys.residuals.Rdata")

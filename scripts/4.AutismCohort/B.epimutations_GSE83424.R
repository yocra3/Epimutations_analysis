#'#################################################################################
#'#################################################################################
#' Detect epimutations in GSE83424
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(meffil)
library(epimutacions)
library(tidyverse)
library(cowplot)

## Load data
load("results/preprocess/GSE83424/GSE83424.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Compute residuals of pcs
beta <- getBeta(gset)
beta[beta <= 1e-3] <- 1e-3
beta[beta >= 1 - 1e-3] <- 1 - 1e-3
assay(gset) <- beta

beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)
m <- getM(gset)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(gset) <- beta
save(gset, file = "results/preprocess/GSE83424/GSE83424.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

## Leave-one-out ####
gse83424.residuals.loo <- epimutations_one_leave_out(methy = gset, method = "quantile",  BPPARAM = BiocParallel::MulticoreParam(5))
save(gse83424.residuals.loo, file = "results/epimutations/GSE83424.epimutations.loo.residuals.Rdata")

## case-control ####
gse83424.residuals.cc <- epimutations(gset[, gset$status== "Case"], gset[, gset$status== "Control"], method = "quantile")
save(gse83424.residuals.cc, file = "results/epimutations/GSE83424.epimutations.casecontrol.residuals.Rdata")

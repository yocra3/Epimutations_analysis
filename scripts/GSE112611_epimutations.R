#'#################################################################################
#'#################################################################################
#' Detect epimutations in GSE112611
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(meffil)
library(epimutacions)
library(tidyverse)
library(cowplot)

## Load data
load("data/GSE112611/GSE112611.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

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
save(gset, file = "results/preprocess/GSE112611/GSE112611.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

## case-control ####
gset$diagnosis <- gset$`diagnosis:ch1`

gse112611.residuals.cc <- epimutations(gset[, gset$diagnosis== "Crohn's disease"], gset[, gset$diagnosis== "non-IBD control"], method = "quantile")
save(gse112611.residuals.cc, file = "results/epimutations/GSE112611.epimutations.casecontrol.residuals.Rdata")

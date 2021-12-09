#'#################################################################################
#'#################################################################################
#' Detect epimutations in CHARGE/Kabuki dataset
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)

load("results/preprocess/GSE97362/GSE97362.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Compute residuals of pcs
# beta <- meffil:::impute.matrix(getBeta(gse97362), margin = 1)
# ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
# pcs <- meffil.methylation.pcs(getBeta(gse97362), probe.range = 40000)
# m <- getM(gse97362)
# res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
# beta <- ilogit2(res)
# assay(helix) <- beta
# save(helix, file = "results/epimutations/HELIX.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
# 

gse97362 <- gse97362[rowMeans(is.na(getBeta(gse97362))) < 0.1, ]

methods <- c("beta", "quantile", "mlm")
names(methods) <- methods

case <- gse97362[, gse97362$disease != "Control"]
control <- gse97362[, gse97362$disease == "Control"]

res.gse97362.list <- lapply(methods, epimutations, case_samples = case, 
                            control_panel = control)
save(res.helix.residuals.list, file = "results/epimutations/HELIX.epimutations.allSamples.residuals.Rdata")

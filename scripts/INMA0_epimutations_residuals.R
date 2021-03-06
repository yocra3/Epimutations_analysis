#'#################################################################################
#'#################################################################################
#' Define epimutations in INMA samples at 0 years
#' - Stratify by boys and girls
#' - Stratify by batch
#'#################################################################################
#'#################################################################################

library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)
library(limma)

load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Remove duplicates
inma0 <- gset[, !(gset$dup & gset$Batch == "MeDALL")]

## Compute residuals of pcs
pcs <- meffil.methylation.pcs(getBeta(inma0), probe.range = 40000)
m <- getM(inma0)
res <- residuals(lmFit(m, pcs[, 1:10]), m)
beta <- ilogit2(res)
assay(inma0) <- beta
save(inma0, file = "INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")


methods <- c("beta", "barbosa", "mlm", "manova")
names(methods) <- methods

res.inma0.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.residuals.list, file = "INMA_comb.epimutations.allSamples.residuals.Rdata")

## Sex ####
inma0.boys <- inma0[, inma0$Sex == "M"]
res.inma0.boys.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.boys, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.boys.residuals.list, file = "INMA_comb.epimutations.boys.residuals.Rdata")


inma0.girls <- inma0[, inma0$Sex == "F"]
res.inma0.girls.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.girls, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.girls.residuals.list, file = "INMA_comb.epimutations.girls.residuals.Rdata")


## Batch ####
inma0.esteller <- inma0[, inma0$Batch == "Esteller"]
res.inma0.esteller.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.esteller, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.esteller.residuals.list, file = "INMA_comb.epimutations.esteller.residuals.Rdata")


inma0.medall <- inma0[, inma0$Batch == "MeDALL"]
res.inma0.medall.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.medall, 
                               BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.medall.residuals.list, file = "INMA_comb.epimutations.medall.residuals.Rdata")


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

load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Remove duplicates
inma0 <- gset[, !(gset$dup & gset$Batch == "MeDALL")]

methods <- c("beta", "barbosa", "mlm", "manova")
names(methods) <- methods

res.inma0.list <- lapply(methods, epimutations_one_leave_out, methy = inma0, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.list, file = "INMA_comb.epimutations.allSamples.Rdata")

## Sex ####
inma0.boys <- inma0[, inma0$Sex == "M"]
res.inma0.boys.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.boys, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.boys.list, file = "INMA_comb.epimutations.boys.Rdata")


inma0.girls <- inma0[, inma0$Sex == "F"]
res.inma0.girls.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.girls, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.girls.list, file = "INMA_comb.epimutations.girls.Rdata")


## Batch ####
inma0.esteller <- inma0[, inma0$Batch == "Esteller"]
res.inma0.esteller.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.esteller, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.esteller.list, file = "INMA_comb.epimutations.esteller.Rdata")


inma0.medall <- inma0[, inma0$Batch == "MeDALL"]
res.inma0.medall.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.medall, 
                               BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.medall.list, file = "INMA_comb.epimutations.medall.Rdata")


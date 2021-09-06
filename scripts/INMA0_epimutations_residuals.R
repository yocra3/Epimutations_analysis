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
### Same batch replicates
gset.un <- gset[ , !grepl("duplicate|Rep1", colnames(gset))]

### Different batch replicates
dups <- gset.un$idnum[duplicated(gset.un$idnum)]
gset.un$dup2 <- gset.un$idnum %in% dups 
inma0 <- gset.un[, !(gset.un$dup2 & gset.un$Batch == "MeDALL")]

## Compute residuals of pcs
beta <- meffil:::impute.matrix(getBeta(inma0), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(inma0), probe.range = 40000)
m <- getM(inma0)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(inma0) <- beta
save(inma0, file = "INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

methods <- names(epi_parameters())
names(methods) <- methods
methods <- methods[methods != "mahdistmcd"]

res.inma0.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0, 
                         BPPARAM = MulticoreParam(7, progressbar = TRUE))
save(res.inma0.residuals.list, file = "results/epimutations/INMA_comb.epimutations.allSamples.residuals.Rdata")

## Sex ####
inma0.boys <- inma0[, inma0$Sex == "M"]
res.inma0.boys.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.boys, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.boys.residuals.list, file = "results/epimutations/INMA_comb.epimutations.boys.residuals.Rdata")


inma0.girls <- inma0[, inma0$Sex == "F"]
res.inma0.girls.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.girls, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.girls.residuals.list, file = "results/epimutations/INMA_comb.epimutations.girls.residuals.Rdata")


## Batch ####
inma0.esteller <- inma0[, inma0$Batch == "Esteller"]
res.inma0.esteller.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.esteller, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.esteller.residuals.list, file = "INMA_comb.epimutations.esteller.residuals.Rdata")


inma0.medall <- inma0[, inma0$Batch == "MeDALL"]
res.inma0.medall.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0.medall, 
                               BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma0.medall.residuals.list, file = "INMA_comb.epimutations.medall.residuals.Rdata")


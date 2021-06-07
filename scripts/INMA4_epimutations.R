#'#################################################################################
#'#################################################################################
#' Define epimutations in INMA samples at 4 years
#' - Stratify by boys and girls
#' - Stratify by batch
#'#################################################################################
#'#################################################################################

library(minfi)
library(epimutacions)
library(BiocParallel)

load("MeDALL_all.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Select 4 years samples
gset.4 <- gset[, gset$Sample_Group == "Age_4"]
              
## Select INMA SAB samples
gset.sab <- gset.4[, grep("^04|SAB", gset.4$Sample_Name)]

## Remove duplicates
inma4 <- gset.sab[, !duplicated(gset.sab$idnum)]

methods <- c("beta", "barbosa", "mlm", "manova")
names(methods) <- methods

res.inma4.list <- lapply(methods, epimutations_one_leave_out, methy = inma4, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.list, file = "INMA4.epimutations.allSamples.Rdata")

## Sex ####
inma4.boys <- inma4[, inma4$Sex == "M"]
res.inma4.boys.list <- lapply(methods, epimutations_one_leave_out, methy = inma4.boys, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.boys.list, file = "INMA4.epimutations.boys.Rdata")


inma4.girls <- inma4[, inma4$Sex == "F"]
res.inma4.girls.list <- lapply(methods, epimutations_one_leave_out, methy = inma4.girls, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.girls.list, file = "INMA4.epimutations.girls.Rdata")
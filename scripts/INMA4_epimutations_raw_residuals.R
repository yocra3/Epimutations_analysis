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
library(meffil)

load("MeDALL_all.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Select 4 years samples
gset.4 <- gset[, gset$Sample_Group == "Age_4"]
              
## Select INMA SAB samples
gset.sab <- gset.4[, grep("^04|SAB", gset.4$Sample_Name)]

## Remove duplicates
inma4 <- gset.sab[, !duplicated(gset.sab$idnum)]

beta <- meffil:::impute.matrix(getBeta(inma4), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(inma4), probe.range = 40000)
m <- getM(inma4)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(inma4) <- beta
save(inma4, file = "INMA4.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

methods <- c("quantile", "beta", "mlm")
names(methods) <- methods

res.inma4.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma4, 
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.residuals.list, file = "results/epimutations/INMA4.epimutations.raw.allSamples.residuals.Rdata")

## Sex ####
inma4.boys <- inma4[, inma4$Sex == "M"]
res.inma4.boys.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma4.boys,
                         BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.boys.residuals.list, file = "results/epimutations/INMA4.epimutations.boys.residuals.Rdata")


inma4.girls <- inma4[, inma4$Sex == "F"]
res.inma4.girls.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma4.girls, 
                              BPPARAM = MulticoreParam(10, progressbar = TRUE))
save(res.inma4.girls.residuals.list, file = "results/epimutations/INMA4.epimutations.girls.residuals.Rdata")


### Use as reference the other sex
res.inma4.boys.girlsref.residuals.list <- mclapply(methods, epimutations, 
                                                   case_samples = inma4.boys, 
                                                   control_panel = inma4.girls, mc.cores = 3)
save(res.inma4.boys.girlsref.residuals.list, file = "results/epimutations/INMA4combined.raw.epimutations.boys.girlsref.residuals.Rdata")

res.inma4.girls.boysref.residuals.list <- mclapply(methods, epimutations, 
                                                   case_samples = inma4.girls, 
                                                   control_panel = inma4.boys, mc.cores = 3)
save(res.inma4.girls.boysref.residuals.list, file = "results/epimutations/INMA4combined.raw.epimutations.girls.boysref.residuals.Rdata")


## Smoking ####
res.inma4.smoking.residuals.list <- mclapply(methods, epimutations, 
                                             case_samples = inma4[ , !is.na(inma4$msmk) & inma4$msmk != "no smoking"], 
                                             control_panel =  inma4[ , !is.na(inma4$msmk) &  inma4$msmk == "no smoking"],
                                             mc.cores = 3)
save(res.inma4.smoking.residuals.list, file = "results/epimutations/INMA4combined.raw.epimutations.smoking.residuals.Rdata")

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

load("INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

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
save(inma0, file = "INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

methods <- c("quantile", "beta", "mlm")
names(methods) <- methods

res.inma0.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = inma0, 
                         BPPARAM = MulticoreParam(7, progressbar = TRUE))
save(res.inma0.residuals.list, file = "results/epimutations/INMA0combined.raw.epimutations.allSamples.residuals.Rdata")


### Use as reference the other sex
res.inma0.boys.girlsref.residuals.list <- mclapply(methods, epimutations, 
                                                   case_samples = inma0.boys, 
                                                   control_panel = inma0.girls, mc.cores = 3)
save(res.inma0.boys.girlsref.residuals.list, file = "results/epimutations/INMA0combined.raw.epimutations.boys.girlsref.residuals.Rdata")

res.inma0.girls.boysref.residuals.list <- mclapply(methods, epimutations, 
                                                   case_samples = inma0.girls, 
                                                   control_panel = inma0.boys, mc.cores = 3)
save(res.inma0.girls.boysref.residuals.list, file = "results/epimutations/INMA0combined.raw.epimutations.girls.boysref.residuals.Rdata")

## Smoking ####
res.inma0.df$smoking <- ifelse(colData(inma0)[res.inma0.df$sample, "msmk"] == "no smoking", "no", "yes")

res.inma0.smoking.residuals.list <- mclapply(methods, epimutations, 
                                                   case_samples = inma0[ , !is.na(inma0$msmk) & inma0$msmk != "no smoking"], 
                                                   control_panel =  inma0[ , !is.na(inma0$msmk) &  inma0$msmk == "no smoking"],
                                             mc.cores = 3)
save(res.inma0.smoking.residuals.list, file = "results/epimutations/INMA0combined.raw.epimutations.smoking.residuals.Rdata")

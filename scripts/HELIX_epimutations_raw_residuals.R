#'#################################################################################
#'#################################################################################
#' Define epimutations in HELIX
#' - Apply residuals
#' - Stratify by boys and girls
#' - Stratify by cohort
#'#################################################################################
#'#################################################################################

library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)

load("results/preprocess/HELIX/HELIX.raw.withNA.GenomicRatioSet.Rdata")

## Remove MOBA and non-European
helix <- gset[, gset$cohort != "MOBA" & gset$h_ethnicity_3cat == "WhiteEur_WhiteOther"]

## Compute residuals of pcs
beta <- meffil:::impute.matrix(getBeta(helix), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(helix), probe.range = 40000)
m <- getM(helix)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(helix) <- beta
save(helix, file = "results/epimutations/HELIX.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")


methods <- c("quantile", "beta", "mlm")
names(methods) <- methods

res.helix.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = helix, 
                         BPPARAM = MulticoreParam(5))
save(res.helix.residuals.list, file = "results/epimutations/HELIX.epimutations.raw.allSamples.residuals.Rdata")

## Sex ####
helix.boys <- helix[, helix$e3_sex == "male"]
res.helix.boys.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = helix.boys, 
                              BPPARAM = MulticoreParam(10))
save(res.helix.boys.residuals.list, file = "results/epimutations/HELIX.epimutations.raw.boys.residuals.Rdata")


helix.girls <- helix[, helix$e3_sex == "female"]
res.helix.girls.residuals.list <- lapply(methods, epimutations_one_leave_out, methy = helix.girls, 
                               BPPARAM = MulticoreParam(10))
save(res.helix.girls.residuals.list, file = "results/epimutations/HELIX.epimutations.raw.girls.residuals.Rdata")

## Cohort
cohort <- unique(helix$cohort)
names(cohort) <- cohort

res.helix.cohort.residuals.list <- lapply(cohort, function(x){
  helix.sub <- helix[, helix$cohort == x]
  lapply(methods, epimutations_one_leave_out, methy = helix.sub, 
         BPPARAM = MulticoreParam(10))
})
save(res.helix.cohort.residuals.list, file = "results/epimutations/HELIX.epimutations.raw.cohort.residuals.Rdata")

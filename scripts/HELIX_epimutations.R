#'#################################################################################
#'#################################################################################
#' Define epimutations in HELIX
#' - Stratify by boys and girls
#' - Stratify by cohort
#'#################################################################################
#'#################################################################################

library(minfi)
library(epimutacions)
library(BiocParallel)

load("results/preprocess/HELIX/HELIX.withNA.GenomicRatioSet.Rdata")

## Remove MOBA and non-European
helix <- gset[, gset$cohort != "MOBA" & gset$h_ethnicity_3cat == "WhiteEur_WhiteOther"]

methods <- c("beta", "barbosa", "mlm", "manova")
names(methods) <- methods

res.helix.list <- lapply(methods, epimutations_one_leave_out, methy = helix, 
                         BPPARAM = MulticoreParam(3))
save(res.helix.list, file = "results/epimutations/HELIX.epimutations.allSamples.Rdata")

## Sex ####
helix.boys <- helix[, helix$e3_sex == "male"]
res.helix.boys.list <- lapply(methods, epimutations_one_leave_out, methy = helix.boys, 
                              BPPARAM = MulticoreParam(3))
save(res.helix.boys.list, file = "results/epimutations/HELIX.epimutations.boys.Rdata")


helix.girls <- helix[, helix$e3_sex == "female"]
res.helix.girls.list <- lapply(methods, epimutations_one_leave_out, methy = helix.girls, 
                               BPPARAM = MulticoreParam(10))
save(res.helix.girls.list, file = "results/epimutations/HELIX.epimutations.girls.Rdata")

## Cohort
cohort <- unique(helix$cohort)
names(cohort) <- cohort

res.helix.cohort.list <- lapply(cohort, function(x){
  helix.sub <- helix[, helix$cohort == x]
  lapply(methods, epimutations_one_leave_out, methy = helix.sub, 
         BPPARAM = MulticoreParam(10))
})
save(res.helix.cohort.list, file = "results/epimutations/HELIX.epimutations.cohort.Rdata")

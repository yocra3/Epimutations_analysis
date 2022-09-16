#'#################################################################################
#'#################################################################################
#' Apply ComBat to GSE168739
#' ComBat will be applied to M-values
#' The result will be stored as a GenomicRatioSet
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(minfi)
library(sva)

load("results/preprocess/GSE168739/GSE168739.normalizedRaw.GenomicRatioSet.Rdata")

pheno <- colData(gset)
batch <- pheno$Batch

## Create ComBat model
modcombat <- model.matrix(~ Disease + Sex + age, data = pheno)

## Run comBat
m <- getM(gset)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots = FALSE)

beta <- ilogit2(combat_M)
assay(gset) <- beta
save(gset, file = "results/preprocess/GSE168739/GSE168739.normalizedComBat.GenomicRatioSet.Rdata")

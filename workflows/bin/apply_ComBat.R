#'#################################################################################
#'#################################################################################
#' Apply ComBat to a GenomicRatioSet
#' ComBat will be applied to M-values
#' The result will be stored as a GenomicRatioSet
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
configFile <- args[1]

## Load libraries ####
library(minfi)
library(sva)

source(configFile)
load(GRSfile)

pheno <- colData(gset)
batch <- pheno[[batch.var]]

## Create ComBat model
modcombat <- model.matrix(baseModel, data = pheno)

## Run comBat
m <- getM(gset)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots = FALSE)

beta <- ilogit2(combat_M)
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".normalizedComBat.GenomicRatioSet.Rdata"))

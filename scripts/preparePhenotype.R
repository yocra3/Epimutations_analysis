#'#################################################################################
#'#################################################################################
#' Prepare phenotype data for meffil pipeline
#'#################################################################################
#'#################################################################################

## Load libraries
library(foreign)

pheno <- read.dta("data/INMA_meth_pheno.dta")
pheno$Sex <- pheno$sex
pheno$SampleID <- pheno$idnum
save(pheno, file = "results/preprocess/phenotypes/INMA_phenotypes.Rdata")
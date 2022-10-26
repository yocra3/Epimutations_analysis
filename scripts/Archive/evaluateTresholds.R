#'#################################################################################
#'#################################################################################
#' Evaluate the thresholds in samples with replicates in Esteller
#'#################################################################################
#'#################################################################################


## Load libraries ####
library(minfi)
library(epimutacions)
library(parallel)
library(meffil)

load("INMA.commonControlSamples.Rdata")
load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

# Crude model ####
methods <- names(epi_parameters())
names(methods) <- methods

gset.sel <- gset[, gset$Batch == "Esteller" | gset$dup]
## Select CpGs in candidate regions
gset.sel <- subsetByOverlaps(gset.sel, candRegsGR)
  
### Make cases and controls datasets
cases <- gset.sel[, gset.sel$idnum %in% c(423, 636)]
controls <- gset.sel[, samps]
  
methods1 <- c("mlm", "manova", "isoforest")
params2 <- epi_parameters()
params2$manova$pvalue_cutoff <- 1
params2$mlm$pvalue_cutoff <- 1
params2$isoforest$outlier_score_cutoff <- 0

res_thres1 <- mclapply(methods1, epimutations, case_samples = cases, control_panel = controls, 
                  epi_params = params2, mc.cores = 3)
names(res_thres1) <- methods1

## mlm
thres <- c(0.005, 0.01, seq(0.05, 0.80, 0.05))
data.frame(thres, n = sapply(thres, function(x) sum(res_thres1$mlm$adj_pvalue < x)))
data.frame("-log10" = 1:9, n = sapply(1:9, function(x) sum(res_thres1$mlm$pvalue < 10^-x)))
### Seleccionar 0.4

## manova
data.frame(thres, n = sapply(thres, function(x) sum(res_thres1$manova$adj_pvalue < x)))
data.frame("-log10" = 1:9, n = sapply(1:9, function(x) sum(res_thres1$manova$pvalue < 10^-x)))
### Seleccionar 0.4

## isoforest
thresi <-  seq(0.5, 1, 0.05)
data.frame(thresi, n = sapply(thresi, function(x) sum(res_thres1$isoforest$outlier_score > x)))


params2$quantile$offset_abs <- 0.05
res_quant <- list()
res_quant$ori <-  epimutations(method = "quantile", case_samples = cases, control_panel = controls)
res_quant$d05 <- epimutations(method = "quantile", case_samples = cases, control_panel = controls, 
                       epi_params = params2)
params2$quantile$offset_abs <- 0.1
res_quant$d10 <- epimutations(method = "quantile", case_samples = cases, control_panel = controls, 
                              epi_params = params2)

sum(res_quant$d05$chromosome != 0)
sum(res_quant$d10$chromosome != 0)
sum(res_quant$ori$chromosome != 0)
## Seleccionar 0.1

## beta
params2$beta$pvalue_cutoff <- 1e-5
res_beta <- list()
res_beta$ori <-  epimutations(method = "beta", case_samples = cases, control_panel = controls)
res_beta$t5 <- epimutations(method = "beta", case_samples = cases, control_panel = controls, 
                              epi_params = params2)

sum(res_beta$ori$chromosome != 0)
sum(res_beta$t5$chromosome != 0)
## Seleccionar 1e-5
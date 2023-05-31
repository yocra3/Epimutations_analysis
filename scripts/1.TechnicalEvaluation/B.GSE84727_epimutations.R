#'#################################################################################
#'#################################################################################
#' Define epimutations in GSE84727
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(BiocParallel)
library(epimutacions)
library(readxl)
library(tidyverse)
library(ramr)

### Load data
load("data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")


methods <- names(epi_parameters())
names(methods) <- methods

### Run epimutation detection
res.GSE84727.list <- lapply(methods, epimutations_one_leave_out, methy = gset, 
                            BPPARAM = MulticoreParam(3))
save(res.GSE84727.list, file = "results/epimutations/GSE84727.epimutations.allSamples.Rdata")

### Run epimutation detection with ramr
## Convert GRSet to ramr format
data.ranges <- granges(gset)
data.betas  <- getBeta(gset)
sample.ids  <- colnames(gset)
mcols(data.ranges) <- data.betas

res.GSE84727.ramr.list <- lapply(c("IQR", "beta", "wbeta"), function(met){
  getAMR(data.ranges, ramr.method = met, min.cpgs = 3,
         qval.cutoff = ifelse(met == "beta", 1e-3, 1e-6), iqr.cutoff = 3,
         merge.window = 1000, cores = 10)
})
save(res.GSE84727.ramr.list, file = "results/epimutations/GSE84727.epimutations.allSamples.ramr.Rdata")

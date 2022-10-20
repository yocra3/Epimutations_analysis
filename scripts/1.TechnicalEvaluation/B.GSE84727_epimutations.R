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

### Load data
load("data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")


methods <- names(epi_parameters())
names(methods) <- methods

### Run epimutation detection
res.GSE84727.list <- lapply(methods, epimutations_one_leave_out, methy = gset, 
                            BPPARAM = MulticoreParam(3))
save(res.GSE84727.list, file = "results/epimutations/GSE84727.epimutations.allSamples.Rdata")
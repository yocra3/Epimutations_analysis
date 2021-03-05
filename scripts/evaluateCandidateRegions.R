#'#################################################################################
#'#################################################################################
#' Evaluate candidate regions defined based on 450K annotation
#' 
#' We will use two estimates
#' 1. Correlation between CpGs inside candidate regions vs other CpGs
#' 2. Overlap with literature epimutations
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)

load("results/candidateRegions/candidateRegions.450K.GRanges.Rdata")

# 1. CpGs correlation ####
load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Remove duplicates
gset.filt <- gset[, !gset$dup | gset$Batch == "Esteller"]

## Compute correlations
beta <- getBeta(gset.filt)

regsCorrs <- lapply(seq_len(length(candRegsGR)), function(x){
  set <- subsetByOverlaps(gset.filt, candRegsGR[x])
  cors <- cor(t(getBeta(set)))
  cors[upper.tri(cors)]
})
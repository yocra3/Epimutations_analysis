#'#################################################################################
#'#################################################################################
#' Add NAs to CpGs with high detection p-values
#'#################################################################################
#'#################################################################################

# Load libraries
library(minfi)

load("data/HELIX.detectionPvals.Rdata")
load("data/HELIX.GenomicRatioSet.Rdata")


gset <- methylome_subcohort_ComBatSlide_6cells

map <- colData(gset)[, c("arrayName", "SampleID")]
rownames(map) <- map$arrayName
colnames(detP) <- map[colnames(detP), "SampleID"]
dp.f <- detP[rownames(gset), colnames(gset)]

beta <- getBeta(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = "results/preprocess/HELIX/HELIX.withNA.GenomicRatioSet.Rdata")

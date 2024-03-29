#'#################################################################################
#'#################################################################################
#' Prepare GSE42861 dataset
#'#################################################################################
#'#################################################################################

# GSE42861 ####
library(GEOquery)
library(minfi)
library(meffil)
library(tidyverse)

options(timeout = 10000)

geo.full <- getGEO("GSE42861", GSEMatrix = TRUE)[[1]]



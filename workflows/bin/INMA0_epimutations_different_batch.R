#! Rscript
#'#################################################################################
#'#################################################################################
#' Run methylation QC data
#' This code performs quality control to methylation data. 
#' Input: 
#' - path to config file with the parameters for the QC
#' - Name of the dataset to name the output files
#' Important decisions:
#' - Remove samples based on QC
#' - Use values from meffil vignette in all parameters - see config file
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetPath <- args[1]
outPrefix <- args[2]

## Load libraries
library(minfi)
library(epimutacions)
library(tidyverse)

load(gsetPath)

## Select samples from Esteller and duplicates
gset.sel <- gset[, gset$Batch == "Esteller" | gset$dup]

### Make cases and controls datasets
cases <- gset.sel[, gset.sel$dup]
controls <- gset.sel[, !gset.sel$dup]

methods <- names(epi_parameters())
names(methods) <- methods
res <- lapply(methods, epimutations, case_samples = cases, control_panel = controls)
save(res, file = paste0(outPrefix, ".epimutations.INMA0.duplicates.Rdata"))


res.df <- Reduce(rbind, res) %>%
  mutate(method = rep(methods[methods != "barbosa"], sapply(res, nrow))) %>%
  left_join(colData(cases) %>% data.frame() %>% select(Sample_Name, idnum, batch) %>% mutate(sample = Sample_Name), by = "sample") %>%
  group_by(batch, idnum, method) %>%
  summarize(n = n())

lapply(unique(res.df$method), function(x) {
  
  mini <- subset(res.df, method == x) %>% arrange(idnum)
  
  t.test(subset(mini, batch == "Esteller")$n, subset(mini, batch == "MeDALL")$n, 
         paired = TRUE)
})

spread(res.df, batch, n) %>% 
  mutate(diff = MeDALL - Esteller) %>%
  ggplot(aes(x = method, y = diff)) +
  geom_boxplot()

ggplot(res.df, aes(x = batch, y = n, color = factor(idnum), group = factor(idnum))) +
  geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~ method, scales = "free_y")
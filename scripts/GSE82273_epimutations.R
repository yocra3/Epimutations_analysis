#'#################################################################################
#'#################################################################################
#' Define epimutations in GSE82273
#'#################################################################################
#'#################################################################################

library(minfi)
library(BiocParallel)
library(epimutacions)
library(readxl)
library(tidyverse)

load("data/GSE82273/GSE82273.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")


methods <- names(epi_parameters())
names(methods) <- methods

res.gse82273.list <- lapply(methods, epimutations_one_leave_out, methy = gset, 
                            BPPARAM = MulticoreParam(3))
save(res.gse82273.list, file = "results/epimutations/gse82273.epimutations.allSamples.Rdata")


gset$samp <- gsub("Blood of newborn infant, ", "", gset$title)

epi_lit <- read_excel("data/Epimutations.PMID32937144.xlsx", skip = 2)
epi_lit <- epi_lit[, 1:6]
colnames(epi_lit) <- c("chr", "start", "end", "cohort", "ids", "direction")

epi_geo <- subset(epi_lit, grepl("GSE82273", cohort)) %>%
  mutate(id_list = strsplit(ids, ","),
         dir_list = strsplit(direction, ",")) %>%
  unnest(c(id_list, dir_list)) %>%
  filter(id_list %in% gset$samp)


epi_geo_GR <- makeGRangesFromDataFrame(epi_geo, keep.extra.columns = TRUE)
overlap <- lapply(methods[-5], function(x) {
  print(x)
  ## Remove samples without epimutations
  tab <- subset(res.gse82273.list[[x]], start != 0)
  tab$id <- colData(gset)[tab$sample, "samp"]
  if (x %in% c("mlm", "manova")){
    tab <- subset(tab, is.na(pvalue) | pvalue < 0.05/40408)
  }
  tab_GR <- makeGRangesFromDataFrame(tab, keep.extra.columns = TRUE)
  
  sapply(seq_len(length(epi_geo_GR)), function(i){
    gr <- epi_geo_GR[i]
    tab_GRm <- tab_GR[tab_GR$id == gr$id_list]
    length(findOverlaps(gr, tab_GRm)) > 0
    
  })
  
})


overlap_samps <- lapply(methods[-5], function(x) {
  print(x)
  ## Remove samples without epimutations
  tab <- subset(res.gse82273.list[[x]], start != 0)
  tab$id <- colData(gset)[tab$sample, "samp"]
  if (x %in% c("mlm", "manova")){
    tab <- subset(tab, is.na(pvalue) | pvalue < 0.05/40408)
  }
  subset(tab, id %in% epi_geo$id_list)
})
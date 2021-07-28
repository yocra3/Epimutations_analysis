#'#################################################################################
#'#################################################################################
#' Define epimutations in GSE84727
#'#################################################################################
#'#################################################################################

library(minfi)
library(BiocParallel)
library(epimutacions)
library(readxl)
library(tidyverse)

load("data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")


methods <- names(epi_parameters())
names(methods) <- methods

res.GSE84727.list <- lapply(methods, epimutations_one_leave_out, methy = gset, 
                            BPPARAM = MulticoreParam(3))
save(res.GSE84727.list, file = "results/epimutations/GSE84727.epimutations.allSamples.Rdata")


gset$samp <- gset$`sentrixids:ch1`

epi_lit <- read_excel("data/Epimutations.PMID32937144.xlsx", skip = 2)
epi_lit <- epi_lit[, 1:6]
colnames(epi_lit) <- c("chr", "start", "end", "cohort", "ids", "direction")

epi_geo <- subset(epi_lit, grepl("GSE84727", cohort)) %>%
  mutate(id_list = strsplit(ids, ","),
         dir_list = strsplit(direction, ",")) %>%
  unnest(c(id_list, dir_list)) %>%
  filter(id_list %in% gset$samp)


epi_geo_GR <- makeGRangesFromDataFrame(epi_geo, keep.extra.columns = TRUE)
epi_geo_GR.filt <- subsetByOverlaps(epi_geo_GR, rowRanges(gset))
overlap <- lapply(methods, function(x) {
  print(x)
  ## Remove samples without epimutations
  tab <- subset(res.GSE84727.list[[x]], start != 0)
  tab$id <- colData(gset)[tab$sample, "samp"]
  if (x %in% c("mlm", "manova")){
    tab <- subset(tab, is.na(pvalue) | pvalue < 0.05/40408)
  }
  tab_GR <- makeGRangesFromDataFrame(tab, keep.extra.columns = TRUE)
  
  sapply(seq_len(length(epi_geo_GR.filt)), function(i){
    gr <- epi_geo_GR.filt[i]
    tab_GRm <- tab_GR[tab_GR$id == gr$id_list]
    length(findOverlaps(gr, tab_GRm)) > 0
    
  })
  
})


overlap_method <- lapply(methods, function(x) {
  print(x)
  ## Remove samples without epimutations
  tab <- subset(res.GSE84727.list[[x]], start != 0)
  tab$id <- colData(gset)[tab$sample, "samp"]
  if (x %in% c("mlm", "manova")){
    tab <- subset(tab, is.na(pvalue) | pvalue < 0.05/40408)
  }
  tab_GR <- makeGRangesFromDataFrame(tab, keep.extra.columns = TRUE)
  
  sapply(seq_len(length(tab_GR)), function(i){
    gr <- tab_GR[i]
    geo_grm <- epi_geo_GR[epi_geo_GR$id_list == gr$id]
    length(findOverlaps(gr, geo_grm)) > 0
    
  })
})

samps <- lapply(res.GSE84727.list, function(x) unique(subset(x, start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))$sample))



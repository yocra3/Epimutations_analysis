library(minfi)
library(meffil)
library(tidyverse)
library(cowplot)
library(epimutacions)

load("../INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
INMA_ind <- gset
load("../INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
INMA_comb <- gset
load("../INMA0combined.epimutations.INMA0.duplicates.Rdata")
res_ind <- res
load("../INMA_comb.epimutations.INMA0.duplicates.Rdata")
res_comb <- res
ind.Esteller <- INMA_ind[, INMA_ind$Batch == "Esteller" | INMA_ind$dup]
comb.gset.Esteller <- INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup]
nMethod <- sapply(res_ind, nrow)
nMethod[sapply(nMethod, is.null)] <- 0
ind.ids <- ind.Esteller$idnum[duplicated(ind.Esteller$idnum)]
ind.res.df <- Reduce(rbind, res_ind) %>%
  mutate(method = rep(names(res_ind), unlist(nMethod))) %>%
  full_join(colData(ind.Esteller) %>% 
              data.frame() %>% 
              filter(idnum %in% ind.ids) %>% 
              select(Sample_Name, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent")
comb.ids <- comb.gset.Esteller$idnum[duplicated(comb.gset.Esteller$idnum)]
comb.res.df <- Reduce(rbind, res_comb) %>%
  mutate(method = rep(names(res_comb), sapply(res_comb, nrow))) %>%
  full_join(colData(comb.gset.Esteller) %>% 
              data.frame() %>% 
              filter(idnum %in% comb.ids) %>% 
              select(Sample_Name, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Joint")

all.res.df <- rbind(ind.res.df, comb.res.df) %>%
  filter(method != "barbosa")

ids <- unique(all.res.df$idnum)
methods <- unique(all.res.df$method)
norms <- unique(all.res.df$Normalization)
combids <- data.frame(expand_grid(ids, methods))
esteller.df <- subset(all.res.df, Batch == "Esteller")
overlapsEsteller <- sapply(seq_len(nrow(combids)), function(x) {
  sub <- subset(esteller.df, idnum == combids[x, 1] & method == combids[x, 2] &
                  (pvalue < 0.05/40408 | is.na(pvalue))  )
  
  if (sum(sub$Normalization == "Independent") == 0 | sum(sub$Normalization == "Joint") == 0){
    return(c(common = 0, independent = 0, joint = 0))
  }
  if (nrow(sub) == 0){
    return(c(common = 0, independent = 0, joint = 0))
  }
  
  indep <- makeGRangesFromDataFrame(subset(sub, Normalization == "Independent"))
  joint <- makeGRangesFromDataFrame(subset(sub, Normalization == "Joint"))
  com <- sum(countOverlaps(indep, joint))
  c(common = com, independent = length(indep) - com, joint = length(joint) - com)
})
overlapsEsteller.df <- overlapsEsteller %>%
  t() %>%
  cbind(combids) %>%
  tibble() %>%
  filter(!(common == 0 & independent == 0 & joint == 0)) %>%
  filter(!(ids %in% c("204", "484", "655", "339")))


## Epimutations samp 16 - mlm/manova
sub <- subset(esteller.df, idnum == "16" & method %in% c("mlm", "manova") &
                (pvalue < 0.05/40408 | is.na(pvalue)) )

sub[, c("chromosome", "start", "end", "method", "pvalue", "Normalization")] %>% 
  arrange(chromosome, start)

sub.uniq <- distinct(sub, chromosome, start, end, .keep_all = TRUE)

ind.set <- ind.Esteller[, (!ind.Esteller$dup | ind.Esteller$idnum == "16") & ind.Esteller$Batch == "Esteller"]
comb.set <- comb.gset.Esteller[, (!comb.gset.Esteller$dup | comb.gset.Esteller$idnum == "16") & comb.gset.Esteller$Batch == "Esteller"]

## Edit plot_epimutations to only show methylation values
plot_list <- lapply(seq_len(nrow(sub.uniq)), function(i){
  p1 <- plot_epimutations(sub.uniq[i,], ind.set) + ggtitle(paste(sub.uniq[i, c("chromosome", "start")], "Independent"))
  p2 <- plot_epimutations(sub.uniq[i,], comb.set) + ggtitle(paste(sub.uniq[i, c("chromosome", "start")], "Combined"))
  list(p1, p2)
  })

pdf("Epimutations_id16.pdf")
lapply(plot_list, function(x) plot_grid(plotlist = x, nrow = 2))
dev.off()


#' Overview of results ####
#' Reg 1 (chr10-132.91Mb)
#'  - Detected in both datasets
#'  - Detected only with manova
#'  - likely deletion - a lot of control samples have similar values
#' Reg 2 (chr10-132.90Mb) 
#'  - Detected in both datasets
#'  - Detected only with manova
#'  - likely deletion - a lot of control samples have similar values
#' Reg 3 (chr10-101.87Mb) 
#'  - Detected in both datasets
#'  - Detected with manova and mlm
#'  - Polymorphic? - Other control samples have similar values
#' Reg 4 (chr15-101.87Mb) 
#'  - Detected only with independent normalization - more outlier values
#'  - Detected with manova and mlm
#'  - Another control samples have similar values
#' Reg 5 (chr16-0.63Mb) 
#'  - Detected only with joint normalization - more outlier values
#'  - Detected with manova and mlm
#'  - No samples with similar values
#' Reg 6 (chr16-54.10Mb) 
#'  - Detected only with joint normalization - more outlier values
#'  - Detected with manova
#'  - False positive, no epimutation is observed in the plot
#' Reg 7 (chr2-5.81Mb) 
#'  - Detected only with joint normalization - more outlier values
#'  - Detected with manova
#'  - False positive, no epimutation is observed in the plot
#'#################################################################################
#'#################################################################################
#' Figures for poster
#'#################################################################################
#'#################################################################################

# Load data and libraries ####
library(minfi)
library(meffil)
library(tidyverse)
library(cowplot)
library(ggh4x)
library(readxl)


getMeanQuantile <- function(cpgs, sampid, set){
  idnum <- colData(set)[sampid, "idnum"]
  samps <- colnames(set[, set$idnum == idnum])
  samp2 <- samps[samps != sampid]
  betas <- getBeta(set[cpgs, ])
  quant <- apply(betas, 1, function(x) {
    f <- ecdf(x)
    f(x[colnames(betas) == samp2])
  })
  
  mean(quant)
}

# Normalization PCAs ####
## Independent normalization ####
load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
INMA_ind <- gset

ind.pcs <- meffil.methylation.pcs(getBeta(INMA_ind), probe.range = 40000, full.obj = TRUE)
ind.pcs.df <- ind.pcs$x %>%
  data.frame() %>%
  select(PC1, PC2) %>%
  mutate(Sample_Name = rownames(.)) %>%
  left_join(colData(INMA_ind) %>% data.frame() %>% select(Sample_Name, Sex, dup, Batch, idnum), by = "Sample_Name") %>%
  mutate(Batch2 = ifelse(Batch == "Esteller", "Reference", "Alternative")) %>%
  as_tibble()

idnum.tab.indep <- ind.pcs.df %>%
  filter(dup) %>%
  group_by(idnum) %>%
  summarize(n = length(unique(Batch))) %>%
  mutate(type = ifelse(n == 1, "Replicate same batch", "Replicate different batch")) %>%
  data.frame()
rownames(idnum.tab.indep) <- idnum.tab.indep$idnum  
ind.pcs.df$type <- factor(ifelse(ind.pcs.df$dup, idnum.tab.indep[as.character(ind.pcs.df$idnum), "type"], 
                                 ind.pcs.df$Batch2),
                          levels = c("Reference", "Alternative", "Replicate same batch", "Replicate different batch"))

## Joint normalization ####
load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
INMA_comb <- gset

joint.pcs <- meffil.methylation.pcs(getBeta(INMA_comb), probe.range = 40000, full.obj = TRUE)
joint.pcs.df <- joint.pcs$x %>%
  data.frame() %>%
  select(PC1, PC2) %>%
  mutate(Sample_Name = rownames(.)) %>%
  left_join(colData(INMA_comb) %>% data.frame() %>% select(Sample_Name, Sex, dup, Batch, idnum), by = "Sample_Name") %>%
  mutate(Batch2 = ifelse(Batch == "Esteller", "Reference", "Alternative")) %>%
  as_tibble()

idnum.tab.joint <- joint.pcs.df %>%
  filter(dup) %>%
  group_by(idnum) %>%
  summarize(n = length(unique(Batch))) %>%
  mutate(type = ifelse(n == 1, "Replicate same batch", "Replicate different batch")) %>%
  data.frame()
rownames(idnum.tab.joint) <- idnum.tab.joint$idnum  
joint.pcs.df$type <- factor(ifelse(joint.pcs.df$dup, idnum.tab.joint[as.character(joint.pcs.df$idnum), "type"], 
                                   joint.pcs.df$Batch2),
                            levels = c("Reference", "Alternative", "Replicate same batch", "Replicate different batch"))

# Merge epimutations datasets ####
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.Rdata")

ind.res.df <- Reduce(rbind, res_indep) %>%
  mutate(method = rep(names(res_indep), sapply(res_indep, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))


comb.res.df <- Reduce(rbind, res_joint) %>%
  mutate(method = rep(names(res_joint), sapply(res_joint, nrow))) %>%
  left_join(joint.pcs.df %>% 
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))

comSamps <- intersect(unique(ind.res.df$sample), unique(comb.res.df$sample))

all.res.df <- rbind(ind.res.df, comb.res.df) %>%
  filter(sample %in% comSamps)


## Technical replicates  ####
tech.rep.res <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_ind[, INMA_ind$Batch == "Esteller"]))) %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Significant and top", "Replicate-specific epimutation"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both significant", "Significant and top", "Replicate-specific epimutation")), 
         epi_id = paste(idnum, epi_region_id)) 
tech.rep.plot <- tech.rep.res %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  ggtitle("Technical replicates - Same batch") +
  facet_grid(idnum ~ method, scales = "free") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

png("figures/poster/INMA0.techRep.png")
tech.rep.plot
dev.off()

## Batch replicates ####
batch.rep.res <- comb.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup]))) %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Significant and top", "Replicate-specific epimutation"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both significant", "Significant and top", "Replicate-specific epimutation")), 
         epi_id = paste(idnum, epi_region_id)) 


batch.rep.plot <- batch.rep.res %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  ggtitle("Technical replicates - Different batch") +
  facet_grid(idnum ~ method, scales = "free") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

png("figures/poster/INMA0.BatchRep.png")
batch.rep.plot
dev.off()
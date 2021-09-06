#'#################################################################################
#'#################################################################################
#' Figures from INMA0 normalization 
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
ind.pcs.vars <- ind.pcs$sdev^2/sum(ind.pcs$sdev^2)
indep_pc <- ggplot(ind.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  theme_bw() +
  ggtitle("Independent Normalization") +
  scale_x_continuous(name = paste0("PC1 (", round(ind.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(ind.pcs.vars[2]*100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")
summary(lm(PC1 ~ Batch, ind.pcs.df))

# png("figures/INMA0.indepNorm.PCA.png")
# indep_pc
# dev.off()



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
joint.pcs.vars <- joint.pcs$sdev^2/sum(joint.pcs$sdev^2)
joint_pc <- ggplot(joint.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  scale_x_continuous(name = paste0("PC1 (", round(joint.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(joint.pcs.vars[2]*100, 1), "%)")) +
  theme_bw() +
  ggtitle("Joint Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

summary(lm(PC1 ~ Batch, joint.pcs.df))

# png("figures/INMA0.jointNorm.PCA.png")
# joint_pc
# dev.off()
# 

png("figures/INMA0.PCA_panel.png", width = 600, height = 300)
plot_grid(indep_pc, joint_pc, labels = "AUTO", nrow = 1, rel_widths = c(5, 8))
dev.off()



# Merge epimutations datasets ####
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.liberal.Rdata")
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.liberal.Rdata")

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


ind.res.lib.df <- Reduce(rbind, res_indep_lib) %>%
  mutate(method = rep(names(res_indep_lib), sapply(res_indep_lib, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))


# Technical replicates in Esteller ####
## Total epimutations ####
ind.res.df %>%
  filter(type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(n = sum(chromosome != 0)) %>%
  filter(!is.na(epi_region_id)) %>%
  group_by(method, idnum) %>%
  summarize(N = length(unique(epi_region_id))) %>%
  spread(idnum, N)


ind.lib.replicates <- ind.res.lib.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(method, idnum, epi_region_id))

#'###############################################################################
#'### Borrar? ####
# tech.num.plot <- ind.res.df %>%
#   filter(chromosome != 0 & type == "Replicate same batch") %>%
#   group_by(method, idnum, epi_region_id) %>%
#   mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
#          epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
#   ungroup() %>%
#   count(method, epi_type) %>% 
#   complete(method, epi_type, fill = list(n = 0)) %>%
#   ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected") +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "")
# 
# png("figures/INMA0.techRep.num.png")
# tech.num.plot
# dev.off()
# 
# tech.num.lib.plot <- ind.res.df %>%
#   filter(chromosome != 0 & type == "Replicate same batch") %>%
#   group_by(method, idnum, epi_region_id) %>%
#   summarize(epi_type = ifelse(length(Batch) == 2, "Both significant", "Replicate-specific epimutation")) %>%
#   mutate(epi_name = paste(method, idnum, epi_region_id), 
#          epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% ind.lib.replicates$epi_name, "Significant and suggestive", epi_type),
#          epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation")),
#          method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist")) %>%
#   ungroup() %>%
#   count(method, epi_type) %>% 
#   complete(method, epi_type, fill = list(n = 0)) %>%
#   ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected") +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "")
# 
# png("figures/INMA0.techRep.num.lib.png")
# tech.num.lib.plot
# dev.off()
#'###############################################################################


## Epimutations per replicate ####
tech.ind.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  group_by(idnum) %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ .)

png("figures/INMA0.techRep.ind.png", height = 300)
tech.ind.plot
dev.off()



tech.ind.lib.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% ind.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(method, idnum, epi_type) %>% 
  complete(method, idnum,  epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ .)


png("figures/INMA0.techRep.ind.lib.png", height = 300)
tech.ind.lib.plot
dev.off()



tech.ind.top.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_ind[, INMA_ind$Batch == "Esteller"]))) %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Significant and top", "Replicate-specific epimutation"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both significant", "Significant and top", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(method, idnum, epi_type) %>% 
  complete(method, idnum,  epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ .)

png("figures/INMA0.techRep.ind.top.png", height = 300)
tech.ind.top.plot
dev.off()





## CpGs replicability ####
rel <- read_xlsx("data/PMID32885222_reliability.xlsx", skip = 2) %>%
  mutate(cpg = `Illumina Probe ID`) %>%
  dplyr::select(cpg, Reliability) 

tech.ind.cpg <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation"),
        epi_name = paste(method, idnum, epi_region_id), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% ind.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  mutate(cpg = strsplit(cpg_ids, ",")) %>%
  unnest(cpg) %>%
  select(method, idnum, epi_region_id, cpg, epi_type) %>%
  distinct() %>%
  left_join(rel, by = "cpg")

tech.ind.cpg.plot <- tech.ind.cpg %>%
  select(method, cpg, epi_type, Reliability) %>%
  distinct() %>%
  ggplot(aes(x = epi_type, y = Reliability)) +
  geom_boxplot() +
  facet_wrap(~ method) +
  theme_bw()

png("figures/INMA0.techRep.cpg.rel.png")
tech.ind.cpg.plot
dev.off()

tech.ind.cpg.relprop.plot <- tech.ind.cpg %>%
  group_by(method, idnum, epi_region_id, epi_type) %>%
  summarize(p = mean(Reliability > 0.4)) %>%
  ggplot(aes(x = epi_type, y = p)) +
  geom_boxplot() +
  facet_wrap(~ method) +
  theme_bw()

png("figures/INMA0.techRep.cpg.proprel.png")
tech.ind.cpg.relprop.plot
dev.off()

tech.ind.cpg.relnum.plot <- tech.ind.cpg %>%
  group_by(method, idnum, epi_region_id, epi_type) %>%
  summarize(n = sum(Reliability > 0.4)) %>%
  ggplot(aes(x = epi_type, y = n)) +
  geom_boxplot() +
  facet_wrap(~ method) +
  theme_bw()

png("figures/INMA0.techRep.cpg.numrel.png")
tech.ind.cpg.relnum.plot
dev.off()

beta.epimut.same <- subset(ind.res.df, method == "beta" & chromosome != 0 & type == "Replicate same batch")
pdf("figures/epimut.beta.techReps.pdf")
lapply(seq_len(nrow(beta.epimut.same)), function(i){
  row <- beta.epimut.same[i, ]
  p1 <- plot_epimutations(row, INMA_ind)
  samp <- row$sample
  idnum <- colData(INMA_ind)[samp, "idnum"]
  samps <- colnames(INMA_ind[, INMA_ind$idnum == idnum])
  samp2 <- samps[samps != samp]
  row$sample <- samp2
  p2 <- plot_epimutations(row, INMA_ind)
  plot_grid(p1, p2)
})
dev.off()


manova.epimut.same <- subset(ind.res.df, method == "manova" & chromosome != 0 & type == "Replicate same batch")
pdf("figures/epimut.manova.techReps.pdf")
lapply(seq_len(nrow(manova.epimut.same)), function(i){
  row <- manova.epimut.same[i, ]
  p1 <- plot_epimutations(row, INMA_ind)
  samp <- row$sample
  idnum <- colData(INMA_ind)[samp, "idnum"]
  samps <- colnames(INMA_ind[, INMA_ind$idnum == idnum])
  samp2 <- samps[samps != samp]
  row$sample <- samp2
  p2 <- plot_epimutations(row, INMA_ind)
  plot_grid(p1, p2)
})
dev.off()


# png("figures/INMA0.techRep.panel.png", width = 1200, height = 500)
# plot_grid(tech.num.plot, tech.ind.plot, labels = "AUTO")
# dev.off()

# Diff. norm. - tech rep. ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.Rdata")
load("results/epimutations/Esteller.epimutations.normalization.liberal.Rdata")
norm.res.df <- Reduce(rbind, 
                      lapply(INMA_norm, function(x) {Reduce(rbind, x) %>%
                                      mutate(method = rep(names(x), sapply(x, nrow)))}
                             )) %>%
  mutate(Normalization = rep(names(INMA_norm), 
                             sapply(INMA_norm, function(x) sum(sapply(x, nrow)))),
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist")) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(Sample_Name %in% unique(ind.res.df$sample))
norm.res.df <- rbind(norm.res.df, 
   mutate(ind.res.df, Normalization = "meffil") %>% filter(sample %in% norm.res.df$sample)) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization"))) 


norm.res.lib.df <- Reduce(rbind, 
                      lapply(INMA_norm_lib, function(x) {Reduce(rbind, x) %>%
                          mutate(method = rep(names(x), sapply(x, nrow)))}
                      )) %>%
  mutate(Normalization = rep(names(INMA_norm_lib), 
                             sapply(INMA_norm_lib, function(x) sum(sapply(x, nrow)))),
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist")) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(Sample_Name %in% unique(ind.res.df$sample))
norm.res.lib.df <- rbind(norm.res.lib.df, 
                         mutate(ind.res.lib.df, Normalization = "meffil") %>% filter(sample %in% ind.res.lib.df$sample)) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization"))) 

#'##############################################################################
#' Borrar ?
#' ## Total number of epimutations per normalization algorithm ####
# tech.totnum.norm.plot <- norm.res.df %>%
#   filter(type == "Replicate same batch") %>%
#   group_by(method, Normalization) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
#   geom_bar(stat = "identity") + theme_bw() +
#   facet_wrap(~ method)
# png("figures/INMA0.techRep.Norm.totnum.png")
# tech.totnum.norm.plot
# dev.off()
#'##############################################################################

## Technical replicates per normalization algorithm ####
tech.num.norm.plot <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  mutate(epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(y = n + 1, label = n)) +
  theme_bw() +
  facet_grid(method ~ Normalization) +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


png("figures/INMA0.techRep.Norm.num.png", width = 1000, height = 700)
tech.num.norm.plot
dev.off()

norm.res.lib.replicates <- norm.res.lib.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization))


tech.num.norm.lib.plot <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(y = n + 1, label = n)) +
  theme_bw() +
  facet_grid(method ~ Normalization) +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png("figures/INMA0.techRep.Norm.num.lib.png", width = 1000, height = 700)
tech.num.norm.lib.plot
dev.off()


tech.num.norm.epi.plot <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation")),
         epi_id = paste(idnum, epi_region_id)) %>%
 ggplot(aes(x = Normalization, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(~ method) + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))
  

png("figures/INMA0.techRep.Norm.num.epi.png", width = 1000, height = 700)
tech.num.norm.epi.plot
dev.off()


norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  group_by(Normalization,method) %>%
  summarize(n = n()) %>%
  spread(Normalization, n) %>%
  data.frame()

norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  group_by(Normalization,method) %>%
  summarize(p = mean(epi_type == "Both significant")) %>%
  spread(Normalization, p) %>%
  data.frame()

norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  group_by(Normalization,method) %>%
  summarize(p = mean(epi_type != "Replicate-specific epimutation")) %>%
  spread(Normalization, p) %>%
  data.frame()


# png("figures/INMA0.techRep.Norm.panel.png", width = 1200, height = 500)
# plot_grid(tech.totnum.norm.plot, tech.num.norm.plot, labels = "AUTO")
# dev.off()
# 
#'##############################################################################
#' Borrar ?
#'# Total number of epimutations per normalization algorithm in all Esteller samples ####
# tech.totnum.all.norm.plot <- norm.res.df %>%
#   filter(Batch == "Esteller") %>%
#   group_by(method, Normalization) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
#   geom_bar(stat = "identity") + theme_bw() +
#   facet_wrap(~ method)
# png("figures/INMA0.techRep.Norm.totnum.allsamps.png")
# tech.totnum.all.norm.plot
# dev.off()
#'##############################################################################


# Diff. norm. same samp. ####
tech.samesamp.num.norm.plot <- norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Default parameters")

png("figures/INMA0.sameSamp.Norm.propOverlap.png", height = 300)
tech.samesamp.num.norm.plot
dev.off()

norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  select(-n) %>% spread(n_norms, prop)


tech.samesamp.num.norm.plot.lib <- norm.res.lib.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  mutate(epi_id = paste(method, Sample_Name, epi_region_id)) %>%
  filter(epi_id %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id)) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Liberal parameters")

png("figures/INMA0.sameSamp.Norm.propOverlap2.png", height = 300)
tech.samesamp.num.norm.plot.lib
dev.off()


tech.samesamp.norm.epi.plot <- norm.res.lib.df %>%
  filter(chromosome != 0) %>%
  mutate(epi_name = paste(method, Sample_Name, epi_region_id),
         epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
  filter(epi_name %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id)) %>%
  mutate(epi_cat = ifelse(epi_name2 %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id, norm.res.df$Normalization),
                          "Significant", "Suggestive"), 
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_grid(~ method) + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.sameSamp.Norm.num.epi.png", width = 1000, height = 700)
tech.samesamp.norm.epi.plot
dev.off()

# 
# a <- norm.res.lib.df %>%
#   filter(chromosome != 0) %>%
#   mutate(epi_name = paste(method, Sample_Name, epi_region_id),
#          epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
#   filter(epi_name %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id)) %>%
#   mutate(epi_cat = ifelse(epi_name2 %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id, norm.res.df$Normalization),
#                           "Significant", "Suggestive"), 
#          epi_id = paste(Sample_Name, epi_region_id),
#          epi_meth = paste(method, Normalization)) %>%
#   count(epi_id, method, Normalization) %>%
#   complete(epi_id, method, Normalization, fill = list(n = 0)) %>%
#   spread(Normalization, n)
# a_mat <- data.matrix(a[, -c(1:2)])
# rownames(a_mat) <- a$epi_id
# 
# png("figures/cot.png", width = 1200, height = 500)
# heatmap(t(a_mat))
# dev.off()
# 
# png("figures/INMA0.sameSamp.Norm.panel.png", width = 1200, height = 500)
# plot_grid(tech.totnum.all.norm.plot, tech.samesamp.num.norm.plot2, labels = "AUTO")
# dev.off()


# Diff batch tech replicates ####
#'###############################################################################
#'### Borrar? ####
#'# Total epimutations ####
#' batch.totnum.norm.plot <- all.res.df %>%
#'   filter(type == "Replicate different Batch") %>%
#'   group_by(method, Normalization) %>%
#'   summarize(n = sum(chromosome != 0)) %>%
#'   ggplot(aes(x = method, y = n, fill = Normalization)) +
#'   geom_bar(stat = "identity", position = "dodge") + theme_bw() +
#'   scale_y_continuous(name = "Total epimutations detected")
#'   
#' 
#' png("figures/INMA0.Batch.totnum.png")
#' batch.totnum.norm.plot
#' dev.off()
#' batch.ind.plot <- all.res.df %>%
# filter(chromosome != 0 & type == "Replicate different Batch") %>%
#   group_by(Normalization, method, idnum, epi_region_id) %>%
#   mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
#          epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
#                            ifelse(grepl("04", epi_type), "Alternative", "Reference")),
#          epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
#   group_by(Normalization) %>%
#   count(method, epi_type) %>% 
#   complete(method, epi_type, fill = list(n = 0)) %>%
#   ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected") +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "") +
#   facet_grid(Normalization ~ ., scales = "free")
# 
# png("figures/INMA0.batch.overlap.png")
# batch.ind.plot
# dev.off()

#' #'###############################################################################

## Prepare data ####
comb.res.lib.df <- Reduce(rbind, res_joint_lib) %>%
  mutate(method = rep(names(res_joint_lib), sapply(res_joint_lib, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))

all.res.lib <- rbind(comb.res.lib.df, ind.res.lib.df) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name)) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(Normalization, method, idnum, epi_region_id))


## Epimutations per replicate ####
batch.ind.boxplot <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name)) %>%
  mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 55)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "blue", "red")) +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.batch.overlap.boxplot.png")
batch.ind.boxplot
dev.off()


batch.ind.boxplot.lib <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name)) %>%
  mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(paste(Normalization, method, idnum, epi_region_id) %in% all.res.lib$epi_name, 
                                  "Significant and suggestive", 
                                  ifelse(grepl("04", epi_type), "Alternative", "Reference"))),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Significant and suggestive", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 55)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "purple1", "blue", "red")) +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.batch.overlap.boxplot.lib.png")
batch.ind.boxplot.lib
dev.off()







set_list <- list(Joint = INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup], 
                 Independent = INMA_ind[, INMA_ind$Batch == "Esteller" | INMA_ind$dup])
batch.ind.boxplot.top <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, set_list[[.[i, ]$Normalization]]))) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name), 
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(rep_quant > 0.95 | rep_quant < 0.05, "Significant and top", 
                                  ifelse(grepl("04", epi_type), "Alternative", "Reference"))),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Significant and top", "Reference", "Alternative"))) %>%
  select(Normalization, method, idnum, epi_region_id, epi_type) %>%
  distinct() %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 55)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "purple1", "blue", "red")) +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.batch.overlap.boxplot.top.png")
batch.ind.boxplot.top
dev.off()


batch.ind.epi.plot <- rbind(comb.res.lib.df, ind.res.lib.df) %>%
  filter(chromosome != 0 & type == "Replicate different batch")%>%
  mutate(epi_name = paste(method, idnum, epi_region_id),
         epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
  filter(epi_name %in% paste(all.res.df$method, all.res.df$idnum, all.res.df$epi_region_id)) %>%
  mutate(epi_cat = ifelse(epi_name2 %in% paste(all.res.df$method, all.res.df$Sample_Name, all.res.df$epi_region_id, all.res.df$Normalization),
                          "Significant", "Suggestive"), 
         epi_id = paste(idnum, epi_region_id),
         Batch2 = recode(Batch, "Esteller" = "Ref.", "MeDALL" = "Alt."),
         Batch2 = factor(Batch2, levels = c("Ref.", "Alt.")),
         norm = recode(Normalization, "Independent" = "Indep."),
         epi_norm = paste(norm, "-", Batch2)) %>%
  ggplot(aes(x = Batch2, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_nested(idnum ~ method + norm, scales = "free_y", space="free") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
png("figures/INMA0.batch.overlap.epi.png", width = 600, height = 800)
batch.ind.epi.plot
dev.off()




# png("figures/INMA0.batch.panel.png", width = 1200, height = 500)
# plot_grid(batch.totnum.norm.plot, batch.ind.plot, labels = "AUTO")
# dev.off()
# 
# png("figures/INMA0.batch.panel2.png", width = 1200, height = 500)
# plot_grid(batch.totnum.norm.plot, batch.ind.boxplot, labels = "AUTO")
# dev.off()

## CpGs replicability ####
tech.batch.cpg <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name),
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(paste(Normalization, method, idnum, epi_region_id) %in% all.res.lib$epi_name, 
                                  "Significant and suggestive", 
                                  ifelse(grepl("04", epi_type), "Alternative", "Reference"))),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Significant and suggestive", "Reference", "Alternative"))) %>%
  ungroup() %>%
  mutate(cpg = strsplit(cpg_ids, ",")) %>%
  unnest(cpg) %>%
  select(method, Normalization, idnum, epi_region_id, cpg, epi_type) %>%
  distinct() %>%
  left_join(rel, by = "cpg")

b1_samps <- colnames(INMA_comb[, INMA_comb$Batch == "Esteller" & INMA_comb$idnum %in% tech.batch.cpg$idnum])
b2_samps <- colnames(INMA_comb[, INMA_comb$Batch == "MeDALL" & INMA_comb$idnum %in% tech.batch.cpg$idnum])
beta_diff <- lapply(b1_samps, function(b1) Reduce(rbind, lapply(b2_samps, function(b2) { 
  vec <- getBeta(INMA_comb[, b1]) - getBeta(INMA_comb[, b2])
  tibble(cpg = rownames(vec), diff = as.vector(vec), samp1 = b1, samp2 = b2)
  })))
beta_diff <- Reduce(rbind, beta_diff)


png("figures/diff.beta.techReps.png", width = 600, height = 800)
beta_diff %>%
  ggplot(aes(x = diff)) +
  geom_density() +
  facet_grid(samp1 ~ samp2) +
  theme_bw()
dev.off()

png("figures/diff.beta.techReps2.png")
beta_diff %>% 
  group_by(samp1, samp2) %>% 
  summarize(p = mean(abs(diff) < 0.05, na.rm = TRUE)) %>%
  ggplot(aes(x = samp2, y = p, color = samp2)) +
  geom_point() +
  facet_wrap(~ samp1) +
theme_bw()
dev.off()

png("figures/diff.beta.techReps3.png", width = 600, height = 800)
beta_diff %>%
  filter(cpg %in% tech.batch.cpg$cpg) %>%
  ggplot(aes(x = diff)) +
  geom_density() +
  facet_grid(samp1 ~ samp2) +
  theme_bw()
dev.off()

png("figures/diff.beta.techReps4.png")
beta_diff %>% 
  filter(cpg %in% tech.batch.cpg$cpg) %>%
  group_by(samp1, samp2) %>% 
  summarize(p = mean(abs(diff) < 0.05, na.rm = TRUE)) %>%
  ggplot(aes(x = samp2, y = p, color = samp2)) +
  geom_point() +
  facet_wrap(~ samp1) +
  theme_bw()
dev.off()


png("figures/diff.beta.techReps5.png", width = 800, height = 800)
beta_diff %>% 
  left_join(colData(INMA_comb)[, c("Sample_Name", "idnum")] %>% data.frame() %>% mutate(samp1 = Sample_Name), by = "samp1") %>%
  mutate(cpg_id = paste(cpg, idnum),
         cpg_type = ifelse(cpg_id %in% paste(tech.batch.cpg$cpg, tech.batch.cpg$idnum),
                           "Individual epimutation", 
                           ifelse(cpg %in% tech.batch.cpg$cpg, "Any Epimutation", "No epimutation"))) %>%
  group_by(cpg_type, samp1, samp2) %>% 
  summarize(p = mean(abs(diff) < 0.05, na.rm = TRUE)) %>%
  ggplot(aes(x = cpg_type, y = p, color = cpg_type)) +
  geom_point() +
  facet_grid(samp2 ~ samp1) +
  theme_bw()
dev.off()

png("figures/diff.beta.techReps6.png", width = 800, height = 800)
beta_diff %>% 
  left_join(colData(INMA_comb)[, c("Sample_Name", "idnum")] %>% data.frame() %>% mutate(samp1 = Sample_Name), by = "samp1") %>%
  mutate(cpg_id = paste(cpg, idnum),
         cpg_type = ifelse(cpg_id %in% paste(tech.batch.cpg$cpg, tech.batch.cpg$idnum),
                           "Individual epimutation", 
                           ifelse(cpg %in% tech.batch.cpg$cpg, "Any Epimutation", "No epimutation"))) %>%
  ggplot(aes(x = cpg_type, y = diff, color = cpg_type)) +
  geom_boxplot() +
  facet_grid(samp2 ~ samp1) +
  theme_bw()
dev.off()

png("figures/diff.beta.techReps7.png", width = 800, height = 800)
beta_diff %>% 
  left_join(colData(INMA_comb)[, c("Sample_Name", "idnum")] %>% data.frame() %>% mutate(samp1 = Sample_Name), by = "samp1") %>%
  mutate(cpg_id = paste(cpg, idnum),
         cpg_type = ifelse(cpg_id %in% paste(tech.batch.cpg$cpg, tech.batch.cpg$idnum),
                           "Individual epimutation", 
                           ifelse(cpg %in% tech.batch.cpg$cpg, "Any Epimutation", "No epimutation"))) %>%
  group_by(cpg_type, samp1, samp2) %>% 
  summarize(propNA = mean(is.na(diff))) %>%
  ggplot(aes(x = cpg_type, y = propNA, color = cpg_type)) +
  geom_point() +
  facet_grid(samp2 ~ samp1) +
  theme_bw()
dev.off()




ranges <- apply(getBeta(INMA_comb), 1, function(x) {
  a <- quantile(x, c(0.99, 0.01))
  a[1] - a[2]
})
ranges_df <- tibble(cpg = names(ranges), range = ranges)

beta_diff_annot <- beta_diff %>% 
  left_join(colData(INMA_comb)[, c("Sample_Name", "idnum")] %>% data.frame() %>% mutate(samp1 = Sample_Name), by = "samp1") %>%
  left_join(colData(INMA_comb)[, c("Sample_Name", "idnum")] %>% data.frame() %>% mutate(samp2 = Sample_Name), by = "samp2") %>%
  select(-starts_with("Sample")) %>%
  left_join(ranges_df) %>%
  left_join(rel, by = "cpg")

png("figures/diff.beta.range.png")
beta_diff_annot %>%
  filter(idnum.x == idnum.y) %>%
  ggplot(aes(x = range, y = diff)) +
  geom_point() +
  facet_grid(~ idnum.x) +
  theme_bw()
dev.off()
  
beta_diff_annot %>%
  filter(idnum.x == idnum.y) %>%
  lm(diff ~ range + factor(idnum.x), .) %>%
  summary()

png("figures/tech.batch.range.png", width = 600, height = 800)

tech.batch.cpg %>%
  filter(Normalization == "Joint") %>%
  select(cpg, epi_type, method) %>%
  distinct() %>%
  right_join(ranges_df, by = "cpg") %>%
  mutate(epi_type = as.character(epi_type),
         epi_type = ifelse(is.na(epi_type), "Not Epimutations", epi_type)) %>%
  ggplot(aes(x = epi_type, y = range, fill = epi_type)) +
  geom_boxplot() +
  facet_grid(~ method) +
  theme_bw()
dev.off()




beta_diff_rel <- beta_diff %>%
  mutate()
  left_join(tech.batch.cpg)


b <- subset(all.res.df, method == "quantile" & chromosome != 0 & type == "Replicate different batch" & Normalization == "Joint")

png("figures/epi.beta.batch1.png", width = 600, height = 800)
plot_epimutations(b, INMA_comb)
dev.off()

b$sample <- "04_17_0"
png("figures/epi.beta.batch2.png", width = 600, height = 800)
plot_epimutations(b, INMA_comb)
dev.off()

beta.epimut <- subset(all.res.df, method == "beta" & chromosome != 0 & type == "Replicate different batch" & Normalization == "Joint")
pdf("figures/epimut.beta.joint.pdf")
lapply(seq_len(nrow(beta.epimut)), function(i){
  row <- beta.epimut[i, ]
  p1 <- plot_epimutations(row, INMA_comb)
  samp <- row$sample
  idnum <- colData(INMA_comb)[samp, "idnum"]
  samps <- colnames(INMA_comb[, INMA_comb$idnum == idnum])
  samp2 <- samps[samps != samp]
  row$sample <- samp2
  p2 <- plot_epimutations(row, INMA_comb)
  plot_grid(p1, p2)
})
dev.off()

manova.epimut <- subset(all.res.df, method == "manova" & chromosome != 0 & type == "Replicate different batch" & Normalization == "Joint" & idnum != 339)
pdf("figures/epimut.manova.joint.pdf")
lapply(seq_len(nrow(manova.epimut)), function(i){
  row <- manova.epimut[i, ]
  p1 <- plot_epimutations(row, INMA_comb)
  samp <- row$sample
  idnum <- colData(INMA_comb)[samp, "idnum"]
  samps <- colnames(INMA_comb[, INMA_comb$idnum == idnum])
  samp2 <- samps[samps != samp]
  row$sample <- samp2
  p2 <- plot_epimutations(row, INMA_comb)
  plot_grid(p1, p2)
})
dev.off()

# Same batch rep. - residuals ####
## Prepare data ####
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.Rdata")
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.liberal.Rdata")
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.liberal.Rdata")
load("INMA0combined.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
load("INMA_comb.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

ind.resid.df <- Reduce(rbind, res_indep_resid) %>%
  mutate(method = rep(names(res_indep_resid), sapply(res_indep_resid, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))


ind.comb.df <- rbind(mutate(ind.resid.df, QC = "Residuals"), 
                     mutate(ind.res.df, QC = "Raw"))


ind.resid.lib.df <- Reduce(rbind, res_indep_resid_lib) %>%
  mutate(method = rep(names(res_indep_resid_lib), sapply(res_indep_resid_lib, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))

ind.comb.resid.df <- rbind(mutate(ind.resid.lib.df, QC = "Residuals"), 
                     mutate(ind.res.lib.df, QC = "Raw"))


#'###############################################################################
# Borrar?
#'# Total epimutations ####
# tech.resid.num.tot.plot <- ind.comb.df %>%
#   filter(type == "Replicate same batch") %>%
#   group_by(QC, method) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = method, y = n, fill = QC)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected") +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "Preprocessing")
# 
# png("figures/INMA0.techRep.resid.totnum.png")
# tech.resid.num.tot.plot
# dev.off()
# Borrar?

## Epimutations per replicate ####
tech.resid.num.plot <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(QC, method, epi_type, idnum) %>% 
  complete(QC, method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 12, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ QC)


png("figures/INMA0.techRep.resid.num.png", width = 700)
tech.resid.num.plot
dev.off()


ind.resid.lib.replicates <- ind.comb.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id, QC) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, QC))

tech.resid.num.lib.plot <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id, QC) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, QC), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% ind.resid.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(QC, method, epi_type, idnum) %>% 
  complete(QC, method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 12, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ QC)

png("figures/INMA0.techRep.resid.num.lib.png", width = 700)
tech.resid.num.lib.plot
dev.off()

tech.resid.epi.plot <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, QC), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% ind.resid.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation")),
         epi_id = paste(idnum, epi_region_id)) %>%
  ggplot(aes(x = QC, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free_y", space = "free_y") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))


png("figures/INMA0.techRep.resid.num.epi.png", width = 1000, height = 700)
tech.resid.epi.plot
dev.off()


set_list2 <- list(Residuals = gset0_res[, gset0_res$Batch == "Esteller"], 
                 Raw = INMA_ind[, INMA_ind$Batch == "Esteller"])
tech.resid.top.plot <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, set_list2[[.[i, ]$QC]]))) %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Significant and top", "Replicate-specific epimutation"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both significant", "Significant and top", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(QC, method, epi_type, idnum) %>% 
  complete(QC, method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 12, 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ QC)

png("figures/INMA0.techRep.resid.top.png", height = 300)
tech.resid.top.plot
dev.off()




# Diff. norm. - tech rep. - residuals ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")
load("results/epimutations/Esteller.epimutations.normalization.residuals.liberal.Rdata")


norm.resid.df <- Reduce(rbind, 
                      lapply(INMA_norm_resid, function(x) {Reduce(rbind, x) %>%
                          mutate(method = rep(names(x), sapply(x, nrow)))}
                      )) %>%
  mutate(Normalization = rep(names(INMA_norm_resid), 
                             sapply(INMA_norm_resid, function(x) sum(sapply(x, nrow))))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist")) %>%
  rbind(., mutate(ind.resid.df, Normalization = "meffil")) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization"))) %>%
  filter(Sample_Name %in% unique(ind.resid.df$sample))



norm.resid.lib.df <- Reduce(rbind, 
                          lapply(INMA_norm_resid_lib, function(x) {Reduce(rbind, x) %>%
                              mutate(method = rep(names(x), sapply(x, nrow)))}
                          )) %>%
  mutate(Normalization = rep(names(INMA_norm_resid_lib), 
                             sapply(INMA_norm_resid_lib, function(x) sum(sapply(x, nrow)))),
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist")) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(Sample_Name %in% unique(ind.res.df$sample))
norm.resid.lib.df <- rbind(norm.resid.lib.df, 
                         mutate(ind.resid.lib.df, Normalization = "meffil") %>% filter(sample %in% norm.resid.lib.df$sample)) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization"))) 


#'###############################################################################
# Borrar?
#'# Total number of epimutations per normalization algorithm ####
# tech.totnum.norm.resid.plot <- norm.resid.df %>%
#   filter(type == "Replicate same batch") %>%
#   group_by(method, Normalization) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
#   geom_bar(stat = "identity") + theme_bw() +
#   facet_wrap(~ method)
# png("figures/INMA0.resid.techRep.Norm.totnum.png")
# tech.totnum.norm.resid.plot
# dev.off()
#'###############################################################################


## Technical replicates per normalization algorithm ####
tech.num.norm.resid.plot <- norm.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  mutate(epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(y = n + 1, label = n)) +
  theme_bw() +
  facet_grid(method ~ Normalization) +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "")


png("figures/INMA0.resid.techRep.Norm.num.png", width = 1000, height = 700)
tech.num.norm.resid.plot
dev.off()

norm.resid.lib.replicates <- norm.resid.lib.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", "Replicate-specific epimutation")) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization))


tech.num.norm.resid.lib.plot <- norm.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.resid.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(y = n + 1, label = n)) +
  theme_bw() +
  facet_grid(method ~ Normalization) +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png("figures/INMA0.resid.techRep.Norm.num.lib.png", width = 1000, height = 700)
tech.num.norm.resid.lib.plot
dev.off()


norm.raw.sum <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.res.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation")),
         epi_id = paste(idnum, epi_region_id)) 


norm.resid.sum <- norm.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both significant", "Replicate-specific epimutation")) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization), 
         epi_type = ifelse(epi_type == "Replicate-specific epimutation" & epi_name %in% norm.resid.lib.replicates$epi_name, "Significant and suggestive", epi_type),
         epi_type = factor(epi_type, levels = c("Both significant", "Significant and suggestive", "Replicate-specific epimutation")),
         epi_id = paste(idnum, epi_region_id)) 


tech.num.norm.resid.lib.epi <- rbind(mutate(norm.raw.sum, QC = "Raw"), 
      mutate(norm.resid.sum, QC = "Residuals")) %>%
  ggplot(aes(x = QC, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(Normalization ~ method, scales = "free_y", space = "free_y") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.resid.techRep.Norm.num.epi.png", width = 1000, height = 700)
tech.num.norm.resid.lib.epi
dev.off()


# 
# 
#'# Total number of epimutations per normalization algorithm in all Esteller samples ####
# tech.totnum.all.norm.resid.plot <- norm.resid.df %>%
#   filter(Batch == "Esteller") %>%
#   group_by(method, Normalization) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
#   geom_bar(stat = "identity") + theme_bw() +
#   facet_wrap(~ method)
# png("figures/INMA0.resid.techRep.Norm.totnum.allsamps.png")
# tech.totnum.all.norm.resid.plot
# dev.off()

# Diff. norm. - same samp- resid ####
tech.samesamp.num.norm.resid.plot <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Default parameters")


png("figures/INMA0.resid.sameSamp.Norm.propOverlap.png", height = 300)
tech.samesamp.num.norm.resid.plot
dev.off()


tech.samesamp.num.norm.resid.lot.lib <- norm.resid.lib.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  mutate(epi_id = paste(method, Sample_Name, epi_region_id)) %>%
  filter(epi_id %in% paste(norm.resid.df$method, norm.resid.df$Sample_Name, norm.resid.df$epi_region_id)) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Liberal parameters")

png("figures/INMA0.resid.sameSamp.Norm.propOverlap2.png", height = 300)
tech.samesamp.num.norm.resid.lot.lib
dev.off()

norm.raw.samp <- norm.res.lib.df %>%
  filter(chromosome != 0) %>%
  mutate(epi_name = paste(method, Sample_Name, epi_region_id),
         epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
  filter(epi_name %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id)) %>%
  mutate(epi_cat = ifelse(epi_name2 %in% paste(norm.res.df$method, norm.res.df$Sample_Name, norm.res.df$epi_region_id, norm.res.df$Normalization),
                          "Significant", "Suggestive"), 
         epi_id = paste(Sample_Name, epi_region_id),
         QC = "Raw")

norm.resid.samp <- norm.resid.lib.df %>%
  filter(chromosome != 0) %>%
  mutate(epi_name = paste(method, Sample_Name, epi_region_id),
         epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
  filter(epi_name %in% paste(norm.resid.df$method, norm.resid.df$Sample_Name, norm.resid.df$epi_region_id)) %>%
  mutate(epi_cat = ifelse(epi_name2 %in% paste(norm.resid.df$method, norm.resid.df$Sample_Name, norm.resid.df$epi_region_id, norm.resid.df$Normalization),
                          "Significant", "Suggestive"), 
         epi_id = paste(Sample_Name, epi_region_id),
         QC = "Residuals") 



tech.samesamp.norm.resid.epi.plot <- rbind(norm.raw.samp, norm.resid.samp) %>%
    ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
    geom_tile() +
    theme_bw() +
    facet_nested( ~ method + QC) + 
    scale_fill_discrete(name = "") +
    scale_y_discrete(name = "Epimutations") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle=90, vjust=0.5),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank())
  
  
png("figures/INMA0.resid.sameSamp.Norm.num.epi.png", width = 1000, height = 700)
tech.samesamp.norm.resid.epi.plot
dev.off()

# 
# 
# 
# tech.samesamp.num.norm.resid.plot3 <- norm.resid.df %>%
#   filter(chromosome != 0 & Batch == "Esteller" & !Normalization %in% c("QuantileNormalization", "IlluminaNormalization")) %>%
#   mutate(Normalization = droplevels(Normalization)) %>%
#   group_by(method, Sample_Name, epi_region_id) %>%
#   summarize(n_norms = length(unique(Normalization))) %>%
#   ungroup() %>%
#   count(method, n_norms) %>% 
#   complete(method, n_norms, fill = list(n = 0)) %>% 
#   group_by(method) %>%
#   mutate(n_norms = factor(n_norms),
#          prop = n/sum(n)*100) %>%
#   ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
#   scale_x_discrete(name = "Method") +
#   scale_fill_discrete("Norms") +
#   scale_color_discrete("Norms")
# 
# png("figures/INMA0.resid.sameSamp.Norm.propOverlap3.png")
# tech.samesamp.num.norm.resid.plot3
# dev.off()
# 
# 
# 
# png("figures/INMA0.resid.sameSamp.Norm.panel.png", width = 1200, height = 500)
# plot_grid(tech.totnum.all.norm.resid.plot, tech.samesamp.num.norm.resid.plot2, labels = "AUTO")
# dev.off()


# Diff. batch tech rep. - residuals ####
comb.resid.df <- Reduce(rbind, res_joint_resid) %>%
  mutate(method = rep(names(res_joint_resid), sapply(res_joint_resid, nrow))) %>%
  left_join(joint.pcs.df %>% 
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))


comSamps.resid <- intersect(unique(ind.resid.df$sample), unique(comb.resid.df$sample))

all.resid.df <- rbind(ind.resid.df, comb.resid.df) %>%
  filter(sample %in% comSamps)

#'# Total epimutations ####
# batch.totnum.norm.resid.plot <- all.resid.df %>%
#   filter(type == "Replicate different Batch") %>%
#   group_by(method, Normalization) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = method, y = n, fill = Normalization)) +
#   geom_bar(stat = "identity", position = "dodge") + theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected")
# 
# 
# png("figures/INMA0.resid.Batch.totnum.png")
# batch.totnum.norm.resid.plot
# dev.off()

# batch.ind.resid.plot <- all.resid.df %>%
#   filter(chromosome != 0 & type == "Replicate different batch") %>%
#   group_by(Normalization, method, idnum, epi_region_id) %>%
#   mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
#          epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
#                            ifelse(grepl("04", epi_type), "Alternative", "Reference")),
#          epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
#   group_by(Normalization) %>%
#   count(method, epi_type) %>% 
#   complete(method, epi_type, fill = list(n = 0)) %>%
#   ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Total epimutations detected") +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "") +
#   facet_grid(Normalization ~ ., scales = "free")
# 
# png("figures/INMA0.resid.batch.overlap.png")
# batch.ind.resid.plot
# dev.off()
# 

## Prepare data ####
comb.resid.lib.df <- Reduce(rbind, res_joint_resid_lib) %>%
  mutate(method = rep(names(res_joint_resid_lib), sapply(res_joint_resid_lib, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")),
         method = recode(method, "isoforest" = "iForest", "mahdistmcd" = "mah-dist"))

all.resid.lib <- rbind(comb.resid.lib.df, ind.resid.lib.df) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name)) %>%
  filter(epi_type == "Shared epimutation") %>%
  mutate(epi_name = paste(Normalization, method, idnum, epi_region_id))


## Epimutations per replicate ####
batch.ind.resid.boxplot <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Shared epimutation", Sample_Name)) %>%
  mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "blue", "red")) +
  facet_grid(Normalization ~ .)
png("figures/INMA0.resid.batch.overlap.boxplot.png")
batch.ind.resid.boxplot
dev.off()


batch.ind.resid.boxplot.lib <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name)) %>%
  mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(paste(Normalization, method, idnum, epi_region_id) %in% all.resid.lib$epi_name, 
                                  "Significant and suggestive", 
                                  ifelse(grepl("04", epi_type), "Alternative", "Reference"))),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Significant and suggestive", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 55)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "purple1", "blue", "red")) +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.resid.batch.overlap.boxplot.lib.png")
batch.ind.resid.boxplot.lib
dev.off()

all.res.comb <- rbind(mutate(all.resid.df, QC = "Residuals"), 
                      mutate(all.res.df, QC = "Raw"))

batch.ind.resid.epi.plot <- rbind(comb.resid.lib.df, ind.resid.lib.df) %>%
  filter(chromosome != 0 & type == "Replicate different batch")%>%
  mutate(epi_name = paste(method, idnum, epi_region_id),
         epi_name2 = paste(method, Sample_Name, epi_region_id, Normalization)) %>%
  filter(epi_name %in% paste(all.resid.df$method, all.resid.df$idnum, all.resid.df$epi_region_id)) %>%
  mutate(epi_cat = ifelse(epi_name2 %in% paste(all.resid.df$method, all.resid.df$Sample_Name, all.resid.df$epi_region_id, all.resid.df$Normalization),
                          "Significant", "Suggestive"), 
         epi_id = paste(idnum, epi_region_id),
         Batch2 = recode(Batch, "Esteller" = "Ref.", "MeDALL" = "Alt."),
         Batch2 = factor(Batch2, levels = c("Ref.", "Alt.")),
         norm = recode(Normalization, "Independent" = "Indep."),
         epi_norm = paste(norm, "-", Batch2)) %>%
  ggplot(aes(x = Batch2, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_nested(idnum ~ method + norm, scales = "free_y", space="free") + 
  scale_fill_discrete(name = "") +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
png("figures/INMA0.resid.batch.overlap.epi.png", width = 600, height = 700)
batch.ind.resid.epi.plot
dev.off()




set_list3 <- list(Joint = gsetcomb_res[, gsetcomb_res$Batch == "Esteller" | gsetcomb_res$dup], 
                 Independent = gset0_res[, gset0_res$Batch == "Esteller" | gset0_res$dup])
batch.resid.boxplot.top <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, set_list3[[.[i, ]$Normalization]]))) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(unique(Batch)) == 2, "Shared epimutation", Sample_Name), 
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(rep_quant > 0.95 | rep_quant < 0.05, "Significant and top", 
                                  ifelse(grepl("04", epi_type), "Alternative", "Reference"))),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Significant and top", "Reference", "Alternative"))) %>%
  select(Normalization, method, idnum, epi_region_id, epi_type) %>%
  distinct() %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 55)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "", values = c("purple4", "purple1", "blue", "red")) +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.resid.batch.overlap.boxplot.top.png")
batch.resid.boxplot.top
dev.off()


# 
# batch.ind.resid.prop.plot <-  all.res.comb %>% 
#   filter(chromosome != 0 & type == "Replicate different Batch") %>%
#   group_by(QC, Normalization, method, idnum, epi_region_id) %>%
#   mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
#          epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
#                            ifelse(grepl("04", epi_type), "Alternative", "Reference")),
#          epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
#   group_by(QC, Normalization, method) %>%
#   summarize(p = mean(epi_type == "Shared epimutation")) %>%
#   ggplot(aes(x = method, y = p*100, color = QC, fill = QC)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_bw() +
#   scale_y_continuous(name = "Proportion of shared epimutations", limits = c(0, 100)) +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "") +
#   facet_grid(Normalization ~ ., scales = "free")
# 
# png("figures/INMA0.resid.batch.overlap.prop.png")
# batch.ind.resid.prop.plot
# dev.off()
# 
# 
# batch.ind.resid.prop.boxplot <-  all.res.comb %>% 
#   filter(chromosome != 0 & type == "Replicate different Batch") %>%
#   group_by(QC, Normalization, method, idnum, epi_region_id) %>%
#   summarize(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name)) %>%
#   mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
#                            ifelse(grepl("04", epi_type), "Alternative", "Reference")),
#          epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
#   group_by(QC, Normalization, method, idnum) %>%
#   summarize(p = mean(epi_type == "Shared epimutation")) %>%
#   ggplot(aes(x = method, y = p*100, color = QC)) +
#   geom_boxplot() +
#   theme_bw() +
#   scale_y_continuous(name = "Proportion of shared epimutations", limits = c(0, 100)) +
#   scale_x_discrete(drop = FALSE) +
#   scale_fill_discrete(name = "") +
#   scale_color_discrete(name = "") +
#   facet_grid(Normalization ~ .)
# 
# png("figures/INMA0.resid.batch.overlap.prop.boxplot.png")
# batch.ind.resid.prop.boxplot
# dev.off()

png("figures/INMA0.resid.batch.panel.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.resid.plot, batch.ind.resid.plot, batch.ind.resid.prop.plot, labels = "AUTO")
dev.off()

png("figures/INMA0.resid.batch.panel2.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.resid.plot, batch.ind.resid.boxplot, batch.ind.resid.prop.boxplot, labels = "AUTO")
dev.off()

# Random #####
selbeta <- subset(ind.comb.df, method == "beta" & type == "Replicate same batch" & cpg_n >= 3)
beta_cpgs <- unique(unlist(strsplit(selbeta$cpg_ids, ",")))
beta_cpgs <- beta_cpgs[!is.na(beta_cpgs)]

load("INMA0combined.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

cot <- INMA_ind[, INMA_ind$idnum %in% c(423, 636)]
cot2 <- gset0_res[, gset0_res$idnum %in% c(423, 636)]
colnames(cot2) <- paste(colnames(cot2), "Resid", sep = "-")
cor(getBeta(cot[beta_cpgs[!is.na(beta_cpgs)], ]))

png("figures/cot.png", width = 1000, height = 700)
pheatmap(getBeta(cbind(cot[beta_cpgs, ], cot2[beta_cpgs, ])))
dev.off()

data.frame(selbeta[, c("Sample_Name", "QC", "epi_region_id", "cpg_ids")]) %>% arrange(Sample_Name)
getBeta(cbind(cot, cot2)[c("cg16526445","cg05038053","cg24419099"), c(1, 5, 4,8, 2,6, 3,7)])

cpgMeans_raw <- rowMeans(getBeta(INMA_ind))
correlsglobal_raw <- cor(getBeta(cot) - cpgMeans_raw)

cpgMeans_resid <- rowMeans(getBeta(gset0_res), na.rm = TRUE)
correlsglobal_resid <- cor(getBeta(cot2) - cpgMeans_resid, use = "complete")

load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
pc <- prcomp(getBeta(gset[1:500,]))
pc <- prcomp(t(meffil:::impute.matrix(getBeta(gset[1:500,]), margin = 1)))

getPCs <- function (normalized.beta, probe.range = 5000, autosomal = T, 
                    verbose = F) {
if (is.matrix(normalized.beta)) 
  sites <- rownames(normalized.beta)
else {
  gds.filename <- normalized.beta
  stopifnot(file.exists(gds.filename))
  gds.file <- openfn.gds(gds.filename)
  on.exit(closefn.gds(gds.file))
  sites <- read.gdsn(index.gdsn(gds.file, "row.names"))
}
subset <- sites
if (autosomal) {
  featureset <- meffil:::guess.featureset(sites)
  autosomal.sites <- meffil.get.autosomal.sites(featureset)
  subset <- intersect(autosomal.sites, subset)
}
meffil:::msg("Calculating variances", verbose = verbose)
if (is.matrix(normalized.beta)) {
  var.sites <- meffil.most.variable.cpgs(normalized.beta[subset, 
  ], n = probe.range)
}
else {
  cores <- getOption("mc.cores", 1)
  cl <- parallel::makeCluster(cores)
  vars <- clusterApply.gdsn(cl = cl, gds.fn = gds.filename, 
                            node.name = "matrix", margin = 1, as.is = "double", 
                            FUN = function(x) var(x, na.rm = T))
  vars <- vars[match(subset, sites)]
  var.sites <- subset[order(vars, decreasing = T)[1:probe.range]]
}
var.idx <- match(var.sites, sites)
meffil:::msg("Calculating beta PCs", verbose = verbose)
if (is.matrix(normalized.beta)) 
  mat <- normalized.beta[var.idx, ]
else {
  matrix.node <- index.gdsn(gds.file, "matrix")
  mat <- t(sapply(var.idx, function(idx) read.gdsn(matrix.node, 
                                                   start = c(idx, 1), count = c(1, -1))))
}
prcomp(t(meffil:::impute.matrix(mat, margin = 1)))
}
pc <- getPCs(getBeta(gset), probe.range = 40000)
m <- getM(gset)

beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)

ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1

gset.resid10 <- gset
gset.residNdim <- gset
gset.residCovars <- gset

res10 <- residuals(lmFit(m, pc$x[, 1:10]), m)
assay(gset.resid10) <- ilogit2(res10)

resNdim <- residuals(lmFit(m, pc$x[, 1:ndim]), m)
assay(gset.residNdim) <- ilogit2(resNdim)

library(epimutations)
load("INMA.commonControlSamples.Rdata")
par <- epi_parameters()
par$beta$pvalue_cutoff <- 0.01

beta_resid10 <- epimutations(case_samples = gset.resid10[, gset.resid10$idnum %in% c(423, 636)], 
                             control_panel = gset.resid10[, samps], method = "beta")
beta_resid10.all <- epimutations(case_samples = gset.resid10[, gset.resid10$idnum %in% c(423, 636)], 
                             control_panel = gset.resid10[, samps], method = "beta", epi_params = par)
beta_residNdim <- epimutations(case_samples = gset.residNdim[, gset.residNdim$idnum %in% c(423, 636)], 
                             control_panel = gset.residNdim[, samps], method = "beta")
beta_residNdim.all <- epimutations(case_samples = gset.residNdim[, gset.residNdim$idnum %in% c(423, 636)], 
                                 control_panel = gset.residNdim[, samps], method = "beta", epi_params = par)

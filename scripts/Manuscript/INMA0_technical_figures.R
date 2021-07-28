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

load("INMA0combined.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
INMA_ind <- gset
load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
INMA_comb <- gset

load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.Rdata")
load("results/epimutations/Esteller.epimutations.normalization.Rdata")
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.Rdata")
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.Rdata")
load("results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")


# Normalization PCAs ####
## Independent normalization ####
ind.pcs <- meffil.methylation.pcs(getBeta(INMA_ind), probe.range = 40000)
ind.pcs.df <- ind.pcs %>%
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
  mutate(type = ifelse(n == 1, "Replicate same Batch", "Replicate different Batch")) %>%
  data.frame()
rownames(idnum.tab.indep) <- idnum.tab.indep$idnum  
ind.pcs.df$type <- factor(ifelse(ind.pcs.df$dup, idnum.tab.indep[as.character(ind.pcs.df$idnum), "type"], 
                          ind.pcs.df$Batch2),
                          levels = c("Reference", "Alternative", "Replicate same Batch", "Replicate different Batch"))

indep_pc <- ggplot(ind.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  theme_bw() +
  ggtitle("Independent Normalization") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")

png("figures/INMA0.indepNorm.PCA.png")
indep_pc
dev.off()



## Joint normalization ####
joint.pcs <- meffil.methylation.pcs(getBeta(INMA_comb), probe.range = 40000)
joint.pcs.df <- joint.pcs %>%
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
  mutate(type = ifelse(n == 1, "Replicate same Batch", "Replicate different Batch")) %>%
  data.frame()
rownames(idnum.tab.joint) <- idnum.tab.joint$idnum  
joint.pcs.df$type <- factor(ifelse(joint.pcs.df$dup, idnum.tab.joint[as.character(joint.pcs.df$idnum), "type"], 
                                 joint.pcs.df$Batch2),
                          levels = c("Reference", "Alternative", "Replicate same Batch", "Replicate different Batch"))

joint_pc <- ggplot(joint.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  theme_bw() +
  ggtitle("Joint Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

png("figures/INMA0.jointNorm.PCA.png")
joint_pc
dev.off()


png("figures/INMA0.PCA_panel.png", width = 600, height = 300)
plot_grid(indep_pc, joint_pc, labels = "AUTO", nrow = 1, rel_widths = c(5, 8))
dev.off()

# Merge epimutations datasets ####
ind.res.df <- Reduce(rbind, res_indep) %>%
  mutate(method = rep(names(res_indep), sapply(res_indep, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")))

comb.res.df <- Reduce(rbind, res_joint) %>%
  mutate(method = rep(names(res_joint), sapply(res_joint, nrow))) %>%
  left_join(joint.pcs.df %>% 
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")))

comSamps <- intersect(unique(ind.res.df$sample), unique(comb.res.df$sample))

all.res.df <- rbind(ind.res.df, comb.res.df) %>%
  filter(sample %in% comSamps)


# Technical replicates in Esteller ####
## Total epimutations ####
ind.res.df %>%
  filter(type == "Replicate same Batch") %>%
  group_by(Normalization, method, idnum) %>%
  summarize(n = sum(chromosome != 0)) %>%
  group_by(method) %>%
  summarize(N = sum(n > 0)) 


tech.num.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "")


png("figures/INMA0.techRep.num.png")
tech.num.plot
dev.off()

## Epimutations per replicate ####
tech.ind.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("Rep1", epi_type), "Replicate 2", "Replicate 1")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate 1", "Replicate 2"))) %>%
  group_by(idnum) %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(idnum ~ .)

png("figures/INMA0.techRep.ind.png")
tech.ind.plot
dev.off()

png("figures/INMA0.techRep.panel.png", width = 1200, height = 500)
plot_grid(tech.num.plot, tech.ind.plot, labels = "AUTO")
dev.off()

# Merge epimutations datasets normalization ####
norm.res.df <- Reduce(rbind, 
                      lapply(INMA_norm, function(x) {Reduce(rbind, x) %>%
                                      mutate(method = rep(names(x), sapply(x, nrow)))}
                             )) %>%
  mutate(Normalization = rep(names(INMA_norm), 
                             sapply(INMA_norm, function(x) sum(sapply(x, nrow))))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  rbind(., mutate(ind.res.df, Normalization = "meffil")) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization")),
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd"))) %>%
  filter(Sample_Name %in% unique(ind.res.df$sample))

## Total number of epimutations per normalization algorithm ####
tech.totnum.norm.plot <- norm.res.df %>%
  filter(type == "Replicate same Batch") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
  geom_bar(stat = "identity") + theme_bw() +
  facet_wrap(~ method)
png("figures/INMA0.techRep.Norm.totnum.png")
tech.totnum.norm.plot
dev.off()



## Technical replicates per normalization algorithm ####
tech.num.norm.plot <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  facet_grid(Normalization ~ method) +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "")


png("figures/INMA0.techRep.Norm.num.png", width = 1000, height = 700)
tech.num.norm.plot
dev.off()

norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  group_by(Normalization,method) %>%
  summarize(p = mean(epi_type == "Shared epimutation")) %>%
  spread(Normalization, p) %>%
  data.frame()

png("figures/INMA0.techRep.Norm.panel.png", width = 1200, height = 500)
plot_grid(tech.totnum.norm.plot, tech.num.norm.plot, labels = "AUTO")
dev.off()


## Total number of epimutations per normalization algorithm in all Esteller samples ####
tech.totnum.all.norm.plot <- norm.res.df %>%
  filter(Batch == "Esteller") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
  geom_bar(stat = "identity") + theme_bw() +
  facet_wrap(~ method)
png("figures/INMA0.techRep.Norm.totnum.allsamps.png")
tech.totnum.all.norm.plot
dev.off()



## Same sample epimutations per normalization algorithm ####
tech.samesamp.num.norm.plot <- norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  mutate(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("Norms") +
  scale_color_discrete("Norms")

png("figures/INMA0.sameSamp.Norm.propOverlap.png")
tech.samesamp.num.norm.plot
dev.off()


tech.samesamp.num.norm.plot2 <- norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("Norms") +
  scale_color_discrete("Norms")

png("figures/INMA0.sameSamp.Norm.propOverlap2.png")
tech.samesamp.num.norm.plot2
dev.off()

png("figures/INMA0.sameSamp.Norm.panel.png", width = 1200, height = 500)
plot_grid(tech.totnum.all.norm.plot, tech.samesamp.num.norm.plot2, labels = "AUTO")
dev.off()

norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  select(-n) %>% spread(n_norms, prop)

# Technical replicates different batch ####
## Total epimutations ####
batch.totnum.norm.plot <- all.res.df %>%
  filter(type == "Replicate different Batch") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = method, y = n, fill = Normalization)) +
  geom_bar(stat = "identity", position = "dodge") + theme_bw() +
  scale_y_continuous(name = "Total epimutations detected")
  

png("figures/INMA0.Batch.totnum.png")
batch.totnum.norm.plot
dev.off()


## Epimutations per replicate ####
batch.ind.plot <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ ., scales = "free")

png("figures/INMA0.batch.overlap.png")
batch.ind.plot
dev.off()


batch.ind.boxplot <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name)) %>%
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
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ ., scales = "free")
png("figures/INMA0.batch.overlap.boxplot.png")
batch.ind.boxplot
dev.off()



png("figures/INMA0.batch.panel.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.plot, batch.ind.plot, labels = "AUTO")
dev.off()

png("figures/INMA0.batch.panel2.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.plot, batch.ind.boxplot, labels = "AUTO")
dev.off()

# Same batch replicates - residuals ####
ind.resid.df <- Reduce(rbind, res_indep_resid) %>%
  mutate(method = rep(names(res_indep_resid), sapply(res_indep_resid, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")))


ind.comb.df <- rbind(mutate(ind.resid.df, QC = "Residuals"), 
                     mutate(ind.res.df, QC = "Raw"))

## Total epimutations ####
tech.resid.num.tot.plot <- ind.comb.df %>%
  filter(type == "Replicate same Batch") %>%
  group_by(QC, method) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = method, y = n, fill = QC)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "Preprocessing")

png("figures/INMA0.techRep.resid.totnum.png")
tech.resid.num.tot.plot
dev.off()


tech.resid.num.plot <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch") %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(QC, method, epi_type) %>% 
  complete(QC, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(QC ~ .)


png("figures/INMA0.techRep.resid.num.png")
tech.resid.num.plot
dev.off()

png("figures/INMA0.techRep.resid.panel.png", width = 1200, height = 500)
plot_grid(tech.resid.num.tot.plot, tech.resid.num.plot, labels = "AUTO")
dev.off()


# Merge epimutations datasets normalization - residuals ####
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
  rbind(., mutate(ind.resid.df, Normalization = "meffil")) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("RawNormalization", "meffil", 
                                           "FunctionalNormalization", "IlluminaNormalization",
                                           "NoobNormalization", "QuantileNormalization",
                                           "SWANNormalization")),
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd"))) %>%
  filter(Sample_Name %in% unique(ind.resid.df$sample))


## Total number of epimutations per normalization algorithm ####
tech.totnum.norm.resid.plot <- norm.resid.df %>%
  filter(type == "Replicate same Batch") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
  geom_bar(stat = "identity") + theme_bw() +
  facet_wrap(~ method)
png("figures/INMA0.resid.techRep.Norm.totnum.png")
tech.totnum.norm.resid.plot
dev.off()



## Technical replicates per normalization algorithm ####
tech.num.norm.resid.plot <- norm.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same Batch" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  count(Normalization, method, epi_type) %>% 
  complete(Normalization, method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = epi_type, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  facet_grid(Normalization ~ method, scales = "free_y") +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "")


png("figures/INMA0.resid.techRep.Norm.num.png", width = 1000, height = 700)
tech.num.norm.resid.plot
dev.off()

norm.res.comb <- rbind(mutate(norm.resid.df, QC = "Residuals"), 
                       mutate(norm.res.df, QC = "Raw"))

norm.prop.rep.plot <- norm.res.comb %>%
  filter(chromosome != 0 & type == "Replicate same Batch" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(QC, Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", "Replicate-specific epimutation"),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Replicate-specific epimutation"))) %>%
  ungroup() %>%
  group_by(QC, Normalization,method) %>%
  summarize(p = mean(epi_type == "Shared epimutation")) %>%
  complete(QC, method, Normalization, fill = list(p = 0)) %>%
  ggplot(aes(x = method, y = p*100, color = QC, fill = QC)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  facet_grid(Normalization ~ .) +
  scale_y_continuous(name = "Proportion of shared epimutation (%)") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "")

png("figures/INMA0.resid.techRep.Norm.propRep.png", width = 1000, height = 700)
norm.prop.rep.plot
dev.off()



png("figures/INMA0.resid.techRep.Norm.panel.png", width = 1200, height = 500)
plot_grid(plot_grid(tech.totnum.norm.resid.plot, tech.num.norm.resid.plot, labels = "AUTO", nrow = 1),
          norm.prop.rep.plot, labels = c("", "C"), nrow = 2)
dev.off()


## Total number of epimutations per normalization algorithm in all Esteller samples ####
tech.totnum.all.norm.resid.plot <- norm.resid.df %>%
  filter(Batch == "Esteller") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = Normalization, y = n, fill = Normalization)) +
  geom_bar(stat = "identity") + theme_bw() +
  facet_wrap(~ method)
png("figures/INMA0.resid.techRep.Norm.totnum.allsamps.png")
tech.totnum.all.norm.resid.plot
dev.off()

## Same sample epimutations per normalization algorithm ####
tech.samesamp.num.norm.resid.plot <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  mutate(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("Norms") +
  scale_color_discrete("Norms")

png("figures/INMA0.resid.sameSamp.Norm.propOverlap.png")
tech.samesamp.num.norm.resid.plot
dev.off()


tech.samesamp.num.norm.resid.plot2 <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & Normalization != "QuantileNormalization") %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("Norms") +
  scale_color_discrete("Norms")

png("figures/INMA0.resid.sameSamp.Norm.propOverlap2.png")
tech.samesamp.num.norm.resid.plot2
dev.off()


tech.samesamp.num.norm.resid.plot3 <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller" & !Normalization %in% c("QuantileNormalization", "IlluminaNormalization")) %>%
  mutate(Normalization = droplevels(Normalization)) %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("Norms") +
  scale_color_discrete("Norms")

png("figures/INMA0.resid.sameSamp.Norm.propOverlap3.png")
tech.samesamp.num.norm.resid.plot3
dev.off()



png("figures/INMA0.resid.sameSamp.Norm.panel.png", width = 1200, height = 500)
plot_grid(tech.totnum.all.norm.resid.plot, tech.samesamp.num.norm.resid.plot2, labels = "AUTO")
dev.off()


# Technical replicates different batch - residuals ####
comb.resid.df <- Reduce(rbind, res_joint_resid) %>%
  mutate(method = rep(names(res_joint_resid), sapply(res_joint_resid, nrow))) %>%
  left_join(joint.pcs.df %>% 
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  mutate(Normalization = "Joint",
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd")))

comSamps.resid <- intersect(unique(ind.resid.df$sample), unique(comb.resid.df$sample))

all.resid.df <- rbind(ind.resid.df, comb.resid.df) %>%
  filter(sample %in% comSamps)

## Total epimutations ####
batch.totnum.norm.resid.plot <- all.resid.df %>%
  filter(type == "Replicate different Batch") %>%
  group_by(method, Normalization) %>%
  summarize(n = sum(chromosome != 0)) %>%
  ggplot(aes(x = method, y = n, fill = Normalization)) +
  geom_bar(stat = "identity", position = "dodge") + theme_bw() +
  scale_y_continuous(name = "Total epimutations detected")


png("figures/INMA0.resid.Batch.totnum.png")
batch.totnum.norm.resid.plot
dev.off()

## Epimutations per replicate ####
batch.ind.resid.plot <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ ., scales = "free")

png("figures/INMA0.resid.batch.overlap.png")
batch.ind.resid.plot
dev.off()



batch.ind.resid.boxplot <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name)) %>%
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
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ .)
png("figures/INMA0.resid.batch.overlap.boxplot.png")
batch.ind.resid.boxplot
dev.off()


all.res.comb <- rbind(mutate(all.resid.df, QC = "Residuals"), 
                                      mutate(all.res.df, QC = "Raw"))

batch.ind.resid.prop.plot <-  all.res.comb %>% 
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(QC, Normalization, method, idnum, epi_region_id) %>%
  mutate(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name),
         epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(QC, Normalization, method) %>%
  summarize(p = mean(epi_type == "Shared epimutation")) %>%
  ggplot(aes(x = method, y = p*100, color = QC, fill = QC)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of shared epimutations", limits = c(0, 100)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ ., scales = "free")

png("figures/INMA0.resid.batch.overlap.prop.png")
batch.ind.resid.prop.plot
dev.off()


batch.ind.resid.prop.boxplot <-  all.res.comb %>% 
  filter(chromosome != 0 & type == "Replicate different Batch") %>%
  group_by(QC, Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(Batch) == 2, "Shared epimutation", Sample_Name)) %>%
  mutate(epi_type = ifelse(epi_type == "Shared epimutation", "Shared epimutation",
                           ifelse(grepl("04", epi_type), "Alternative", "Reference")),
         epi_type = factor(epi_type, levels = c("Shared epimutation", "Reference", "Alternative"))) %>%
  group_by(QC, Normalization, method, idnum) %>%
  summarize(p = mean(epi_type == "Shared epimutation")) %>%
  ggplot(aes(x = method, y = p*100, color = QC)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Proportion of shared epimutations", limits = c(0, 100)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(name = "") +
  scale_color_discrete(name = "") +
  facet_grid(Normalization ~ .)

png("figures/INMA0.resid.batch.overlap.prop.boxplot.png")
batch.ind.resid.prop.boxplot
dev.off()

png("figures/INMA0.resid.batch.panel.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.resid.plot, batch.ind.resid.plot, batch.ind.resid.prop.plot, labels = "AUTO")
dev.off()

png("figures/INMA0.resid.batch.panel2.png", width = 1200, height = 500)
plot_grid(batch.totnum.norm.resid.plot, batch.ind.resid.boxplot, batch.ind.resid.prop.boxplot, labels = "AUTO")
dev.off()

## Random #####
selbeta <- subset(ind.comb.df, method == "beta" & type == "Replicate same Batch")
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

gset.resid10 <- gset
gset.resid90p <- gset

res10 <- residuals(lmFit(m, pc$x[, 1:10]), m)
assay(gset.resid10) <- ilogit2(res10)

res90p <- residuals(lmFit(m, pc$x[, 1:285]), m)
assay(gset.resid90p) <- ilogit2(res90p)

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

getMeanQuantile2 <- function(cpgs, sampid, set, rep_quant){
  if (!is.na(rep_quant)){
    return(rep_quant)
  }
  betas <- getBeta(set[cpgs, ])
  quant <- apply(betas, 1, function(x) {
    f <- ecdf(x)
    f(x[colnames(betas) == sampid])
  })
  
  mean(quant)
}

# Normalization PCAs ####
## Independent normalization ####
load("INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
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
load("INMA_comb.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
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


## Epimutations per replicate ####
tech.ind.plot <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", "One replicate")) %>%
  group_by(idnum) %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(name = "Algorithm", drop = FALSE) +
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "brown")) +
  facet_grid(idnum ~ .)

png("figures/INMA0.techRep.ind.png", height = 300)
tech.ind.plot
dev.off()


ind.res.tech.rep <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_ind[, INMA_ind$Batch == "Esteller"]))) 

tech.ind.top.plot <- ind.res.tech.rep %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
  ungroup() %>%
  count(method, idnum, epi_type) %>% 
  complete(method, idnum,  epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(name = "Algorithm", drop = FALSE) +
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  facet_grid(idnum ~ .)

png("figures/INMA0.techRep.ind.top.png", height = 300)
tech.ind.top.plot
dev.off()


ind.res.tech.rep %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
  ungroup() %>%
  count(method, epi_type) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  spread(epi_type, n)

ind.res.tech.rep %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
  ungroup() %>%
  count(method, idnum) %>% 
  complete(method, epi_type, fill = list(n = 0)) %>%
  spread(epi_type, n)


tech.rep.epi.plot <- ind.res.tech.rep %>%
group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate")), 
  epi_id = paste(idnum, epi_region_id)) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  scale_x_discrete(name = "Algorithm") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

png("figures/INMA0.techRep.ind.epi.png")
tech.rep.epi.plot
dev.off()

# Diff. norm. - tech rep. ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.Rdata")
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


## Technical replicates per normalization algorithm ####
paths <- dir(pattern = "Esteller.*Normalization.*autosomic.*withNA.*")
names(paths) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))
norm_set_list <- lapply(paths[-1], function(x) {
  load(x)
  gset
})
norm_set_list$meffil <- INMA_ind[, INMA_ind$Batch == "Esteller"]

norm.res.out <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, norm_set_list[[.[i, ]$Normalization]]))) 


tech.norm.plot <- norm.res.out %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_name = paste(method, idnum, epi_region_id, Normalization),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
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
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_color_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png("figures/INMA0.techRep.Norm.png", width = 1000, height = 700)
tech.norm.plot
dev.off()


tech.num.norm.epi.plot <- norm.res.out %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
 ggplot(aes(x = Normalization, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scale = "free_y") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))
  

png("figures/INMA0.techRep.Norm.epi.png", width = 1000, height = 700)
tech.num.norm.epi.plot
dev.off()


# Diff. norm. same samp. ####
## Bar plot ####
tech.samesamp.norm.quant <- norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  select(method, Sample_Name, epi_region_id, Normalization, cpg_ids) %>%
  mutate(rep_quant = 0) %>%
  group_by(method, Sample_Name, epi_region_id, cpg_ids) %>%
  complete(Normalization, fill = list(rep_quant = NA)) %>%
  ungroup() %>%
  mutate(rep_quant2 = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile2(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$Sample_Name, 
                     norm_set_list[[.[i, ]$Normalization]],
                     .[i, ]$rep_quant))) 


tech.samesamp.norm.plot <- norm.res.df %>%
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
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Only epimutations detected") + 
  theme(plot.title = element_text(hjust = 0.5))

png("figures/INMA0.sameSamp.Norm.propOverlap.png", height = 300)
tech.samesamp.norm.plot
dev.off()


tech.samesamp.norm.outlier.plot <- tech.samesamp.norm.quant %>%
  group_by(method, Sample_Name, epi_region_id, Normalization) %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm == min(rep_norm)) %>%
  select(method, Sample_Name, epi_region_id, rep_norm, Normalization) %>%
  distinct() %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = sum(rep_norm < 0.05)) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Epimutations detected and outlier signals") + 
  theme(plot.title = element_text(hjust = 0.5))

png("figures/INMA0.sameSamp.Norm.outlier.propOverlap.png", height = 300)
tech.samesamp.norm.outlier.plot
dev.off()


png("figures/INMA0.sameSamp.Norm.propOverlap.panel.png", width = 1100)
plot_grid(tech.samesamp.norm.plot,tech.samesamp.norm.outlier.plot, ncol = 2,
          labels = c("A", "B")) 
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


tech.samesamp.norm.quant %>%
  group_by(method, Sample_Name, epi_region_id, Normalization) %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm == min(rep_norm)) %>%
  select(method, Sample_Name, epi_region_id, rep_norm, Normalization) %>%
  distinct() %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = sum(rep_norm < 0.05)) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  select(-n) %>% spread(n_norms, prop)

## Epimutations plot ####
tech.samesamp.norm.epi.plot <- tech.samesamp.norm.quant %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm < 0.05) %>%
  mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
         epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_grid(Sample_Name ~ method, scales = "free_y", space = "free") + 
  scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "white")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.sameSamp.Norm.epi.png", height = 1000, width = 800)
tech.samesamp.norm.epi.plot
dev.off()


tech.samesamp.norm.epi2.plot <- tech.samesamp.norm.quant %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm < 0.05) %>%
  mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
         epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_wrap( ~ method, scales = "free_y") + 
  scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "white")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.sameSamp.Norm.epi2.png")
tech.samesamp.norm.epi2.plot
dev.off()

# Diff batch tech replicates ####
#'###############################################################################

## Prepare data ####
batch_list <- list(Joint = INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup], 
                 Independent = INMA_ind[, INMA_ind$Batch == "Esteller" | INMA_ind$dup])

batch.quant <- all.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, batch_list[[.[i, ]$Normalization]]))) 

## Boxplot ####
batch.boxplot <- batch.quant %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  mutate(n = pmin(n, 30)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "darkgreen", "blue")) +
  facet_grid(Normalization ~ ., scales = "free")

png("figures/INMA0.batch.boxplot.png")
batch.boxplot
dev.off()


batch.quant %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  group_by(Normalization, method, idnum) %>%
  mutate(p = n/sum(n)) %>%
  filter(!is.na(p)) %>%
  group_by(Normalization, method, epi_type) %>%
  summarize(mean = median(p)) %>%
  spread(epi_type, mean) %>%
  ungroup() %>%
  mutate(rep = as.vector(.[, 5] + .[, 6])) %>%
  data.frame()

## Epimutations plot ####
batch.ind.epi.plot <- batch.quant %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free_y", space="free") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "darkgreen", "blue")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
png("figures/INMA0.batch.epi.png", width = 600, height = 800)
batch.ind.epi.plot
dev.off()



# Same batch rep. - residuals ####
## Prepare data ####
load("results/epimutations/INMA_comb.epimutations.INMA0.duplicates.residuals.Rdata")
load("results/epimutations/INMA0combined.epimutations.INMA0.duplicates.residuals.Rdata")
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


## Epimutations per replicate ####

esteller_resid_list <- list(Residuals = gset0_res[, gset0_res$Batch == "Esteller"], 
                  Raw = INMA_ind[, INMA_ind$Batch == "Esteller"])

tech.resid.quant <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, esteller_resid_list[[.[i, ]$QC]]))) 


tech.resid.plot <- tech.resid.quant %>%
   group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
  ungroup() %>%
  count(QC, method, epi_type, idnum) %>% 
  complete(QC, method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 12, 2)) +
  scale_x_discrete(name = "Algorithm", drop = FALSE) +
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_color_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  facet_grid(idnum ~ QC)

png("figures/INMA0.techRep.resid.png", height = 300, width = 1000)
tech.resid.plot
dev.off()



## Epimutation plot ####
tech.resid.epi.plot <- tech.resid.quant %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate")),
         epi_id = paste(idnum, epi_region_id)) %>%
  ggplot(aes(x = QC, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free_y", space = "free_y") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))


png("figures/INMA0.techRep.resid.epi.png", width = 1000, height = 700)
tech.resid.epi.plot
dev.off()

# Diff. norm. - tech rep. - residuals ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")
load("Esteller.allminfiNormalizations.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

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

gset.residuals$meffil <- gset0_res[, gset0_res$Batch == "Esteller"]

norm.resid.out <- norm.resid.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, gset.residuals[[.[i, ]$Normalization]]))) 




## Technical replicates per normalization algorithm ####
tech.norm.resid.epi <- rbind(mutate(norm.res.out, QC = "Raw"), 
      mutate(norm.resid.out, QC = "Residuals")) %>%
  group_by(QC, method, idnum, epi_region_id, Normalization) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
  ggplot(aes(x = QC, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(Normalization ~ method, scales = "free_y", space = "free_y") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.resid.techRep.Norm.epi.png", width = 1000, height = 700)
tech.norm.resid.epi
dev.off()

tech.samesamp.norm.resid.quant <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  select(method, Sample_Name, epi_region_id, Normalization, cpg_ids) %>%
  mutate(rep_quant = 0) %>%
  group_by(method, Sample_Name, epi_region_id, cpg_ids) %>%
  complete(Normalization, fill = list(rep_quant = NA)) %>%
  ungroup() %>%
  mutate(rep_quant2 = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile2(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$Sample_Name, 
                     gset.residuals[[.[i, ]$Normalization]],
                     .[i, ]$rep_quant))) 


# Diff. norm. - same samp- resid ####
## Overlap ####
tech.samesamp.norm.resid.plot <- norm.resid.df %>%
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
  ggtitle("Only epimutations detected") + 
  theme(plot.title = element_text(hjust = 0.5))


png("figures/INMA0.resid.sameSamp.Norm.propOverlap.png", height = 300)
tech.samesamp.norm.resid.plot
dev.off()



tech.samesamp.norm.resid.outlier.plot <- tech.samesamp.norm.resid.quant %>%
  group_by(method, Sample_Name, epi_region_id, Normalization) %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm == min(rep_norm)) %>%
  select(method, Sample_Name, epi_region_id, rep_norm, Normalization) %>%
  distinct() %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = sum(rep_norm < 0.05)) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_discrete("N Norm. Algorithm") +
  scale_color_discrete("N Norm. Algorithm") +
  ggtitle("Epimutations detected and outlier signals") + 
  theme(plot.title = element_text(hjust = 0.5))

png("figures/INMA0.resid.sameSamp.Norm.outlier.propOverlap.png", height = 300)
tech.samesamp.norm.resid.outlier.plot
dev.off()


png("figures/INMA0.resid.sameSamp.Norm.propOverlap.panel.png", width = 1100)
plot_grid(tech.samesamp.norm.resid.plot,tech.samesamp.norm.resid.outlier.plot, ncol = 2,
          labels = c("A", "B")) 
dev.off()


## Epimutations plot ####
tech.samesamp.resid.norm.epi.plot <- rbind(mutate(tech.samesamp.norm.resid.quant, QC = "Residuals"),
                                           mutate(tech.samesamp.norm.quant, QC = "Raw")) %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm < 0.05) %>%
  mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
         epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_nested( ~ method + QC) + 
  scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "white")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.resid.sameSamp.Norm.epi.png", width = 1000, height = 700)
tech.samesamp.resid.norm.epi.plot
dev.off()



# Diff. batch tech rep. - residuals ####
## Normalization PCAs ####
### Independent normalization ####
ind.resid.pcs <- meffil.methylation.pcs(getBeta(gset0_res), probe.range = 40000, full.obj = TRUE)
ind.resid.pcs.df <- ind.resid.pcs$x %>%
  data.frame() %>%
  select(PC1, PC2) %>%
  mutate(Sample_Name = rownames(.)) %>%
  left_join(colData(gset0_res) %>% data.frame() %>% select(Sample_Name, Sex, dup, Batch, idnum), by = "Sample_Name") %>%
  mutate(Batch2 = ifelse(Batch == "Esteller", "Reference", "Alternative")) %>%
  as_tibble()

idnum.resid.tab.indep <- ind.resid.pcs.df %>%
  filter(dup) %>%
  group_by(idnum) %>%
  summarize(n = length(unique(Batch))) %>%
  mutate(type = ifelse(n == 1, "Replicate same batch", "Replicate different batch")) %>%
  data.frame()
rownames(idnum.resid.tab.indep) <- idnum.resid.tab.indep$idnum  
ind.resid.pcs.df$type <- factor(ifelse(ind.resid.pcs.df$dup, idnum.resid.tab.indep[as.character(ind.resid.pcs.df$idnum), "type"], 
                                 ind.resid.pcs.df$Batch2),
                          levels = c("Reference", "Alternative", "Replicate same batch", "Replicate different batch"))
ind.resid.pcs.vars <- ind.resid.pcs$sdev^2/sum(ind.resid.pcs$sdev^2)
indep_resid_pc <- ggplot(ind.resid.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  theme_bw() +
  ggtitle("Independent Normalization") +
  scale_x_continuous(name = paste0("PC1 (", round(ind.resid.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(ind.resid.pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")
summary(lm(PC1 ~ Batch, ind.resid.pcs.df))
summary(lm(PC2 ~ Batch, ind.resid.pcs.df))

## Joint normalization ####
joint.resid.pcs <- meffil.methylation.pcs(getBeta(gsetcomb_res), probe.range = 40000, full.obj = TRUE)
joint.resid.pcs.df <- joint.resid.pcs$x %>%
  data.frame() %>%
  select(PC1, PC2) %>%
  mutate(Sample_Name = rownames(.)) %>%
  left_join(colData(gsetcomb_res) %>% data.frame() %>% select(Sample_Name, Sex, dup, Batch, idnum), by = "Sample_Name") %>%
  mutate(Batch2 = ifelse(Batch == "Esteller", "Reference", "Alternative")) %>%
  as_tibble()

idnum.resid.tab.joint <- joint.resid.pcs.df %>%
  filter(dup) %>%
  group_by(idnum) %>%
  summarize(n = length(unique(Batch))) %>%
  mutate(type = ifelse(n == 1, "Replicate same batch", "Replicate different batch")) %>%
  data.frame()
rownames(idnum.resid.tab.joint) <- idnum.resid.tab.joint$idnum  
joint.resid.pcs.df$type <- factor(ifelse(joint.resid.pcs.df$dup, idnum.resid.tab.joint[as.character(joint.resid.pcs.df$idnum), "type"], 
                                   joint.resid.pcs.df$Batch2),
                            levels = c("Reference", "Alternative", "Replicate same batch", "Replicate different batch"))
joint.resid.pcs.vars <- joint.resid.pcs$sdev^2/sum(joint.resid.pcs$sdev^2)
joint_resid_pc <- ggplot(joint.resid.pcs.df, aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("blue", "red", "grey", "black")) +
  scale_x_continuous(name = paste0("PC1 (", round(joint.resid.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(joint.resid.pcs.vars[2]*100, 1), "%)")) +
  theme_bw() +
  ggtitle("Joint Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

summary(lm(PC1 ~ Batch, joint.resid.pcs.df))

png("figures/INMA0.resid.PCA_panel.png", width = 600, height = 300)
plot_grid(indep_resid_pc, joint_resid_pc, labels = "AUTO", nrow = 1, rel_widths = c(5, 8))
dev.off()



## Prepare data ####
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

res_batch_list <- list(Joint = gsetcomb_res[, gsetcomb_res$Batch == "Esteller" | gsetcomb_res$dup], 
                   Independent = gset0_res[, gset0_res$Batch == "Esteller" | gset0_res$dup])

batch.resid.quant <- all.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, res_batch_list[[.[i, ]$Normalization]]))) 



## Epimutations per replicate ####
batch.resid.boxplot <- batch.resid.quant %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "darkgreen", "blue")) +
  facet_grid(Normalization ~ .)
png("figures/INMA0.resid.batch.boxplot.png")
batch.resid.boxplot
dev.off()

batch.resid.quant %>%
  group_by(Normalization, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(Normalization) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  group_by(Normalization, method, idnum) %>%
  mutate(p = n/sum(n)) %>%
  filter(!is.na(p)) %>%
  group_by(Normalization, method, epi_type) %>%
  summarize(mean = median(p)) %>%
  spread(epi_type, mean) %>%
  ungroup() %>%
  mutate(rep = as.vector(.[, 5] + .[, 6])) %>%
  data.frame()


## Epimutations ####
batch.resid.epi.plot1 <-  rbind(mutate(batch.resid.quant, QC = "Residuals"), 
                      mutate(batch.quant, QC = "Raw")) %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  filter(idnum %in% c(17, 16, 339)) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(idnum, method, QC, epi_region_id, Normalization) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative")),
         Batch = paste(QC, Normalization)) %>%
  ggplot(aes(x = Batch, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free_y", space="free") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "darkgreen", "blue")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

batch.resid.epi.plot2 <-  rbind(mutate(batch.resid.quant, QC = "Residuals"), 
                                mutate(batch.quant, QC = "Raw")) %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  filter(!idnum %in% c(17, 16, 339)) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(idnum, method, QC, epi_region_id, Normalization) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative")),
         Batch = paste(QC, Normalization)) %>%
  ggplot(aes(x = Batch, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(idnum ~ method, scales = "free_y", space="free") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkgoldenrod2", "gray", "darkgreen", "blue")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
png("figures/INMA0.resid.batch.epi.png", width = 1200, height = 700)
plot_grid(batch.resid.epi.plot1, batch.resid.epi.plot2)
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

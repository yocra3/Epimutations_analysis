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
library(viridis)

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

getMeanDifference <- function(cpgs, ref, case){
  cpgs <- cpgs[cpgs %in% rownames(ref)]
  betas <- getBeta(ref[cpgs, ])
  means <- rowMedians(betas, na.rm = TRUE)
  diff  <- getBeta(case[cpgs, ]) - means
  mean(diff, na.rm = TRUE)
}


# Normalization PCAs (Sup Figure 11) ####
## Independent normalization ####
load("INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
load("INMA.commonControlSamples.Rdata")  
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
  scale_color_manual(name = "Batch", values = c("#F0E68C", "#228B22", "#4169E1", "#DC143C")) +
  geom_line(data = subset(ind.pcs.df, dup), aes(group = idnum)) +
  theme_bw() +
  ggtitle("Independent Normalization") +
  scale_x_continuous(name = paste0("PC1 (", round(ind.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(ind.pcs.vars[2]*100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")
summary(lm(PC1 ~ Batch, ind.pcs.df))

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
joint_pc <- joint.pcs.df %>%
  filter(!(type == "Replicate same batch" & Batch2 == "Alternative")) %>%
  ggplot(aes(x = PC1, y = PC2, color = type)) +
  geom_point() +
  scale_color_manual(name = "Batch", values = c("#F0E68C", "#228B22", "#4169E1", "#DC143C")) +
  geom_line(data = subset(joint.pcs.df, dup & !(type == "Replicate same batch" & Batch2 == "Alternative")), aes(group = idnum)) +
  scale_x_continuous(name = paste0("PC1 (", round(joint.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(joint.pcs.vars[2]*100, 1), "%)")) +
  theme_bw() +
  ggtitle("Combined Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

summary(lm(PC1 ~ Batch, joint.pcs.df))

## Sup Figure 11
png("figures/INMA0.PCA_panel.png", width = 2400, height = 1200, res = 300)
plot_grid(indep_pc, joint_pc, labels = "AUTO", nrow = 1, rel_widths = c(5, 8))
dev.off()



# Merge epimutations datasets ####
load("results/epimutations/INMA0combined.raw.epimutations.INMA0.duplicates.Rdata")

ind.res.df <- Reduce(rbind, res_indep) %>%
  mutate(method = rep(names(res_indep), sapply(res_indep, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "mlm")))

inma_ids <- sprintf("I%02d", seq_len(length(unique(ind.res.df$idnum))))
names(inma_ids) <- as.character(sort(unique(ind.res.df$idnum)))

ind.res.df$ID <- inma_ids[as.character(ind.res.df$idnum)]

# Technical replicates in Lab 1 ####
ind.res.tech.rep <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_ind[, INMA_ind$Batch == "Esteller"]))) 

## Figure 4A
tech.ind.top.plot <- ind.res.tech.rep %>%
  group_by(method, ID, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal"))) %>%
  ungroup() %>%
  count(method, ID, epi_type) %>% 
  complete(method, ID,  epi_type, fill = list(n = 0)) %>%
  ggplot(aes(x = method, y = n, fill = epi_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Total epimutations detected", breaks = seq(0, 10, 2)) +
  scale_x_discrete(name = "Algorithm", drop = FALSE) +
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1")) +
  facet_grid(. ~ ID)



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

## Sup Figure 8
tech.rep.epi.plot <- ind.res.tech.rep %>%
group_by(method, ID, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate")), 
  epi_id = paste(ID, epi_region_id)) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ ., scales = "free") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  scale_x_discrete(name = "Algorithm") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

png("figures/INMA0.techRep.ind.epi.png", height = 800, width = 1600, res = 300)
tech.rep.epi.plot
dev.off()

# Diff. norm. - tech rep. ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.Rdata")
load("results/epimutations/Esteller.epimutations.normalizationBMIQ.Rdata")

INMA_norm$BMIQ <- lapply(INMA_norm_BMIQ, function(x) x[, colnames(x) != "delta_beta"])
norm.res.df <- Reduce(rbind, 
                      lapply(INMA_norm, function(x) {Reduce(rbind, x) %>%
                                      mutate(method = rep(names(x), sapply(x, nrow)))}
                             )) %>%
  mutate(Normalization = rep(names(INMA_norm), 
                             sapply(INMA_norm, function(x) sum(sapply(x, nrow))))) %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(Normalization = gsub("Normalization", "", Normalization),
         Normalization = recode(Normalization, "Functional" = "Func-minfi"),
         method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(Sample_Name %in% unique(ind.res.df$sample)) 
norm.res.df$ID <- inma_ids[as.character(norm.res.df$idnum)]

norm.res.df <- rbind(norm.res.df, 
   mutate(ind.res.df, Normalization = "Func-meffil") %>% filter(sample %in% norm.res.df$sample)) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("Func-meffil", "Func-minfi", "Raw", "Illumina",
                                           "Noob", "Quantile", "SWAN", "BMIQ"))) 


## Technical replicates per normalization algorithm ####
paths <- dir(pattern = "Esteller.*Normalization.*normalizedComBat.autosomic.*withNA.*")
names(paths) <- gsub(".normalizedComBat.*$", "", gsub("Esteller.minfi", "", paths))
norm_set_list <- lapply(paths, function(x) {
  load(x)
  gset
})
load("Esteller.minfiRawNormalization.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
norm_set_list$BMIQ <- gset
norm_set_list$meffil <- INMA_ind[, INMA_ind$Batch == "Esteller"]
names(norm_set_list) <- gsub("Normalization", "", names(norm_set_list))
names(norm_set_list)[names(norm_set_list) == "meffil"] <- "Func-meffil"
names(norm_set_list)[names(norm_set_list) == "Functional"] <- "Func-minfi"
colData(norm_set_list$BMIQ) <- colData(norm_set_list$SWAN )

norm.res.out <- norm.res.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, norm_set_list[[as.character(.[i, ]$Normalization)]])
  ) )

### Sup Figure 9
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
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "grey")) +
  scale_color_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "grey")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

png("figures/INMA0.techRep.Norm.png", width = 3000, height = 1000, res = 300)
tech.norm.plot
dev.off()

## Sup Figure 10
tech.num.norm.epi.plot <- norm.res.out %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  group_by(Normalization, method, ID, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_id = paste(ID, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate"))) %>%
 ggplot(aes(x = Normalization, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ method, scale = "free_y") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "grey")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))
  

png("figures/INMA0.techRep.Norm.epi.png", width = 2400, height = 1400, res = 300)
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

norm.res.df_mean <-  norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  mutate(diff_median = sapply(seq_len(nrow(.)), function(i) 
      getMeanDifference(strsplit(.[i, ]$cpg_ids, ",")[[1]], 
                        ref = norm_set_list[[.[i, ]$Normalization]][, samps], 
                        case = norm_set_list[[.[i, ]$Normalization]][, .[i, ]$Sample_Name])))


### Sup Figure 9
tech.samesamp.norm.plot <- norm.res.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms, levels = 8:1),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_viridis("N Norm. Algorithm", discrete = TRUE) +
  scale_color_viridis("N Norm. Algorithm", discrete = TRUE) +
  ggtitle("Only epimutations detected") + 
  theme(plot.title = element_text(hjust = 0.5))


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
  mutate(n_norms = factor(n_norms, levels = 8:1),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_viridis("N Norm. Algorithm", discrete = TRUE) +
  scale_color_viridis("N Norm. Algorithm", discrete = TRUE) +
  ggtitle("Epimutations detected and outlier signals") + 
  theme(plot.title = element_text(hjust = 0.5))


png("figures/INMA0.sameSamp.Norm.propOverlap.panel.png", width = 3300, height = 1200, res = 300)
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

## Epimutations plot (Sup Figure 7) ####
inma_ids2 <- sprintf("S%02d", seq_len(length(unique(tech.samesamp.norm.quant$Sample_Name))))
names(inma_ids2) <- as.character(sort(unique(tech.samesamp.norm.quant$Sample_Name)))
tech.samesamp.norm.quant$ID <- inma_ids2[tech.samesamp.norm.quant$Sample_Name]

tech.samesamp.norm.epi.plot1 <- tech.samesamp.norm.quant %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm < 0.05) %>%
  filter(Sample_Name %in% c("SAB_C_0017", "SAB_C_0244", "SAB_C_0636", "SAB_C_0636_Rep1",
                            "SAB_C_0016", "SAB_C_0120", "SAB_C_0211")) %>%
  mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
         epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
         Sample_Name = gsub("SAB_C_", "", Sample_Name),
         Sample_Name = ifelse(Sample_Name == "0636", "0636_Rep2", 
                              ifelse(Sample_Name == "0423", "0423_Rep2", Sample_Name)),
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ method, scales = "free_y", space = "free") + 
  scale_fill_manual(name = "", values = c("darkorchid4", "darkorchid1", "white")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none")
tech.samesamp.norm.epi.plot2 <- tech.samesamp.norm.quant %>%
  mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
  filter(rep_norm < 0.05) %>%
  filter(!Sample_Name %in% c("SAB_C_0017", "SAB_C_0244", "SAB_C_0636", "SAB_C_0636_Rep1",
                            "SAB_C_0016", "SAB_C_0120", "SAB_C_0211")) %>%
  mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
         epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
         Sample_Name = gsub("SAB_C_", "", Sample_Name),
         Sample_Name = ifelse(Sample_Name == "0636", "0636_Rep2", 
                              ifelse(Sample_Name == "0423", "0423_Rep2", Sample_Name)),
         epi_id = paste(Sample_Name, epi_region_id)) %>%
  ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ method, scales = "free_y", space = "free") + 
  scale_fill_manual(name = "", values = c("darkorchid4", "darkorchid1", "white")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())


png("figures/INMA0.sameSamp.Norm.epi.png", height = 3200, width = 2800, res = 300)
plot_grid(tech.samesamp.norm.epi.plot1, tech.samesamp.norm.epi.plot2, ncol = 2, rel_widths = c(3, 4))
dev.off()


# Diff batch technical replicates ####
#'###############################################################################

## Prepare data ####
batch_list <- list(Joint = INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup], 
                 Independent = INMA_ind[, INMA_ind$Batch == "Esteller" | INMA_ind$dup])

ind.techbatch.df  <- ind.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, batch_list[[.[i, ]$Normalization]]))) 

# Same batch rep. - residuals ####
## Prepare data ####
load("results/epimutations/INMA0combined.raw.epimutations.INMA0.duplicates.residuals.Rdata")
load("INMA0combined.raw.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
load("INMA_comb.raw.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

ind.resid.df <- Reduce(rbind, res_indep_resid) %>%
  mutate(method = rep(names(res_indep_resid), sapply(res_indep_resid, nrow))) %>%
  left_join(ind.pcs.df %>%
              select(Sample_Name, type, idnum, Batch) %>% 
              mutate(sample = Sample_Name), by = "sample") %>%
  filter(!is.na(Batch)) %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(Normalization = "Independent",
         method = factor(method, levels = c("quantile", "beta", "mlm")))

ind.resid.df$ID <- inma_ids[as.character(ind.resid.df$idnum)]

ind.comb.df <- rbind(mutate(ind.resid.df, QC = "Residuals"), 
                     mutate(ind.res.df, QC = "Raw"))


## Epimutations per replicate ####
esteller_resid_list <- list(Residuals = gset0_res[, gset0_res$Batch == "Esteller"], 
                  Raw = INMA_ind[, INMA_ind$Batch == "Esteller"])

tech.resid.quant <- ind.comb.df %>%
  filter(chromosome != 0 & type == "Replicate same batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, esteller_resid_list[[.[i, ]$QC]]))) 


## Epimutation plot (Sup Figure 15) ####
tech.resid.epi.plot <- tech.resid.quant %>%
  group_by(QC, method, ID, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", "One replicate"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "One replicate")),
         epi_id = paste(ID, epi_region_id)) %>%
  ggplot(aes(x = QC, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ method, scales = "free_y", space = "free_y") + 
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "grey")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))


png("figures/INMA0.techRep.resid.epi.png", width = 2000, height = 1200, res = 300)
tech.resid.epi.plot
dev.off()

# Diff. norm. - tech rep. - residuals ####
## Prepare data ####
load("results/epimutations/Esteller.epimutations.normalization.residuals.Rdata")
load("Esteller.allminfiNormalizations.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
load("results/epimutations/Esteller.epimutations.BMIQ.residuals.Rdata")
load("Esteller.BMIQ.residuals.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

INMA_norm_resid$BMIQ <- lapply(INMA_norm_resid_BMIQ, function(x) x[, colnames(x) != "delta_beta"])

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
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(Normalization = gsub("Normalization", "", Normalization),
         Normalization = recode(Normalization, "Functional" = "Func-minfi"),
         method = factor(method, levels = c("quantile", "beta", "mlm"))) 

norm.resid.df$ID <- inma_ids[as.character(norm.resid.df$idnum)]

norm.resid.df <- rbind(norm.resid.df, mutate(ind.resid.df, Normalization = "Func-meffil")) %>%
  mutate(Normalization = factor(Normalization, 
                                levels = c("Func-meffil", "Func-minfi", "Raw", "Illumina",
                                           "Noob", "Quantile", "SWAN", "BMIQ"))) %>%
  filter(Sample_Name %in% unique(ind.resid.df$sample))


gset.residuals$BMIQ <- gset.residuals.BMIQ
gset.residuals$`Func-meffil` <- gset0_res[, gset0_res$Batch == "Esteller"]
names(gset.residuals) <- gsub("Normalization", "", names(gset.residuals))
names(gset.residuals)[names(gset.residuals) == "Functional"] <- "Func-minfi"
colData(gset.residuals$BMIQ) <- colData(gset.residuals$SWAN )

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
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "gray")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

## Sup Figure 16
png("figures/INMA0.resid.techRep.Norm.epi.png", width = 2800, height = 2400, res = 300)
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
## Overlap (Sup Figure 14) ####
tech.samesamp.norm.resid.plot <- norm.resid.df %>%
  filter(chromosome != 0 & Batch == "Esteller") %>%
  group_by(method, Sample_Name, epi_region_id) %>%
  summarize(n_norms = length(unique(Normalization))) %>%
  ungroup() %>%
  count(method, n_norms) %>% 
  complete(method, n_norms, fill = list(n = 0)) %>% 
  group_by(method) %>%
  mutate(n_norms = factor(n_norms, levels = 8:1),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)", limits = c(0, 100)) +
  scale_x_discrete(name = "Method") +
  scale_fill_viridis("N Norm. Algorithm", discrete = TRUE) +
  scale_color_viridis("N Norm. Algorithm", discrete = TRUE) +
  ggtitle("Only epimutations detected") + 
  theme(plot.title = element_text(hjust = 0.5))


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
  mutate(n_norms = factor(n_norms, levels = 8:1),
         prop = n/sum(n)*100) %>%
  ggplot(aes(x = method, y = prop, color = n_norms, fill = n_norms)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of epimutations (%)") +
  scale_x_discrete(name = "Method") +
  scale_fill_viridis("N Norm. Algorithm", discrete = TRUE) +
  scale_color_viridis("N Norm. Algorithm", discrete = TRUE) +
  ggtitle("Epimutations detected and outlier signals") + 
  theme(plot.title = element_text(hjust = 0.5))

##  Sup Figure 14
png("figures/INMA0.resid.sameSamp.Norm.propOverlap.panel.png", width = 3300, height = 1200, res = 300)
plot_grid(tech.samesamp.norm.resid.plot,tech.samesamp.norm.resid.outlier.plot, ncol = 2,
          labels = c("A", "B")) 
dev.off()


## Epimutations plot ####
# tech.samesamp.resid.norm.epi.plot <- rbind(mutate(tech.samesamp.norm.resid.quant, QC = "Residuals"),
#                                            mutate(tech.samesamp.norm.quant, QC = "Raw")) %>%
#   mutate(rep_norm = pmin(rep_quant2, 1 - rep_quant2)) %>%
#   filter(rep_norm < 0.05) %>%
#   mutate(epi_cat = ifelse(!is.na(rep_quant), "Epimutation", "Outlier signal"), 
#          epi_cat = factor(epi_cat, levels = c("Epimutation", "Outlier signal", "No signal")),
#          epi_id = paste(Sample_Name, epi_region_id)) %>%
#   ggplot(aes(x = Normalization, y = epi_id, fill = epi_cat)) +
#   geom_tile() +
#   theme_bw() +
#   facet_nested( ~ method + QC) + 
#   scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "white")) +
#   scale_y_discrete(name = "Epimutations") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.x = element_text(angle=90, vjust=0.5),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank())
# 
# 
# png("figures/INMA0.resid.sameSamp.Norm.epi.png", width = 4000, height = 2800, res = 300)
# tech.samesamp.resid.norm.epi.plot
# dev.off()

norm.resid.df %>%
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


# Diff. batch tech rep. - residuals ####
## Normalization PCAs ####
### Independent normalization (Sup Figure 12) ####
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
  geom_line(data = subset(ind.resid.pcs.df, dup & !(type == "Replicate same batch" & Batch2 == "Alternative")), aes(group = idnum)) +
  scale_color_manual(name = "Batch", values = c("#F0E68C", "#228B22", "#4169E1", "#DC143C")) +
  theme_bw() +
  ggtitle("After residuals extraction") +
  scale_x_continuous(name = paste0("PC1 (", round(ind.resid.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(ind.resid.pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5))
summary(lm(PC1 ~ Batch, ind.resid.pcs.df))
summary(lm(PC2 ~ Batch, ind.resid.pcs.df))

## Sup Figure 11
png("figures/INMA0.indepNorm.resid.PCA.png", height = 1200, width = 1600, res = 300)
indep_resid_pc
dev.off()

## Prepare data ####
res_batch_list <- list(Joint = gsetcomb_res[, gsetcomb_res$Batch == "Esteller" | gsetcomb_res$dup], 
                   Independent = gset0_res[, gset0_res$Batch == "Esteller" | gset0_res$dup])

batch.resid.quant <- ind.resid.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, res_batch_list[[.[i, ]$Normalization]]))) 



## Figure 4B
batch.resid.epiprop.plot <- rbind(mutate(batch.resid.quant, QC = "Residuals"), 
                             mutate(ind.techbatch.df, QC = "Raw")) %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(QC) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  group_by(idnum, QC, method) %>%
  filter(any(n > 0)) %>%
  mutate(p = n/sum(n)) %>%
   ggplot(aes(x = QC, y = p*100, color = epi_type)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ method) +
  scale_y_continuous(name = "Proportion of epimutations per individual (%)") +
  scale_color_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "#F0E68C", "#228B22"))
  

batch.epiprop <- rbind(mutate(batch.resid.quant, QC = "Residuals"), 
                                  mutate(ind.techbatch.df, QC = "Raw")) %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  group_by(QC, method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal", 
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(idnum, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  group_by(QC) %>%
  count(method, epi_type, idnum) %>% 
  complete(method, epi_type, idnum, fill = list(n = 0)) %>%
  group_by(idnum, QC, method) %>%
  filter(any(n > 0)) %>%
  mutate(p = n/sum(n)) %>%
  summarize(p_epi = p[epi_type == "Both replicates"],
            p_signal = sum(p[epi_type %in% c("Both replicates","One replicate and outlier signal")]))


lapply(unique(batch.epiprop$method), function(x) {
  wilcox.test(p_epi ~ QC, subset(batch.epiprop, method == x), conf.int = TRUE, exact = FALSE)
})

# ## Epimutations (Sup Figure 13) ####
batch.resid.epi.plot1 <-  batch.resid.quant %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  filter(!idnum %in% c(244, 484, 339)) %>%
  group_by(ID, method, epi_region_id, Normalization) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates",
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal",
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(ID, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ ., scales = "free_y", space="free") +
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "#F0E68C", "#228B22")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none")

batch.resid.epi.plot2 <-  batch.resid.quant %>%
  mutate(rep_quant = ifelse(is.na(rep_quant), 0.5, rep_quant)) %>%
  filter(idnum %in% c(244, 484, 339)) %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  group_by(ID, method, epi_region_id, Normalization) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Both replicates",
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "One replicate and outlier signal",
                                     ifelse(grepl("04", unique(Sample_Name)), "Alternative", "Reference")))) %>%
  mutate(epi_id = paste(ID, epi_region_id),
         epi_type = factor(epi_type, levels = c("Both replicates", "One replicate and outlier signal", "Reference", "Alternative"))) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  facet_grid(ID ~ ., scales = "free_y", space="free") +
  scale_fill_manual(name = "Epimutation detection", values = c("darkorchid4", "darkorchid1", "#F0E68C", "#228B22")) +
  scale_y_discrete(name = "Epimutations") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

png("figures/INMA0.resid.batch.epi.png", width = 2800, height = 2400, res = 300)
plot_grid(batch.resid.epi.plot1, batch.resid.epi.plot2, rel_widths = c(2, 3))
dev.off()



# Panel of technical results (Figure 4) ####
grid <- plot_grid(tech.ind.top.plot + theme(plot.margin = margin(7, 20, 7, 7)), batch.resid.epiprop.plot + theme(plot.margin = margin(7, 20, 7, 7)), ncol = 2, labels = LETTERS[1:2])
ggsave(file = "figures/INMA0.tech.panel.eps", plot = grid, width = 8000, height = 2700, dpi = 600, units = "px")


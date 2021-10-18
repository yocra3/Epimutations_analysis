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
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Epimutation detected in both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Epimutation detected in one replicate and signal in the other", "Epimutation detected in one replicate. No signal in the other"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Epimutation detected in both replicates", "Epimutation detected in one replicate and signal in the other", "Epimutation detected in one replicate. No signal in the other")), 
         epi_id = paste(idnum, epi_region_id)) 
tech.rep.plot <- tech.rep.res %>%
  mutate(idnum = ifelse(idnum == 423, "Individual 1", "Individual 2")) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  ggtitle("Same batch") +
  facet_grid(idnum ~ method, scales = "free") + 
  scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "Epimutations") +
  scale_x_discrete(name = "Algorithm") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 28),
        axis.title = element_text(size = 24),
        axis.text.x  = element_text(size = 16),
        legend.position = "none") 



## Batch replicates ####
batch.rep.res <- comb.res.df %>%
  filter(chromosome != 0 & type == "Replicate different batch") %>%
  mutate(rep_quant = sapply(seq_len(nrow(.)), function(i) 
    getMeanQuantile(strsplit(.[i, ]$cpg_ids, ",")[[1]], .[i, ]$sample, 
                    INMA_comb[, INMA_comb$Batch == "Esteller" | INMA_comb$dup]))) %>%
  group_by(method, idnum, epi_region_id) %>%
  summarize(epi_type = ifelse(length(unique(Sample_Name)) == 2, "Epimutation detected in both replicates", 
                              ifelse(any(rep_quant > 0.95 | rep_quant < 0.05) , "Epimutation detected in one replicate and signal in the other", "Epimutation detected in one replicate. No signal in the other"))) %>%
  mutate(epi_type = factor(epi_type, levels = c("Epimutation detected in both replicates", "Epimutation detected in one replicate and signal in the other", "Epimutation detected in one replicate. No signal in the other")), 
         epi_id = paste(idnum, epi_region_id)) 

batch.rep.plot <- batch.rep.res %>%
  filter(idnum != 339) %>%
  ungroup() %>%
  mutate(ind = factor(idnum, labels = paste("Ind.", 1:8))) %>%
  ggplot(aes(x = method, y = epi_id, fill = epi_type)) +
  geom_tile() +
  theme_bw() +
  ggtitle("Different batch") +
  facet_grid(ind ~ method, scales = "free") + 
  scale_fill_manual(name = "", values = c("darkgoldenrod2", "gray", "brown")) +
  scale_y_discrete(name = "") +
  scale_x_discrete(name = "Algorithm") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 28),
        axis.title = element_text(size = 24),
        axis.text.x  = element_text(size = 16),
        legend.text = element_text(size = 20),
        legend.position = "none") 

title <- ggdraw() + draw_label("Technical Replicates", fontface='bold', size = 32, 
                               hjust = 0.5)

legend_b <- get_legend(
  batch.rep.plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "right")
)

png("figures/poster/INMA0.Panel.png", width = 1300, height = 650)
p <- plot_grid(tech.rep.plot, batch.rep.plot, labels = "", nrow = 1)
plot_grid(title, p, legend_b, ncol = 1, rel_heights = c(0.1, 1, 0.15)) 
dev.off()

# INMA results ####
## 0 years
load("INMA_comb.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gset.0 <- gset

## 4 years
load("MeDALL_all.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

gset.sab <- gset[, gset$Sample_Group == "Age_4" & grepl("^04|SAB", gset$Sample_Name)]
gset.4  <- gset.sab[, !duplicated(gset.sab$idnum)]

## 8 years
load("results/preprocess/HELIX/HELIX.withNA.GenomicRatioSet.Rdata")
helix <- gset[, gset$cohort != "MOBA" & gset$h_ethnicity_3cat == "WhiteEur_WhiteOther"]

## Define functions
make_res_df <- function(res.list, pheno, dataset){
  nMethod <- sapply(res.list, nrow)
  nMethod[sapply(nMethod, is.null)] <- 0
  methods <- names(res.list)
  names(methods) <- methods
  Reduce(rbind, res.list) %>%
    mutate(method = rep(methods, unlist(nMethod)),
           dataset = dataset) %>%
    left_join(pheno %>% 
                data.frame() %>% 
                select(Sample_Name, Sex, idnum) %>% 
                mutate(sample = Sample_Name,
                       idnum = as.character(idnum)), by = "sample")
}
make_sum_df <- function(res_df){
  res_df %>%
    group_by(method, sample, Sex, idnum, dataset) %>%
    summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
    mutate(n_cat = ifelse(n > 1, "1+", as.character(n)))
}


load("results/epimutations/INMA_comb.epimutations.allSamples.Rdata")
load("results/epimutations/INMA4.epimutations.allSamples.Rdata")
load("results/epimutations/HELIX.epimutations.allSamples.Rdata")

helixp <- colData(helix)
helixp$idnum <- sapply(strsplit(helixp$SampleID, "_"), `[`, 2)
helixp$Sex <- helixp$e3_sex
helixp$Sample_Name <- helixp$SampleID

res.inma0.df <- make_res_df(res.inma0.list, pheno = colData(gset.0), dataset = "CordBlood") 
inma0.sum.df <- make_sum_df(res.inma0.df)


res.inma4.df <- make_res_df(res.inma4.list, pheno = colData(gset.4), dataset = "4years") 
inma4.sum.df <- make_sum_df(res.inma4.df)

res.helix.df <- make_res_df(res.helix.list, pheno = helixp, dataset = "8years") 
helix.sum.df <- make_sum_df(res.helix.df)


sab.res.df <- rbind(res.inma0.df, res.inma4.df, res.helix.df) %>%
  mutate(dataset = factor(dataset, levels = c("CordBlood", "4years", "8years"))) %>%
  filter(dataset != "8years" | grepl("SAB", sample))
methods <- unique(sab.res.df$method)
names(methods ) <- methods

selSamps <- sab.res.df %>%
  select(dataset, idnum) %>%
  distinct() %>%
  group_by(idnum) %>%
  summarize(n = n()) %>%
  filter(n == 3)

sab.rep.res <- sab.res.df %>%
  mutate(method = ifelse(method == "isoforest", "iForest", method)) %>%
  filter(idnum %in% selSamps$idnum) %>%
  filter(chromosome != 0) %>%
  group_by(dataset, method, idnum, epi_region_id) %>%
  summarize(cpg_ids = paste(unique(cpg_ids), collapse = ",")) %>%
  spread(dataset, cpg_ids)
sab.rep.res$cpg_list = lapply(seq_len(nrow(sab.rep.res)), function(i) {
  cpg_string <- paste(sab.rep.res[i, c("CordBlood", "4years", "8years")], collapse = ",")
  cpg_list <- strsplit(cpg_string, ",")[[1]]
  unique(cpg_list[cpg_list != "NA"])
})
  
getMeanQuantile2 <- function(cpgs, idnum, set){
  samps <- colnames(set[, set$idnum == idnum])
  cpgs <- cpgs[cpgs %in% rownames(set)]
  betas <- getBeta(set[cpgs, ])
  quant <- apply(betas, 1, function(x) {
    f <- ecdf(x)
    f(x[colnames(betas) == samps])
  })
  
  mean(quant)
}

getCord <- function(i) {
  if (!is.na(sab.rep.res[i, ]$CordBlood)){
    return(0)
  } else{
    getMeanQuantile2(sab.rep.res[i, ]$cpg_list[[1]], sab.rep.res[i, ]$idnum, gset.0) 
  }
}

get4 <- function(i) {
  if (!is.na(sab.rep.res[i, ]$`4years`)){
    return(0)
  } else{
    getMeanQuantile2(sab.rep.res[i, ]$cpg_list[[1]], sab.rep.res[i, ]$idnum, gset.4) 
  }
}
get8 <- function(i) {
  if (!is.na(sab.rep.res[i, ]$`8years`)){
    return(0)
  } else{
    getMeanQuantile2(sab.rep.res[i, ]$cpg_list[[1]], sab.rep.res[i, ]$idnum, helix_sab) 
  }
}

helix$idnum <- sapply(strsplit(helix$SampleID, "_"), `[`, 2)
helix_sab <- helix[, helix$cohort == "SAB"]
sab.rep.res$rep_quant0 <- sapply(seq_len(nrow(sab.rep.res)), getCord)
sab.rep.res$rep_quant4 <- sapply(seq_len(nrow(sab.rep.res)), get4)
sab.rep.res$rep_quant8 <- sapply(seq_len(nrow(sab.rep.res)), get8)

isSig <- function(x) x < 0.05 | x > 0.95

dats <- c("Cord", "4 years", "8 years")
sab.rep.res2 <- sab.rep.res %>%
 filter(!is.na(rep_quant0) & !is.na(rep_quant4) & !is.na(rep_quant8)) %>%
  filter(idnum != 320) %>%
  mutate(log0 = isSig(rep_quant0),
         log4 = isSig(rep_quant4),
         log8 = isSig(rep_quant8),
         cord_st = ifelse(log0, "Cord blood", ""),
         y4_st = ifelse(log4, "4 years", ""),
         y8_st = ifelse(log8, "8 years", ""),
         sigDatasets = paste(cord_st, y4_st, y8_st, sep = "-"))

## Plot ####
sab.plot <- sab.rep.res2 %>%
  group_by(method, sigDatasets)  %>%
  summarize(n = n()) %>%
  complete(method, n, fill = list(n = 0)) %>%
  mutate(p = n/sum(n), 
         method = factor(method, levels = c("quantile", "beta", "manova", "mlm", "iForest")),
  category = recode(sigDatasets, "--8 years" = "8 years",
  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
  "Cord blood-4 years-" = "Cord blood + 4 years", 
  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
  category = factor(category, levels = c("Cord blood", "4 years", 
                                         "8 years", "Cord blood + 4 years",
                                         "Cord blood + 8 years", "4 years + 8 years",
                                         "Cord blood + 4 years + 8 years"))) %>%
  ggplot(aes(fill = category, color = category, y = p*100, x = method)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Detection time point") +
  scale_color_discrete(name = "Detection time point") +
  scale_y_continuous(name = "Epimutations detected (%)") +
  theme_bw() +
  ggtitle("Biological Replicates") +
  scale_x_discrete(name = "Algorithm") +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
        axis.title = element_text(size = 24),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 24)) 
 


png("figures/poster/sab.rep.png", width = 975, height = 650)
sab.plot
dev.off()



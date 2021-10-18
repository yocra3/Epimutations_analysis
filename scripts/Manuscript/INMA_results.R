#'#################################################################################
#'#################################################################################
#' Figures from INMA results 
#' # Use residuals as main model
#'#################################################################################
#'#################################################################################

# Load data and libraries ####
library(minfi)
library(meffil)
library(tidyverse)
library(cowplot)
library(readxl)
library(robustbase)
library(rexposome)

## Load data ###
load("INMA_comb.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
load("INMA4.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
load("results/epimutations/HELIX.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

load("results/epimutations/INMA_comb.raw.epimutations.allSamples.residuals.Rdata")
load("results/epimutations/INMA4.epimutations.raw.allSamples.residuals.Rdata")
load("results/epimutations/HELIX.epimutations.raw.allSamples.residuals.Rdata")

load("data/HELIX.genexp.Rdata")
load("data/INMA4.genexp.Rdata")

## Exposome data -- Ask permission?
load("data/postExposome.Rdata")
load("data/pregExposome.Rdata")

## Prepare epimutations results ####
### Select quantile, beta and mlm
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
    group_by(method, sample, Sex, idnum, dataset, smoking) %>%
    summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1", 
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))
}

getMeanQuantile <- function(cpgs, idnum, set){
  samps <- colnames(set[, set$idnum == idnum])
  cpgs <- cpgs[cpgs %in% rownames(set)]
  betas <- getBeta(set[cpgs, ])
  quant <- apply(betas, 1, function(x) {
    f <- ecdf(x)
    f(x[colnames(betas) == samps])
  })
  
  mean(quant)
}

getMeanDifference <- function(cpgs, idnum, set){
  samps <- colnames(set[, set$idnum == idnum])
  cpgs <- cpgs[cpgs %in% rownames(set)]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samps, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samps] - means
  mean(diff, na.rm = TRUE)
}

getMeanDifference2 <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}


helixp <- colData(helix)
helixp$idnum <- sapply(strsplit(helixp$SampleID, "_"), `[`, 2)
helixp$Sex <- ifelse(helixp$e3_sex == "female", "F", "M")
helixp$Sample_Name <- helixp$SampleID
colData(helix) <- helixp

res.inma0.df <- make_res_df(res.inma0.residuals.list, pheno = colData(inma0), dataset = "Newborn") 
res.inma0.df$smoking <- ifelse(colData(inma0)[res.inma0.df$sample, "msmk"] == "no smoking", "no", "yes")
inma0.sum.df <- make_sum_df(res.inma0.df)

res.inma4.df <- make_res_df(res.inma4.residuals.list, pheno = colData(inma4), dataset = "4 years") 
res.inma4.df$smoking <- ifelse(colData(inma4)[res.inma4.df$sample, "msmk"] == "no smoking", "no", "yes")
inma4.sum.df <- make_sum_df(res.inma4.df)

## Add smoking data
helix_smk <- read.table("~/data/WS_HELIX/HELIX_analyses/PGRS_smok_GF/db/HELIX_smok.txt", header = TRUE)
helix_smk <- helix_smk[!duplicated(helix_smk$HelixID), ]
helix_smk <- subset(helix_smk, HelixID %in% helix$HelixID)
helix_map <- colData(helix)[, c("HelixID", "SampleID")]
rownames(helix_map) <- helix_map$HelixID
rownames(helix_smk) <- helix_map[helix_smk$HelixID, "SampleID"]

res.helix.df <- make_res_df(res.helix.residuals.list, pheno = helixp, dataset = "8 years") 
res.helix.df$smoking <- ifelse(helix_smk[res.helix.df$sample, "msmok_a"] == 1, "yes", "no")
helix.sum.df <- make_sum_df(res.helix.df)

all.sum.df <- rbind(inma0.sum.df, inma4.sum.df, helix.sum.df) %>%
  mutate(dataset = factor(dataset, levels = c("Newborn", "4 years", "8 years")))

methods <- c("quantile", "beta", "mlm")
names(methods) <- methods

all.res.df <- rbind(res.inma0.df, res.inma4.df, res.helix.df) %>%
  mutate(dataset = factor(dataset, levels = c("Newborn", "4 years", "8 years"))) %>%
  filter(chromosome != 0)

# Descriptives ####
## Proportion of epimutations ####
epi.burden.plot <- all.sum.df %>%
  filter(method %in% methods) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, dataset, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, dataset, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = dataset, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(~ method) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Time point") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample") 
  


png("figures/allINMA.epiburden.png", height = 300, width = 900)
epi.burden.plot
dev.off()

all.sum.df %>%
  filter(method %in% methods) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, dataset, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, dataset, n_cat, fill = list(n = 0, p = 0)) %>% 
  select(-n) %>%
  spread(n_cat, p)

## Epimutation magnitude ####
datasets <- list(Newborn = inma0, `4 years` = inma4, `8 years` = helix)

mean_diffs <- mclapply(seq_len(nrow(all.res.df)), function(i) {
  set <- datasets[[as.character(all.res.df[i, ]$dataset)]]
  getMeanDifference2(all.res.df[i, ]$cpg_ids, all.res.df[i, ]$sample, set = set)
}, mc.cores = 20)
all.res.df$mean_diff <- unlist(mean_diffs)

all.magnitude.plot <- all.res.df %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  ggplot(aes(x = mean_diff, color = method)) +
  geom_density() +
  theme_bw() +
  facet_wrap(~ dataset, ncol = 2) +
  scale_color_discrete(name = "Algorithm") +
  scale_x_continuous("Epimutations magnitude")

png("figures/allINMA.magnitude.png")
all.magnitude.plot
dev.off()

magnitude_tabs <- lapply(unique(all.res.df$dataset), function(dats){
    tab <- all.res.df %>%
      filter(dataset %in% dats) %>%
      mutate(cat = ifelse(mean_diff > 0, "Positive", "Negative")) 
    table(tab$cat, tab$method)
})


magnitude_mod <- lapply(unique(all.res.df$dataset), function(dats){
  tab <- all.res.df %>%
    filter(dataset %in% dats) 
  lm(abs(mean_diff) ~ method, tab)
})
lapply(magnitude_mod, summary)


magnitude_mod_sign <- lapply(unique(all.res.df$dataset), function(dats){
  tab <- all.res.df %>%
    filter(dataset %in% dats) 
  pos <- lm(mean_diff ~ method, tab, mean_diff > 0)
  neg <- lm(mean_diff ~ method, tab, mean_diff < 0)
  list(pos, neg)
})
lapply(unlist(magnitude_mod_sign, recursive = FALSE), summary)
## Age differences
epi_age_tabs <- lapply(methods, function(m){
  tab <- all.sum.df %>%
    filter(method %in% m) %>%
    mutate(cat = ifelse(n == 0, "No epi", "epi")) 
  table(tab$cat, tab$dataset)
})
lapply(epi_age_tabs, chisq.test)

age_mods <- lapply(methods, function(m) {
  tab <- all.sum.df %>%
    filter(method %in% m) %>%
    filter(n > 0) %>%
    mutate(age = ifelse(dataset == "Newborn", 0, ifelse(dataset == "4 years", 4, 8)))
  summary(glmrob(n ~ age, data = tab, family = "poisson"))
})


## Sex differences ####
epi_sex_tabs <- lapply(methods, function(m){
  lapply(list(inma0.sum.df, inma4.sum.df, helix.sum.df), function(dats){
    tab <- dats %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi")) 
    table(tab$cat, tab$Sex)
  })
})
lapply(unlist(epi_sex_tabs, recursive = FALSE), chisq.test)



inma0.sex.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ Sex, data = inma0.sum.df, family = "poisson", subset = method == x & n > 0))
})

inma4.sex.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ Sex, data = inma4.sum.df, family = "poisson", subset = method == x & n > 0))
})

helix.sex.mod <- lapply(methods[-1], function(x) {
  summary(glmrob(n ~ Sex, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})

# 
# sex.burden.plot <- all.sum.df %>%
#   filter(method %in% c("quantile", "beta", "mlm")) %>%
#   mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
#   group_by(method, dataset, Sex, n_cat) %>%
#   summarize(n = n()) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   complete(method, dataset, n_cat, Sex, fill = list(n = 0, p = 0)) %>%
#   ggplot(aes(x = Sex, y = p*100, color = n_cat, fill = n_cat)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   facet_grid(method ~ dataset) +
#   scale_y_continuous(name = "Proportion of individuals") +
#   scale_x_discrete(name = "Sex") +
#   scale_color_discrete(name = "Epimutations per sample") +
#   scale_fill_discrete(name = "Epimutations per sample") 
# 
# png("figures/allINMA.burden.sex.png")
# sex.burden.plot
# dev.off()
# 
# sex.hist.plot <- all.sum.df %>%
#   filter(method %in% c("quantile", "beta", "mlm")) %>%
#   mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
#   ggplot(aes(x = n, fill = Sex)) +
#   geom_histogram(binwidth = 1, position = "dodge") + 
#   theme_bw() +
#   facet_grid(dataset  ~ method, scales = "free_y") +
#   scale_y_continuous(name = "Proportion of individuals") +
#   scale_x_continuous(limits = c(-1, 10)) +
#   scale_color_discrete(name = "Sex") +
#   scale_fill_discrete(name = "Sex") 
# 
# 
# 
# png("figures/allINMA.hist.sex.png")
# sex.hist.plot
# dev.off()


## Smoking differences ####
### Pregnancy smoking ####
epi_smk_tabs <- lapply(methods, function(m){
  lapply(list(inma0.sum.df, inma4.sum.df, helix.sum.df), function(dats){
    tab <- dats %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi"),
             cat = factor(cat, levels = c("No epi", "epi")) )
    table(tab$cat, tab$smoking)
  })
})
lapply(unlist(epi_smk_tabs, recursive = FALSE), fisher.test)

inma0.smk.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ smoking, data = inma0.sum.df, family = "poisson", subset = method == x & n > 0))
})
inma4.smk.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ smoking, data = inma4.sum.df, family = "poisson", subset = method == x & n > 0))
})

helix.smk.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ smoking, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})

### 4 years smoking ####
inma4.sum.df$smk4 <- colData(inma4)[inma4.sum.df$sample, "smk_v4"]

smk4_tabs <- lapply(methods, function(m){
    tab <- inma4.sum.df %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi"),
             cat = factor(cat, levels = c("No epi", "epi")) )
    table(tab$cat, tab$smk4)
})
lapply(smk4_tabs, fisher.test)

inma4.smk4.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ smk4, data = inma4.sum.df, family = "poisson", subset = method == x & n > 0))
})

### 8 years smoking ####
imppost <- toES(imppostnatal, rid = 1)
post_expos <- expos(imppost)
rownames(post_expos) <- pData(imppost)$SampleID
helix.sum.df$smk8 <- post_expos[helix.sum.df$sample, "hs_smk_parents_None"]
helix.sum.df$smk8 <- ifelse(helix.sum.df$smk8 == "neither", "no", "yes")

smk8_tabs <- lapply(methods, function(m){
  tab <- helix.sum.df %>%
    filter(method %in% m) %>%
    mutate(cat = ifelse(n == 0, "No epi", "epi"),
           cat = factor(cat, levels = c("No epi", "epi")) )
  table(tab$cat, tab$smk8)
})
lapply(smk8_tabs, fisher.test)

helix.smk8.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ smk8, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})


smk.burden.plot <- all.sum.df %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  filter(!is.na(smoking)) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, dataset, smoking, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = smoking, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ dataset) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "smoking") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample") 



png("figures/allINMA.burden.smoking.png")
smk.burden.plot
dev.off()

## Batch differences in INMA 0 ####
inma0.sum.df$batch <- colData(inma0)[inma0.sum.df$sample, "Batch"]
batch_tabs <- lapply(methods, function(m){
    tab <- inma0.sum.df %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi")) 
    table(tab$cat, tab$batch)
  })
lapply(batch_tabs, chisq.test)
lapply(batch_tabs, fisher.test)

batch_tabs2 <- lapply(methods, function(m){
  tab <- inma0.sum.df %>%
    filter(method %in% m) 
  table(tab$n_cat, tab$batch)
})

inma0.batch.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ batch, data = inma0.sum.df, family = "poisson", subset = method == x & n > 0))
})

# ### Load esteller results ####
# load("results/epimutations/INMA_comb.epimutations.esteller.residuals.Rdata")
# res.esteller.df <- make_res_df(res.inma0.esteller.residuals.list, pheno = colData(inma0), dataset = "Newborn") 
# res.esteller.df$smoking <- ifelse(colData(inma0)[res.esteller.df$sample, "msmk"] == "no smoking", "no", "yes")
# esteller.sum.df <- make_sum_df(res.esteller.df)
# 
# esteller.sum.df$id <- paste(esteller.sum.df$sample, esteller.sum.df$method)
# 
# comb.esteller <- inma0.sum.df %>%
#   ungroup() %>%
#   mutate(id = paste(sample, method)) %>%
#   select(method, sample, n) %>%
#   right_join(ungroup(esteller.sum.df) %>% select(method, sample, n), by = c("method", "sample"))
# table(all = comb.esteller$n.x, esteller = comb.esteller$n.y, comb.esteller$method)
# 
# ### Load medall results ####
# load("results/epimutations/INMA_comb.epimutations.medall.residuals.Rdata")
# res.medall.df <- make_res_df(res.inma0.medall.residuals.list, pheno = colData(inma0), dataset = "Newborn") 
# res.medall.df$smoking <- ifelse(colData(inma0)[res.medall.df$sample, "msmk"] == "no smoking", "no", "yes")
# medall.sum.df <- make_sum_df(res.medall.df)
# 
# medall.sum.df$id <- paste(medall.sum.df$sample, medall.sum.df$method)
# 
# comb.medall <- inma0.sum.df %>%
#   ungroup() %>%
#   mutate(id = paste(sample, method)) %>%
#   select(method, sample, n) %>%
#   right_join(ungroup(medall.sum.df) %>% select(method, sample, n), by = c("method", "sample"))
# table(all = comb.medall$n.x, medall = comb.medall$n.y, comb.medall$method)

## Cohort differences in HELIX ####
helix.sum.df$cohort <- colData(helix)[helix.sum.df$sample, "cohort"]
cohort_tabs <- lapply(methods, function(m){
  tab <- helix.sum.df %>%
    filter(method %in% m) %>%
    mutate(cat = ifelse(n == 0, "No epi", "epi")) 
  table(tab$cat, tab$cohort)
})
lapply(cohort_tabs, fisher.test)

cohort.burden.plot <- helix.sum.df %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, cohort, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, n_cat, cohort, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = cohort, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ .) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Cohort") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/HELIX.burden.cohort.png")
cohort.burden.plot
dev.off()

helix.cohort.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ cohort, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})


# Replicability across time points ####
## Sabadell results
sab.res.df <- rbind(res.inma0.df, res.inma4.df, res.helix.df) %>%
  mutate(dataset = factor(dataset, levels = c("Newborn", "4 years", "8 years"))) %>%
  filter(dataset != "8 years" | grepl("SAB", sample))


selSamps <- sab.res.df %>%
  select(dataset, idnum) %>%
  distinct() %>%
  group_by(idnum) %>%
  summarize(n = n()) %>%
  filter(n == 3)

sab.rep.res <- sab.res.df %>%
  filter(idnum %in% selSamps$idnum) %>%
  filter(chromosome != 0) %>%
  group_by(dataset, method, idnum, epi_region_id) %>%
  summarize(cpg_ids = paste(unique(cpg_ids), collapse = ",")) %>%
  spread(dataset, cpg_ids)

sab.rep.res$cpg_list = lapply(seq_len(nrow(sab.rep.res)), function(i) {
  cpg_string <- paste(sab.rep.res[i, c("Newborn", "4 years", "8 years")], collapse = ",")
  cpg_list <- strsplit(cpg_string, ",")[[1]]
  unique(cpg_list[cpg_list != "NA"])
})
helix_sab <- helix[, helix$cohort == "SAB"]

getRep <- function(i, col, set) {
  if (!is.na(sab.rep.res[i, col])){
    return(0)
  } else{
    getMeanQuantile(sab.rep.res[i, ]$cpg_list[[1]], sab.rep.res[i, ]$idnum, set) 
  }
}


sab.rep.res$rep_quant0 <- sapply(seq_len(nrow(sab.rep.res)), getRep, col = "Newborn", set = inma0)
sab.rep.res$rep_quant4 <- sapply(seq_len(nrow(sab.rep.res)), getRep, col = "4 years", set = inma4)
sab.rep.res$rep_quant8 <- sapply(seq_len(nrow(sab.rep.res)), getRep, col = "8 years", set = helix_sab)


isSig <- function(x) x < 0.05 | x > 0.95

sab.rep.res2 <- sab.rep.res %>%
  ungroup() %>%
  filter(!is.na(rep_quant0) & !is.na(rep_quant4) & !is.na(rep_quant8)) %>%
  mutate(cord_st = ifelse(!is.na(Newborn), "Cord blood", ""),
         y4_st = ifelse(!is.na(`4 years`), "4 years", ""),
         y8_st = ifelse(!is.na(`8 years`), "8 years", ""),
         sigDatasets = paste(cord_st, y4_st, y8_st, sep = "-"),
         log0 = isSig(rep_quant0),
         log4 = isSig(rep_quant4),
         log8 = isSig(rep_quant8),
         cord_st2 = ifelse(log0, "Cord blood", ""),
         y4_st2 = ifelse(log4, "4 years", ""),
         y8_st2 = ifelse(log8, "8 years", ""),
         sigDatasets2 = paste(cord_st2, y4_st2, y8_st2, sep = "-"))

## Plot 
sab.replic.strict.plot <- sab.rep.res2 %>%
  group_by(method, sigDatasets)  %>%
  filter(method %in% methods) %>%
  summarize(n = n()) %>%
  complete(method, n, fill = list(n = 0)) %>%
  mutate(p = n/sum(n), 
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
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
  ggtitle("Epimutations called") +
  scale_x_discrete(name = "Algorithm") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

sab.replic.signal.plot <- sab.rep.res2 %>%
  group_by(method, sigDatasets2)  %>%
  filter(method %in% methods) %>%
  summarize(n = n()) %>%
  complete(method, n, fill = list(n = 0)) %>%
  mutate(p = n/sum(n), 
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets2, "--8 years" = "8 years",
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
  ggtitle("Epimutations with signal") +
  scale_x_discrete(name = "Algorithm") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

png("figures/INMAsab.replic.png", width = 975, height = 400)
plot_grid(sab.replic.strict.plot, sab.replic.signal.plot, nrow = 1, 
          rel_widths = c(2, 3))
dev.off()

sab.rep.res2 %>%
  group_by(method, sigDatasets)  %>%
  filter(method %in% methods) %>%
  summarize(n = n()) %>%
  complete(method, n, fill = list(n = 0)) %>%
  mutate(p = n/sum(n), 
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
                                  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
                                  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
                                  "Cord blood-4 years-" = "Cord blood + 4 years",
                                  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
         category = factor(category, levels = c("Cord blood", "4 years",
                                                "8 years", "Cord blood + 4 years",
                                                "Cord blood + 8 years", "4 years + 8 years",
                                                "Cord blood + 4 years + 8 years"))) %>%
  select(method, p, category) %>%
  spread(method, p)



sab.rep.res2 %>%
  group_by(method, sigDatasets2)  %>%
  filter(method %in% methods) %>%
  summarize(n = n()) %>%
  complete(method, n, fill = list(n = 0)) %>%
  mutate(p = n/sum(n), 
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets2, "--8 years" = "8 years",
                                  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
                                  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
                                  "Cord blood-4 years-" = "Cord blood + 4 years",
                                  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
         category = factor(category, levels = c("Cord blood", "4 years",
                                                "8 years", "Cord blood + 4 years",
                                                "Cord blood + 8 years", "4 years + 8 years",
                                                "Cord blood + 4 years + 8 years"))) %>%
  select(method, p, category) %>%
  spread(method, p)

## Check effect across time points ####
sab.rep.res2.meth <-  sab.rep.res2 %>%
  filter(method %in% methods)
sab.rep.res2.meth$mean_diff0 <- sapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
  getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = inma0)
})
sab.rep.res2.meth$mean_diff4 <- sapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
  getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = inma4)
})
sab.rep.res2.meth$mean_diff8 <- sapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
  getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = helix_sab)
})

sab.replic.diff.plot <- sab.rep.res2.meth %>% 
  group_by(method, sigDatasets2) %>%
  slice_sample(n = 20) %>% 
  ungroup() %>%
  gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
  mutate(age = gsub("mean_diff", "", age),
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets2, "--8 years" = "8 years",
                                  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
                                  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
                                  "Cord blood-4 years-" = "Cord blood + 4 years", 
                                  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
         category = factor(category, levels = c("Cord blood", "4 years", 
                                                "8 years", "Cord blood + 4 years",
                                                "Cord blood + 8 years", "4 years + 8 years",
                                                "Cord blood + 4 years + 8 years"))) %>%
  ggplot(aes(x = age, y = mean_diff, group = epi_region_id, color = category)) +
  geom_line() + geom_point() +
  facet_wrap(~ method) +
  scale_y_continuous("Epimutation mean methylation difference") +
  scale_color_discrete(name = "Epimutation signal detection") +
  scale_x_discrete("Age (in years)") +
  theme_bw()

png("figures/INMAsab.replic.diff.png")
sab.replic.diff.plot
dev.off()


sab.replic.diff.coher.plot <- sab.rep.res2.meth %>% 
  filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
  group_by(method, sigDatasets) %>%
  slice_sample(n = 20) %>% 
  gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
  mutate(age = gsub("mean_diff", "", age),
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
                                  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
                                  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
                                  "Cord blood-4 years-" = "Cord blood + 4 years", 
                                  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
         category = factor(category, levels = c("Cord blood", "4 years", 
                                                "8 years", "Cord blood + 4 years",
                                                "Cord blood + 8 years", "4 years + 8 years",
                                                "Cord blood + 4 years + 8 years"))) %>%
  ggplot(aes(x = age, y = mean_diff, group = paste(epi_region_id, idnum), color = category)) +
  geom_line() + geom_point() +
  facet_wrap(~ method) +
  scale_y_continuous("Epimutation mean methylation difference") +
  scale_color_discrete(name = "Epimutation detection") +
  scale_x_discrete("Age (in years)") +
  theme_bw()

png("figures/INMAsab.replic.diff.coherent.png")
sab.replic.diff.coher.plot
dev.off()

replic.direc.mod <- lapply(methods, function(x) {
  df <- sab.rep.res2.meth %>%
    filter(method == x) %>%
    filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
    mutate(id = paste(idnum, epi_region_id)) %>%
    gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
    mutate(age = gsub("mean_diff", "", age),
           age = as.numeric(age))
  lmer(abs(mean_diff) ~ age + (1|id), df)
})
lapply(replic.direc.mod, car::Anova)



sab.replic.diff.coher.plot2 <- sab.rep.res2.meth %>% 
  filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
  group_by(method, sigDatasets) %>%
  gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
  mutate(age = gsub("mean_diff", "", age),
         method = factor(method, levels = c("quantile", "beta", "mlm")),
         category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
                                  "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
                                  "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
                                  "Cord blood-4 years-" = "Cord blood + 4 years", 
                                  "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
         category = factor(category, levels = c("Cord blood", "4 years", 
                                                "8 years", "Cord blood + 4 years",
                                                "Cord blood + 8 years", "4 years + 8 years",
                                                "Cord blood + 4 years + 8 years"))) %>%
  ggplot(aes(x = age, y = abs(mean_diff))) +
geom_boxplot() +  facet_wrap(~ method) +
  scale_y_continuous("Epimutation mean methylation difference") +
  scale_color_discrete(name = "Epimutation detection") +
  scale_x_discrete("Age (in years)") +
  theme_bw()

png("figures/INMAsab.replic.diff.coherent2.png")
sab.replic.diff.coher.plot2
dev.off()


### Check batch in INMA0 ####
batch_map <- colData(inma0)[, c("idnum", "Batch")]
rownames(batch_map) <- batch_map$idnum
quant.rep <- sab.rep.res2.meth %>% 
  filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
  filter(method == "quantile") %>%
  mutate(batch = batch_map[as.character(idnum), "Batch"])

png("figures/INMAsab.replic.diff.quantile.png")
quant.rep %>%
  gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
  mutate(age = gsub("mean_diff", "", age)) %>%
  ggplot(aes(x = age, y = mean_diff, group = paste(epi_region_id, idnum), color = batch)) +
  geom_line() + geom_point() +
  scale_y_continuous("Epimutation mean methylation difference") +
  scale_color_discrete(name = "Epimutation detection") +
  scale_x_discrete("Age (in years)") +
  theme_bw()
dev.off()



# Correlation with gene expression ####

getGenesZ <- function(epi_df, row, gexpZ, gexpRange, window){
  
  exp_name <- epi_df[row, "SampleName"]
  if (!exp_name %in% colnames(gexpZ)){
    return(  list(geneName = NA,
                  Zscore = NA,
                  distance = NA))
  }
  
  epi_range <- GRanges(gsub("_", ":", epi_df[row, "epi_region_id"]))
  
  transRanges <- subsetByOverlaps(gexpRange, epi_range + window)
  if (length(transRanges) == 0){
    return(  list(geneName = NA,
                  Zscore = NA,
                  distance = NA))
  }
  
  selTranscripts <- names(transRanges)
  
  zs <- gexpZ[selTranscripts, exp_name]
  list(geneNames = selTranscripts,
       Zscores = zs,
       distances = distance(epi_range, transRanges))
  
}

selectTop <- function(x) {
  if (length(x$Zscore) == 1){
    data.frame(geneName = x$geneName, Zscore = x$Zscore, distance = x$Zscore) 
  }
  else {
    ind <- which.max(abs(x$Zscore))
    data.frame(geneName =  x$geneName[ind], Zscore = x$Zscore[ind], distance = x$distance[ind])
  }
}
selectNearest <- function(x) {
  if (length(x$Zscore) == 1){
    data.frame(geneName = x$geneName, Zscore = x$Zscore, distance = x$Zscore) 
  }
  else {
    ind <- which.min(x$distance)
    data.frame(geneName =  x$geneName[ind], Zscore = x$Zscore[ind], distance = x$distance[ind])
  }
}

## INMA4 ####
inma4.exp.z <- t(apply(exprs(inma.expset), 1, function(x) (x - mean(x))/sd(x)))
inma4.exp.range <- makeGRangesFromDataFrame(fData(inma.expset))

res.inma4.filt <- res.inma4.df %>%
  filter(chromosome != 0 & method %in% methods)

res.inma4.filt$SampleName <- paste0("X04_", res.inma4.filt$idnum)
geneInfo4 <- lapply(seq_len(nrow(res.inma4.filt)), getGenesZ, epi_df = res.inma4.filt,
                    gexpZ = inma4.exp.z, gexpRange = inma4.exp.range, window = 250e3)


geneInfo4df <- Reduce(rbind, lapply(geneInfo4, selectTop))
res.inma4.filt <- cbind(res.inma4.filt, geneInfo4df)

geneInfo4Near <- Reduce(rbind, lapply(geneInfo4, selectNearest))
colnames(geneInfo4Near) <- paste0(colnames(geneInfo4Near), "Near" )
res.inma4.filt <- cbind(res.inma4.filt, geneInfo4Near)

inma4.gexp.plot <- res.inma4.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
  ggplot(aes(color = method, x = Zscore)) +
  geom_density() +
  theme_bw() +
  facet_grid(n_cat ~ .)

png("figures/INMA4.genexp.png")
inma4.gexp.plot
dev.off()

res.inma4.filt$pvals <- pnorm(abs(res.inma4.filt$Zscore), lower = FALSE)
inma4.gexp.pvals.plot <- res.inma4.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
  ggplot(aes(x = method, y = -log10(pvals))) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ n_cat)

png("figures/INMA4.genexp.pvals.png")
inma4.gexp.pvals.plot
dev.off()


cot <- res.inma4.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))

summary(lm(abs(Zscore) ~ n_cat, cot, subset = method == "beta"))  
summary(lm(abs(Zscore) ~ n_cat, cot, subset = method == "mlm"))  


cot <- res.inma4.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(ZscoreNear)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))

summary(lm(abs(ZscoreNear) ~ n_cat, cot, subset = method == "beta"))  
summary(lm(abs(ZscoreNear) ~ n_cat, cot, subset = method == "mlm"))  


## HELIX ####
helix.gexp <- transcriptome_subcohort_f1[, transcriptome_subcohort_f1$cohort != "KAN" & 
                                           transcriptome_subcohort_f1$h_ethnicity_3cat == "WhiteEur_WhiteOther"] 

helix.exp.z <- t(apply(exprs(helix.gexp), 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
helix.exp.range <- makeGRangesFromDataFrame(fData(helix.gexp))

res.helix.filt <- res.helix.df %>%
  filter(chromosome != 0 & method %in% methods)

res.helix.filt$SampleName <- res.helix.filt$sample
geneInfoHelix <- lapply(seq_len(nrow(res.helix.filt)), getGenesZ, epi_df = data.frame(res.helix.filt),
                    gexpZ = helix.exp.z, gexpRange = helix.exp.range, window = 250e3)


geneInfohelixdf <- Reduce(rbind, lapply(geneInfoHelix, selectTop))
res.helix.filt <- cbind(res.helix.filt, geneInfohelixdf)

geneInfohelixNear <- Reduce(rbind, lapply(geneInfoHelix, selectNearest))
colnames(geneInfohelixNear) <- paste0(colnames(geneInfohelixNear), "Near" )
res.helix.filt <- cbind(res.helix.filt, geneInfohelixNear)

helix.gexp.plot <- res.helix.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
  ggplot(aes(color = method, x = Zscore)) +
  geom_density() +
  theme_bw() +
  facet_grid(n_cat ~ .)

png("figures/HELIX.genexp.png")
helix.gexp.plot
dev.off()

res.helix.filt$pvals <- pnorm(abs(res.helix.filt$Zscore), lower = FALSE)
helix.gexp.pvals.plot <- res.helix.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
  ggplot(aes(x = method, y = -log10(pvals))) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ n_cat)

png("figures/HELIX.genexp.pvals.png")
helix.gexp.pvals.plot
dev.off()


pac <- res.helix.filt %>% 
  group_by(method, idnum) %>%
  filter(!is.na(Zscore)) %>%
  mutate(n = n(),
         n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))

summary(lm(abs(Zscore) ~ n_cat, pac, subset = method == "beta"))  
summary(lm(abs(Zscore) ~ n_cat, pac, subset = method == "mlm"))  

## Replicates ####
sab.rep.gexp <- sab.rep.res2 %>%
  filter(method %in% methods)

sab.rep.gexp$SampleName <-  paste0("X04_", sab.rep.gexp$idnum)
geneInfoSab4 <- lapply(seq_len(nrow(sab.rep.gexp)), getGenesZ, epi_df = data.frame(sab.rep.gexp),
                   gexpZ = inma4.exp.z, gexpRange = inma4.exp.range, window = 250e3)

sab.rep.gexp$SampleName <-  colnames(helix_sab)[match(sab.rep.gexp$idnum, helix_sab$idnum)]
geneInfoSab8 <- lapply(seq_len(nrow(sab.rep.gexp)), getGenesZ, epi_df = data.frame(sab.rep.gexp),
                       gexpZ = helix.exp.z, gexpRange = helix.exp.range, window = 250e3)

geneInfoSab4top <- Reduce(rbind, lapply(geneInfoSab4, selectTop))
colnames(geneInfoSab4top) <- paste0(colnames(geneInfoSab4top), "Top4" )
sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab4top)

geneInfoSab4Near <- Reduce(rbind, lapply(geneInfoSab4, selectNearest))
colnames(geneInfoSab4Near) <- paste0(colnames(geneInfoSab4Near), "Near4" )
sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab4Near)

geneInfoSab8top <- Reduce(rbind, lapply(geneInfoSab8, selectTop))
colnames(geneInfoSab8top) <- paste0(colnames(geneInfoSab8top), "Top8" )
sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab8top)

geneInfoSab8Near <- Reduce(rbind, lapply(geneInfoSab8, selectNearest))
colnames(geneInfoSab8Near) <- paste0(colnames(geneInfoSab8Near), "Near8" )
sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab8Near)

sab.gexp4.top.plot <- sab.rep.gexp %>% 
  ggplot(aes(fill = sigDatasets2, x = ZscoreTop4)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y") +
  theme(legend.position = "none")

sab.gexp8.top.plot <- sab.rep.gexp %>% 
  ggplot(aes(fill = sigDatasets2, x = ZscoreTop8)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y")

png("figures/INMAsab.top.genexp.png", width = 1200)
plot_grid(sab.gexp4.top.plot, sab.gexp8.top.plot, nrow = 1, rel_widths = c(3, 5))
dev.off()



sab.gexp4.near.plot <- sab.rep.gexp %>% 
  ggplot(aes(fill = sigDatasets2, x = ZscoreNear4)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y") +
  theme(legend.position = "none")

sab.gexp8.near.plot <- sab.rep.gexp %>% 
  ggplot(aes(fill = sigDatasets2, x = ZscoreNear8)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y")

png("figures/INMAsab.near.genexp.png", width = 1200)
plot_grid(sab.gexp4.near.plot, sab.gexp8.near.plot, nrow = 1, rel_widths = c(3, 5))
dev.off()


sab.rep.gexp %>%
  filter(!grepl("4 years", sigDatasets2)) %>%
  pull(ZscoreTop4) %>%
  abs() %>%
  quantile(., 0.95, na.rm = TRUE)


sab.rep.gexp %>%
  filter(!grepl("8 years", sigDatasets2)) %>%
  pull(ZscoreTop8) %>%
  abs() %>%
  quantile(., 0.95, na.rm = TRUE)


sab.rep.gexp %>%
  filter(!grepl("4 years", sigDatasets2)) %>%
  pull(ZscoreNear4) %>%
  abs() %>%
  quantile(., 0.95, na.rm = TRUE)


sab.rep.gexp %>%
  filter(!grepl("8 years", sigDatasets2)) %>%
  pull(ZscoreNear8) %>%
  abs() %>%
  quantile(., 0.95, na.rm = TRUE)

sab.gexp4.distance.signal.plot <- sab.rep.gexp %>% 
  ggplot(aes(x = unlist(distance4), color = sigDatasets2, y = Zscore4)) +
  geom_point() +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y")

png("figures/INMAsab.distance.signal.genexp.png")
sab.gexp4.distance.signal.plot
dev.off()




png("figures/INMAsab.distanceDistr.signal.genexp.png")

sab.rep.gexp %>% 
  ggplot(aes(y = unlist(distance4), x = sigDatasets2)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(method ~ ., scales = "free_y")
dev.off()



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
library(epimutacions)
library(writexl)

## Load data ###
load("INMA0combined.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
load("INMA4.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
load("results/epimutations/HELIX.normalizedRaw.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

load("results/epimutations/INMA0combined.raw.epimutations.allSamples.residuals.Rdata")
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

helixp$smoking <- helix_smk[rownames(helixp), "msmok_a"]

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
## Data descriptives ####
getSum <- function(vec, type = "continuous"){
  if (type == "continuous"){
    c(median = median(vec, na.rm = TRUE), 
      range = quantile(vec, probs = c(0.25, 0.75), na.rm = TRUE))
  } else if (type == "categorical"){
    t <- table(vec, useNA = "ifany")
    data.frame(names = names(t),  tab = as.vector(t), 
               props = as.vector(prop.table(t)))
  }
}
getSumCohort <- function(tab, varname, type = "continuous"){
  tab$var <- tab[, varname]
  df <- tab %>%
    data.frame() %>%
    group_by(cohort)
  if (type == "continuous"){
    df %>% 
      summarize(median = median(var, na.rm = TRUE), 
                range.25 = quantile(var, probs = 0.25, na.rm = TRUE),
                range.75 = quantile(var, probs = 0.75, na.rm = TRUE)) %>%
      mutate(val = sprintf("%.2f (%.2f-%.2f)", median, range.25, range.75)) %>%
      select(cohort, val)
    
  } else if (type == "categorical"){
    df %>%
      group_by(cohort, var) %>%
      summarize(N = n()) %>%
      group_by(cohort) %>%
      mutate(val = paste0(N, " (", round(N/sum(N)*100, 1), "%)")) %>%
      select(cohort, var, val) %>%
      spread(var, val)
  }
}
catVecs0 <- lapply(c("Sex", "Batch", "msmk"), function(x) colData(inma0)[, x])
catSums0 <- lapply(catVecs0, getSum, type = "categorical")
tab0 <- Reduce(rbind, catSums0) %>%
  mutate(val = sprintf("%i (%.2f%%)", tab, props*100)) %>%
  select(names, val) %>%
rbind(data.frame(names = "Total", val = ncol(inma0)), .)
  
catVecs4 <- lapply(c("Sex", "msmk"), function(x) colData(inma4)[, x])
catSums4 <- lapply(catVecs4, getSum, type = "categorical")

contSums4 <- getSum(inma4$age_v4, type = "continuous")
tab4 <- Reduce(rbind, catSums4) %>%
  mutate(val = sprintf("%i (%.2f%%)", tab, props*100)) %>%
  select(names, val) %>%
  rbind(data.frame(names = "Total", val = ncol(inma4)), .) %>%
  rbind(data.frame(names = "age", val = sprintf("%.2f (%.2f-%.2f)", contSums4[1], contSums4[2], contSums4[3])))

helixp$smoking <- ifelse(helixp$smoking == 0, "no smoking", "yes, any smoking during pregnancy")
catVecs8 <- lapply(c("Sex", "smoking"), function(x) helixp[, x])
catSums8 <- lapply(catVecs8, getSum, type = "categorical")

contSums8 <- getSum(helixp$age_sample_years, type = "continuous")
tab8 <- Reduce(rbind, catSums8) %>%
  mutate(val = sprintf("%i (%.2f%%)", tab, props*100)) %>%
  select(names, val) %>%
  rbind(data.frame(names = "Total", val = nrow(helixp)), .) %>%
  rbind(data.frame(names = "age", val = sprintf("%.2f (%.2f-%.2f)", contSums8[1], contSums8[2], contSums8[3])))


catSumscohort <- lapply(c("Sex", "smoking"), getSumCohort, tab = helixp, type = "categorical") %>%
  Reduce(cbind, .) %>%
  t() %>%
  rbind(sapply(unique(helixp$cohort), function(x) sum(helixp$cohort == x)))

colnames(catSumscohort) <- catSumscohort[1, ]
catSumscohort <- catSumscohort[-c(1, 4), ]

cohortCont <- getSumCohort(tab = helixp, "age_sample_years") %>%
  t()
colnames(cohortCont) <- cohortCont[1, ]
cohortCont <- cohortCont[-1, , drop = FALSE]
tabCohort <- rbind(catSumscohort, cohortCont) %>%
  data.frame()
tabCohort$names <- rownames(tabCohort)
tabCohort$names[5] <- NA
tabCohort$names[6] <- "Total" 
tabCohort$names[7] <- "age" 

full.tab <- Reduce(function(x, y) full_join(x, y, by = "names"),
                   list(tab0, tab4, tab8, tabCohort))
full.tab$names[is.na(full.tab$names)] <- "Missing smoking"
rownames(full.tab) <- full.tab$names

write.table(full.tab[c(1:3, 9, 6:8, 4:5), ], file = "figures/PopDescrip_cat.txt", col.names = TRUE, 
            quote = FALSE, row.names = FALSE, sep = "\t")


## Check differences in variables ####
summary(lm(age_sample_years ~ cohort, helixp))
chisq.test(table(helixp$smoking, helixp$cohort))

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
  facet_grid(~ dataset) +
  scale_color_discrete(name = "Algorithm") +
  scale_x_continuous("Epimutations magnitude")

png("figures/allINMA.magnitude.png", height = 200)
all.magnitude.plot
dev.off()


# all.magnitude.plot2 <- all.res.df %>%
#   mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
#   ggplot(aes(x = method, y = mean_diff, color = method)) +
#   geom_violin() +
#   geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5,
#                binwidth = 0.01) +
#   theme_bw() +
#   facet_grid(~ dataset) +
#   scale_color_discrete(name = "Algorithm") +
#   scale_y_continuous("Epimutations magnitude")
# 
# png("figures/allINMA.magnitude2.png")
# all.magnitude.plot2
# dev.off()
# 




magnitude_tabs <- lapply(unique(all.res.df$dataset), function(dats){
    tab <- all.res.df %>%
      filter(dataset %in% dats) %>%
      mutate(cat = ifelse(mean_diff > 0, "Positive", "Negative")) 
    table(tab$cat, tab$method)
})
lapply(magnitude_tabs, prop.table, margin = 2)

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


## Recurrent epimutations ####
### Load literature regions
epi_lit <- read_excel("data/Epimutations.PMID32937144.xlsx", skip = 2)
epi_litGR <- makeGRangesFromDataFrame(epi_lit, seqnames.field = "Chr, DMR",
                                      start.field = "Start, DMR (hg19)",
                                      end.field = "End, DMR (hg19)",
                                      keep.extra.columns = TRUE)

## Children and European
### Selected cohorts: GSE105018, GSE82273, GSE103657, EGAD00010001461
child_cohorts <- c("Environmental risk Twins_GEO GSE105018", "Infant_GEO GSE82273", "Neonatal_GEO GSE103657", "BAMSE_EGA D00010001461")
total_samps <- 3262

epi_litGR$n_inds <- rowSums(data.matrix(mcols(epi_litGR)[, child_cohorts]))
epi_litGR$litfreq <- epi_litGR$n_inds/total_samps
epi_litGR$n_cohort <- rowSums(data.matrix(mcols(epi_litGR)[, child_cohorts]) > 0)

# epi_litGR$n_cohort <- sapply(strsplit(epi_litGR$Cohort, ","), length)
# epi_litGR$min_freq <- sapply(strsplit(epi_litGR$`Frequency of DMR per 10,000 (95% CI)`, "-"), `[`, 1)
# epi_litGR$max_freq <- sapply(strsplit(epi_litGR$`Frequency of DMR per 10,000 (95% CI)`, "-"), `[`, 2)
# epi_litGR$min_freq <- as.numeric(epi_litGR$min_freq)
# epi_litGR$max_freq <- as.numeric(epi_litGR$max_freq)

over <- findOverlaps(epi_litGR, candRegsGR)
epi_litGR$epi_region_id <- names(candRegsGR)[to(over)]

### Distribution ####
recur.epi <- all.res.df %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(dataset) %>%
  mutate(dataset_n = length(unique(sample))) %>%
  group_by(method, dataset, epi_region_id, dataset_n) %>%
  summarize(n = n(),
            freq = n()/dataset_n) %>%
  ungroup() %>%
  distinct() %>%
  mutate(type = "All samples")

recur.epi.freq.plot <- recur.epi %>%
  ggplot(aes(x = freq*100)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(method ~ dataset) +
  scale_x_log10(name = "Epimutation frequency (%)", 
                breaks = c(0.1, 0.2, 0.5, 1, 2, 5))

png("figures/allINMA.epi.freq.png")
recur.epi.freq.plot
dev.off()

### Correlation with literature ####
recur.epi.comb <- left_join(recur.epi, 
                            select(data.frame(mcols(epi_litGR)), n_cohort, litfreq, epi_region_id),
                            by = "epi_region_id") %>%
  mutate(n_cohort = ifelse(is.na(n_cohort), 0, n_cohort),
         litfreq = ifelse(is.na(litfreq) | litfreq == 0, 3e-4, litfreq)) %>%
  group_by(method, dataset, epi_region_id, freq) %>%
  summarize(n_cohort = max(n_cohort, na.rm = TRUE),
            litfreq = max(litfreq, na.rm = TRUE))


# recur.epi.freq.cor.min.plot <- recur.epi.comb %>%
#   ggplot(aes(x = freq*100, y = min_freq/100)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   facet_grid(method ~ dataset) +
#   scale_x_log10(name = "Epimutation frequency (%)", 
#                 breaks = c(0.1, 0.2, 0.5, 1, 2, 5)) +
#   scale_y_log10(name = "Epimutation frequency (%)", 
#               breaks = c(0.1, 0.2, 0.5, 1, 2, 5))

recur.epi.freq.shared.plot <- recur.epi.comb %>%
  group_by(method, dataset) %>%
  summarize(p = mean(n_cohort > 0)) %>%
  ggplot(aes(x = dataset, y = p*100)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(. ~ method) +
  scale_x_discrete(name = "Dataset") +
  scale_y_continuous(name = "Epimutations in literature (%)", limits = c(0, 100))

png("figures/allINMA.epi.freq.shared.png", height = 300)
recur.epi.freq.shared.plot
dev.off()


recur.epi.freq.cor.plot <- recur.epi.comb %>%
  filter(n_cohort > 0) %>%
  ggplot(aes(x = freq*100, y = litfreq*100)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(method ~ dataset) +
  scale_x_log10(name = "Epimutation frequency (%)", 
                breaks = c(0.1, 0.2, 0.5, 1, 2, 5)) +
  scale_y_log10(name = "Literature frequency (%)", 
                breaks = c(0.1, 0.2, 0.5, 1, 2, 5))


png("figures/allINMA.epi.freq.cor.png")
recur.epi.freq.cor.plot
dev.off()

epi_rec_shared <- lapply(methods, function(m){
    tab <- recur.epi.comb %>%
      filter(method %in% m) %>%
      filter(!is.na(litfreq))
    table(tab$dataset, tab$n_cohort > 0)
})


epi_rec_cors <- lapply(methods, function(m){
  lapply(unique(recur.epi.comb$dataset), function(dats){
    tab <- recur.epi.comb %>%
      filter(method %in% m) %>%
      filter(dataset %in% dats) %>%
      filter(n_cohort > 0)
    cor.test(log10(tab$freq), log10(tab$litfreq))
  })
})
sapply(unlist(epi_rec_cors, recursive = FALSE), function(x) x$estimate)
sapply(unlist(epi_rec_cors, recursive = FALSE), function(x) x$p.value)




epi_rec_list <- spread(recur.epi.comb, dataset, freq) %>%
  replace(is.na(.), 0) 

epi_rec_list <- lapply(methods, function(x) {
  filter(epi_rec_list, method == x) %>%
    arrange(desc(Newborn))
})
writexl::write_xlsx(epi_rec_list, path = "figures/INMA.Epimutations.xlsx")


## Age differences ####
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

inma4.sum.df$age <- colData(inma4)[inma4.sum.df$sample, "age_v4"]

age_mods4 <- lapply(methods, function(m) {
  tab <- inma4.sum.df %>%
    subset(method == m) %>%
    mutate(epi = ifelse(n == 0, 0, 1))
  summary(glm(epi ~ age, data = tab, family = "binomial"))
})

age_mods4_e <- lapply(methods, function(m) {
  summary(glmrob(n ~ age, data = subset(inma4.sum.df, method == m & n > 0), family = "poisson"))
})



helix.sum.df$age <- colData(helix)[helix.sum.df$sample, "age_sample_years"]

age_mods8 <- lapply(methods, function(m) {
  tab <- helix.sum.df %>%
    subset(method == m) %>%
    mutate(epi = ifelse(n == 0, 0, 1))
  summary(glm(epi ~ age, data = tab, family = "binomial"))
})

age_mods8_adj <- lapply(methods, function(m) {
  tab <- helix.sum.df %>%
    subset(method == m) %>%
    mutate(epi = ifelse(n == 0, 0, 1))
  summary(glm(epi ~ age + cohort, data = tab, family = "binomial"))
})

age_mods8_e <- lapply(methods, function(m) {
  summary(glmrob(n ~ age, data = subset(helix.sum.df, method == m & n > 0), family = "poisson"))
})

age_mods8_e_adj <- lapply(methods, function(m) {
  summary(glmrob(n ~ age + cohort, data = subset(helix.sum.df, method == m & n > 0), family = "poisson"))
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

helix.sex.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ Sex, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})

helix.sex.mod.b <- lapply(methods, function(x) {
  tab <- subset(helix.sum.df, method == x) %>%
    mutate(epi = ifelse(n > 0, 1, 0))
  summary(glm(epi ~ Sex, data = tab, family = "binomial"))
})

helix.sex.mod.b.adj <- lapply(methods, function(x) {
  tab <- subset(helix.sum.df, method == x) %>%
    mutate(epi = ifelse(n > 0, 1, 0))
  summary(glm(epi ~ Sex + cohort, data = tab, family = "binomial"))
})


sex.burden.plot <- all.sum.df %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm")),
         Sex = ifelse(Sex == "F", "Girls", "Boys")) %>%
  group_by(method, dataset, Sex, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, dataset, n_cat, Sex, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = Sex, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(dataset ~ method) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Sex") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/allINMA.burden.sex.png")
sex.burden.plot
dev.off()


### Use other sex as reference ####
load("results/epimutations/INMA0combined.raw.epimutations.boys.girlsref.residuals.Rdata")
load("results/epimutations/INMA0combined.raw.epimutations.girls.boysref.residuals.Rdata")

load("results/epimutations/INMA4combined.raw.epimutations.boys.girlsref.residuals.Rdata")
load("results/epimutations/INMA4combined.raw.epimutations.girls.boysref.residuals.Rdata")


load("results/epimutations/HELIX.epimutations.raw.boys.girlsref.residuals.Rdata")
load("results/epimutations/HELIX.epimutations.raw.girls.boysref.residuals.Rdata")

inma0.girls <- make_res_df(res.inma0.girls.boysref.residuals.list, pheno = colData(inma0), 
                           dataset = "Newborn") 
inma0.girls$smoking <- ifelse(colData(inma0)[inma0.girls$sample, "msmk"] == "no smoking", "no", "yes")
inma0.girls.sum <- make_sum_df(inma0.girls)

inma0.boys <- make_res_df(res.inma0.boys.girlsref.residuals.list, pheno = colData(inma0), 
                           dataset = "Newborn") 
inma0.boys$smoking <- ifelse(colData(inma0)[inma0.boys$sample, "msmk"] == "no smoking", "no", "yes")
inma0.boys.sum <- make_sum_df(inma0.boys)

inma0.sex.sum <- rbind(inma0.girls.sum, inma0.boys.sum)

sex.diff0 <- inner_join(select(ungroup(inma0.sex.sum), method, sample, n), 
           select(ungroup(inma0.sum.df), sample, method, n, Sex, dataset), 
           by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)

sex.diff0 %>%
  group_by(method, Sex) %>%
  summarize(n_same = sum(diff == 0),
            n1 = sum(abs(diff) <= 1),
            n2 = sum(abs(diff) <= 2),
            n = n()) %>%
  mutate(p = n_same/n,
         p1 = n1/n,
         p2 = n2/n)


inma4.girls <- make_res_df(res.inma4.girls.boysref.residuals.list, pheno = colData(inma4), 
                           dataset = "4 years") 
inma4.girls$smoking <- ifelse(colData(inma4)[inma4.girls$sample, "msmk"] == "no smoking", "no", "yes")
inma4.girls.sum <- make_sum_df(inma4.girls)

inma4.boys <- make_res_df(res.inma4.boys.girlsref.residuals.list, pheno = colData(inma4), 
                          dataset = "4 years") 
inma4.boys$smoking <- ifelse(colData(inma4)[inma4.boys$sample, "msmk"] == "no smoking", "no", "yes")
inma4.boys.sum <- make_sum_df(inma4.boys)

inma4.sex.sum <- rbind(inma4.girls.sum, inma4.boys.sum)

sex.diff4 <- inner_join(select(ungroup(inma4.sex.sum), method, sample, n), 
                        select(ungroup(inma4.sum.df), sample, method, n, Sex, dataset), 
                        by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)



helix.girls <- make_res_df(res.helix.girls.boysref.residuals.list, pheno = colData(inma4), 
                           dataset = "8 years") 
helix.girls$smoking <- ifelse(helix_smk[helix.girls$sample, "msmok_a"] == 1, "yes", "no")
helix.girls.sum <- make_sum_df(helix.girls)

helix.boys <- make_res_df(res.helix.boys.girlsref.residuals.list, pheno = colData(inma4), 
                           dataset = "8 years") 
helix.boys$smoking <- ifelse(helix_smk[helix.boys$sample, "msmok_a"] == 1, "yes", "no")
helix.boys.sum <- make_sum_df(helix.boys)

helix.sex.sum <- rbind(helix.girls.sum, helix.boys.sum)

sex.diff8 <- inner_join(select(ungroup(helix.sex.sum), method, sample, n), 
                        select(ungroup(helix.sum.df), sample, method, n, Sex, dataset), 
                        by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)
# 
# sex.diff.plot <- rbind(sex.diff0, sex.diff4, sex.diff8) %>%
#   mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
#   mutate(Sex = ifelse(Sex == "F", "Girls", "Boys"),
#          dataset = factor(dataset, levels =  c("Newborn", "4 years", "8 years")),
#          diff = pmin(diff, 20),
#          diff = pmax(diff, -20)) %>%
#   ggplot(aes(x = diff, fill = Sex)) +
#   geom_histogram(binwidth = 1, position = "dodge") +
#   facet_grid(method ~ dataset) +
#   theme_bw() +
#   scale_x_continuous(name = "Differences in Epimutations detected")
# 
# png("figures/allINMA.sexDiffs.png")
# sex.diff.plot
# dev.off()


all.sex.res <- rbind(mutate(all.res.df, type = "All samples") %>% select(-mean_diff),
                     mutate(inma0.girls, type = "Girls"), 
                     mutate(inma0.boys, type = "Boys"), 
                     mutate(inma4.girls, type = "Girls"), 
                     mutate(inma4.boys, type = "Boys"), 
                     mutate(helix.girls, type = "Girls"), 
                     mutate(helix.boys, type = "Boys")) %>%
  mutate(dataset = factor(dataset, levels = c("Newborn", "4 years", "8 years"))) %>%
  filter(chromosome != 0)

  
### Recurrence ####
sex_cat <- read_delim("data/EWAS_catalog_sex.tsv", delim = "\t")
cord_sex <- subset(sex_cat, tissue == "Cord Blood")
blood_sex <- subset(sex_cat, grepl("Whole", tissue))
blood_child_sex <- subset(sex_cat, pmid == "33450751")


recur.sex.epi <- all.sex.res %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(dataset, type) %>%
  mutate(dataset_n = length(unique(sample))) %>%
  group_by(method, dataset, type, epi_region_id, dataset_n) %>%
  summarize(n = n(),
            freq = n()/dataset_n) %>%
  ungroup() %>%
  distinct() 

# 
# recur.epi.comb.freq.plot <- recur.sex.epi %>%
#   ggplot(aes(x = freq*100, fill = type)) +
#   geom_histogram(position = "dodge") +
#   theme_bw() +
#   facet_grid(method ~ dataset) +
#   scale_x_log10(name = "Epimutation frequency (%)", 
#                 breaks = c(0.1, 0.2, 0.5, 1, 2, 5))
# 
# png("figures/allINMA.epi.comb.freq.png")
# recur.epi.comb.freq.plot
# dev.off()
# 

recur_thres <- 0.005
recur.sex.top <- recur.sex.epi %>%
  complete(method, dataset, type, epi_region_id, fill = list(n = 1e-5, freq = 1e-5)) %>%
  group_by(method, dataset, epi_region_id) %>%
  select(-dataset_n, -n) %>% 
  pivot_wider(names_from = c(dataset, type), 
              values_from = "freq") %>%
  replace(is.na(.), 0) %>%
  filter(`Newborn_All samples` == 1e-5 & `4 years_All samples` == 1e-5, `8 years_All samples` == 1e-5) %>%
  filter(Newborn_Boys > recur_thres | Newborn_Girls > recur_thres | `4 years_Boys` > recur_thres |
           `4 years_Girls` > recur_thres | `8 years_Boys` > recur_thres | 
           `8 years_Girls` > recur_thres)


recur.sex.top %>%
  group_by(method) %>%
  summarize(n_girls = sum(Newborn_Girls > recur_thres | `4 years_Girls` > recur_thres | `8 years_Girls` > recur_thres ),
            n_boys = sum(Newborn_Boys > recur_thres | `4 years_Boys` > recur_thres | `8 years_Boys` > recur_thres ),
            n_girls_com = sum(Newborn_Girls > recur_thres & `4 years_Girls` > recur_thres & `8 years_Girls` > recur_thres ),
            n_boys_com = sum(Newborn_Boys > recur_thres & `4 years_Boys` > recur_thres & `8 years_Boys` > recur_thres ))

epi_rec_sex_list <- lapply(methods, function(x) {
  filter(recur.sex.top, method == x) %>%
    arrange(desc(Newborn_Boys))
})
writexl::write_xlsx(epi_rec_sex_list, path = "figures/INMA.Epimutations.recurrenceSex.xlsx")


recur.sex.top.coherent <- recur.sex.top %>%
  group_by(method) %>%
  filter((Newborn_Girls > recur_thres & `4 years_Girls` > recur_thres & `8 years_Girls` > recur_thres ) |
          (Newborn_Boys > recur_thres & `4 years_Boys` > recur_thres & `8 years_Boys` > recur_thres ))


recur.sex.top.coherent$cpgs <- lapply(recur.sex.top.coherent$epi_region_id, function(reg){
  
  gr <- candRegsGR[reg]
  
  cpgs <- names(subsetByOverlaps(rowRanges(inma0), gr))
})

recur.sex.top.coherent$n_cpgs <- lengths(recur.sex.top.coherent$cpgs)
recur.sex.top.coherent$n_cord <- sapply(recur.sex.top.coherent$cpgs, 
                                        function(x) sum(x %in% unique(cord_sex$cpg)))

recur.sex.top.coherent$n_blood <- sapply(recur.sex.top.coherent$cpgs, 
                                        function(x) sum(x %in% unique(blood_sex$cpg)))
recur.sex.top.coherent$n_child <- sapply(recur.sex.top.coherent$cpgs, 
                                         function(x) sum(x %in% unique(blood_child_sex$cpg)))

recur.sex.top.coherent %>% 
  summarize(N = sum(n_cpgs),
    n_cord = sum(n_cord),
    n_blood = sum(n_blood),
    n_child = sum(n_child)) %>%
  mutate(p_cord = n_cord/N,
         p_blood = n_blood/N,
         p_child = n_child/N)
    
plotSex <- function(set, range){
  
  miniset <- subsetByOverlaps(set, range)
  
  df <- getBeta(miniset)

  df <- t(df) %>% data.frame()
  df$Sex <- ifelse(set$Sex == "F", "Girls", "Boys")
  df$id <- colnames(miniset)
  
  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$pos <- start(rowRanges(miniset)[df.gath$cpg])
  
  df.sup <- df.gath %>%
    group_by(pos, Sex) %>%
    summarize(m = median(methylation))
    
  ggplot(df.gath, aes(x = pos, y = methylation, group = id, col = Sex)) +
    geom_point(alpha = 0.15) +
    geom_line(alpha = 0.15) +
    geom_line(data = df.sup, aes(x = pos, y = m, col = Sex, group = Sex), size = 2, linetype = "dashed") +
    scale_y_continuous(name = "DNA methylation", limits = c(0, 1)) +
    scale_x_continuous(name = "Coordinates") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

plotRecurrentReg <- function(reg){
  
  p0 <- plotSex(inma0, reg) + ggtitle("Newborns") + theme(legend.position = "none")
  p4 <- plotSex(inma4, reg) + ggtitle("4 years") + theme(legend.position = "none")
  p8 <- plotSex(helix, reg) + ggtitle("8 years")
  
  plot_grid(p0, p4, p8, ncol = 3, rel_widths = c(2, 2, 3) )
}
rec.sex.quant.plot <- plotRecurrentReg(candRegsGR["chr1_75590483"])
rec.sex.beta.plot <- plotRecurrentReg(candRegsGR["chr6_150346497"])


png("figures/INMA.sex.rec.epis.png", width = 1200)
plot_grid(rec.sex.quant.plot, rec.sex.beta.plot, ncol = 1, labels = c("A", "B"))
dev.off()

# ggplot(aes(x = type, group = epi_region_id, y = freq)) +
#   geom_line() +
#   geom_point() +
#   theme_bw() +
#   facet_grid(method ~ dataset) +
#   scale_y_log10(name = "Epimutation frequency (%)", 
#                 breaks = c(0.1, 0.2, 0.5, 1, 2, 5))



a <- res.helix.df %>% 
  group_by(method, epi_region_id) %>%
  summarize(n = n()) %>%
  mutate(p = n/length(unique(res.helix.df$sample)))


b <- helix.boys %>% 
  group_by(method, epi_region_id) %>%
  summarize(n = n()) %>%
  mutate(p = n/length(unique(helix.boys$sample)))

boys.regs8 <- lapply(methods, function(x) {
  sort(table( subset(helix.boys, method == x)$epi_region_id))
})


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


batch_tabs3 <- lapply(methods, function(m){
  tab <- inma0.sum.df %>%
    filter(method %in% m) %>%
    mutate(cat = ifelse(n >= 20, "No epi", "epi")) 
  table(tab$cat, tab$batch)
})
lapply(batch_tabs3, fisher.test)


inma0.batch.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ batch, data = inma0.sum.df, family = "poisson", subset = method == x & n > 0))
})

batch.burden.plot <- inma0.sum.df %>%
  filter(method %in% c("quantile", "beta", "mlm")) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm")),
         batch = ifelse(batch == "Esteller", "Reference", "Alternative"),
         batch = factor(batch, levels = c("Reference", "Alternative"))) %>%
  group_by(method, batch, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, n_cat, batch, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = batch, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(. ~ method) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Batch") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/INMA0.burden.batch.png", height = 250, width = 500)
batch.burden.plot
dev.off()


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
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm")),
         cohort = factor(ifelse(cohort == "SAB", "INMA", cohort)),
         cohort = relevel(cohort, "INMA")) %>%
  group_by(method, cohort, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, n_cat, cohort, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = cohort, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(. ~ method) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Cohort") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/HELIX.burden.cohort.png", width = 800, height = 300)
cohort.burden.plot
dev.off()

helix.cohort.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ cohort, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})


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

helix.smk.mod.adj <- lapply(methods, function(x) {
  summary(glmrob(n ~ smoking + cohort, data = helix.sum.df, family = "poisson", subset = method == x & n > 0))
})

epi_smk_tabs_8cohort <- lapply(methods, function(m){
  lapply(unique(helix.sum.df$cohort), function(coh){
    tab <- helix.sum.df %>%
      filter(method %in% m & cohort == coh) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi"),
             cat = factor(cat, levels = c("No epi", "epi")) )
    table(tab$cat, tab$smoking)
  })
})
lapply(unlist(epi_smk_tabs_8cohort, recursive = FALSE), fisher.test)


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
  facet_grid(dataset ~ method) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Maternal smoking during pregnancy") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample") 



png("figures/allINMA.burden.smoking.png")
smk.burden.plot
dev.off()

### Smoking vs non-smoking ####
load("results/epimutations/INMA0combined.raw.epimutations.smoking.residuals.Rdata")
load("results/epimutations/INMA4combined.raw.epimutations.smoking.residuals.Rdata")
load("results/epimutations/HELIX.epimutations.raw.smoking.residuals.Rdata")

inma0.smoking <- make_res_df(res.inma0.smoking.residuals.list, pheno = colData(inma0), 
                           dataset = "Newborn") 
inma0.smoking$smoking <- "yes"
inma0.smoking.sum <- make_sum_df(inma0.smoking)
smoking.diff0 <- inner_join(select(ungroup(inma0.smoking.sum), method, sample, n), 
                        select(ungroup(inma0.sum.df), sample, method, n, Sex, dataset), 
                        by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)

smoking.diff0 %>%
  group_by(method) %>%
  summarize(p = mean(diff == 0))


inma4.smoking <- make_res_df(res.inma4.smoking.residuals.list, pheno = colData(inma4), 
                             dataset = "Newborn") 
inma4.smoking$smoking <- "yes"
inma4.smoking.sum <- make_sum_df(inma4.smoking)
smoking.diff4 <- inner_join(select(ungroup(inma4.smoking.sum), method, sample, n), 
                            select(ungroup(inma4.sum.df), sample, method, n, Sex, dataset), 
                            by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)

helix.smoking <- make_res_df(res.helix.smoking.residuals.list, pheno = colData(helix), 
                             dataset = "Newborn") 
helix.smoking$smoking <- "yes"
helix.smoking.sum <- make_sum_df(helix.smoking)
smoking.diff8 <- inner_join(select(ungroup(helix.smoking.sum), method, sample, n), 
                            select(ungroup(helix.sum.df), sample, method, n, Sex, dataset), 
                            by = c("sample", "method")) %>%
  mutate(diff = n.x - n.y)

smoking.diff.plot <- rbind(smoking.diff0, smoking.diff4, smoking.diff8) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  mutate(dataset = factor(dataset, levels =  c("Newborn", "4 years", "8 years")),
         diff = pmin(diff, 20),
         diff = pmax(diff, -20)) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_grid(method ~ dataset) +
  theme_bw() +
  scale_x_continuous(name = "Differences in Epimutations detected")

png("figures/allINMA.smokingDiffs.png")
smoking.diff.plot
dev.off()

### Recurrence ####
all.smoking.res <- rbind(mutate(all.res.df, type = "All samples") %>% select(-mean_diff),
                     mutate(inma0.smoking, type = "Smokers"), 
                     mutate(inma4.smoking, type = "Smokers"), 
                     mutate(helix.smoking, type = "Smokers")) %>%
  mutate(dataset = factor(dataset, levels = c("Newborn", "4 years", "8 years"))) %>%
  filter(chromosome != 0)


recur.smoking.epi <- all.smoking.res %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(dataset, type) %>%
  mutate(dataset_n = length(unique(sample))) %>%
  group_by(method, dataset, type, epi_region_id, dataset_n) %>%
  summarize(n = n(),
            freq = n()/dataset_n) %>%
  ungroup() %>%
  distinct() 

recur_thres <- 0.005
recur.smoking.top <- recur.smoking.epi %>%
  complete(method, dataset, type, epi_region_id, fill = list(n = 1e-5, freq = 1e-5)) %>%
  group_by(method, dataset, epi_region_id) %>%
  select(-dataset_n, -n) %>% 
  pivot_wider(names_from = c(dataset, type), 
              values_from = "freq") %>%
  replace(is.na(.), 0) %>%
  filter(`Newborn_All samples` == 1e-5 & `4 years_All samples` == 1e-5, `8 years_All samples` == 1e-5) %>%
  filter(Newborn_Smokers > recur_thres | `4 years_Smokers` > recur_thres | `8 years_Smokers` > recur_thres)


writexl::write_xlsx(recur.smoking.top, path = "figures/INMA.Epimutations.recurrenceSmoking.xlsx")



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

sab.rep.res$cpg_new <-  strsplit(sab.rep.res$Newborn, ",")
sab.rep.res$cpg_4y <-  strsplit(sab.rep.res$`4 years`, ",")
sab.rep.res$cpg_8y <-  strsplit(sab.rep.res$`8 years`, ",")

interFun <- function(x, y){
  if (is.na(x)[1] | is.na(y)[1]){
    return(TRUE)
  }
  if (length(intersect(x, y)) > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

sab.rep.res$new_4 <-  unlist(Map(interFun, sab.rep.res$cpg_new, sab.rep.res$cpg_4y))
sab.rep.res$new_8 <-  unlist(Map(interFun, sab.rep.res$cpg_new, sab.rep.res$cpg_8y))
sab.rep.res$y4_y8 <-  unlist(Map(interFun, sab.rep.res$cpg_4y, sab.rep.res$cpg_8y))

sab.rep.res <- subset(sab.rep.res, new_4)

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

# ## Check effect across time points ####
# sab.rep.res2.meth <-  sab.rep.res2 %>%
#   filter(method %in% methods)
# sab.rep.res2.meth$mean_diff0 <- unlist(mclapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
#   getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = inma0)
# }, mc.cores = 20))
# sab.rep.res2.meth$mean_diff4 <- unlist(mclapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
#   getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = inma4)
# }, mc.cores = 20))
# sab.rep.res2.meth$mean_diff8 <- unlist(mclapply(seq_len(nrow(sab.rep.res2.meth)), function(i) {
#   getMeanDifference(sab.rep.res2.meth[i, ]$cpg_list[[1]], sab.rep.res2.meth[i, ]$idnum, set = helix_sab)
# }, mc.cores = 20))
# 
# sab.replic.diff.plot <- sab.rep.res2.meth %>% 
#   group_by(method, sigDatasets2) %>%
#   slice_sample(n = 20) %>% 
#   ungroup() %>%
#   gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
#   mutate(age = gsub("mean_diff", "", age),
#          method = factor(method, levels = c("quantile", "beta", "mlm")),
#          category = dplyr::recode(sigDatasets2, "--8 years" = "8 years",
#                                   "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
#                                   "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
#                                   "Cord blood-4 years-" = "Cord blood + 4 years", 
#                                   "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
#          category = factor(category, levels = c("Cord blood", "4 years", 
#                                                 "8 years", "Cord blood + 4 years",
#                                                 "Cord blood + 8 years", "4 years + 8 years",
#                                                 "Cord blood + 4 years + 8 years"))) %>%
#   ggplot(aes(x = age, y = mean_diff, group = epi_region_id, color = category)) +
#   geom_line() + geom_point() +
#   facet_wrap(~ method) +
#   scale_y_continuous("Epimutation mean methylation difference") +
#   scale_color_discrete(name = "Epimutation signal detection") +
#   scale_x_discrete("Age (in years)") +
#   theme_bw()
# 
# png("figures/INMAsab.replic.diff.png")
# sab.replic.diff.plot
# dev.off()
# 
# 
# sab.replic.diff.coher.plot <- sab.rep.res2.meth %>% 
#   filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
#   group_by(method, sigDatasets) %>%
#   slice_sample(n = 20) %>% 
#   gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
#   mutate(age = gsub("mean_diff", "", age),
#          method = factor(method, levels = c("quantile", "beta", "mlm")),
#          category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
#                                   "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
#                                   "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
#                                   "Cord blood-4 years-" = "Cord blood + 4 years", 
#                                   "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
#          category = factor(category, levels = c("Cord blood", "4 years", 
#                                                 "8 years", "Cord blood + 4 years",
#                                                 "Cord blood + 8 years", "4 years + 8 years",
#                                                 "Cord blood + 4 years + 8 years"))) %>%
#   ggplot(aes(x = age, y = mean_diff, group = paste(epi_region_id, idnum), color = category)) +
#   geom_line() + geom_point() +
#   facet_wrap(~ method) +
#   scale_y_continuous("Epimutation mean methylation difference") +
#   scale_color_discrete(name = "Epimutation detection") +
#   scale_x_discrete("Age (in years)") +
# #   theme_bw()
# # 
# # png("figures/INMAsab.replic.diff.coherent.png")
# # sab.replic.diff.coher.plot
# # dev.off()
# # 
# replic.direc.mod <- lapply(methods, function(x) {
#   df <- sab.rep.res2.meth %>%
#     filter(method == x) %>%
#     filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
#     mutate(id = paste(idnum, epi_region_id)) %>%
#     gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
#     mutate(age = gsub("mean_diff", "", age),
#            age = as.numeric(age))
#   lmer(abs(mean_diff) ~ age + (1|id), df)
# })
# lapply(replic.direc.mod, car::Anova)
# # 
# 
# 
# sab.replic.diff.coher.plot2 <- sab.rep.res2.meth %>% 
#   filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
#   group_by(method, sigDatasets) %>%
#   gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
#   mutate(age = gsub("mean_diff", "", age),
#          method = factor(method, levels = c("quantile", "beta", "mlm")),
#          category = dplyr::recode(sigDatasets, "--8 years" = "8 years",
#                                   "-4 years-" = "4 years", "-4 years-8 years" = "4 years + 8 years",
#                                   "Cord blood--"  = "Cord blood", "Cord blood--8 years" = "Cord blood + 8 years",
#                                   "Cord blood-4 years-" = "Cord blood + 4 years", 
#                                   "Cord blood-4 years-8 years" = "Cord blood + 4 years + 8 years"),
#          category = factor(category, levels = c("Cord blood", "4 years", 
#                                                 "8 years", "Cord blood + 4 years",
#                                                 "Cord blood + 8 years", "4 years + 8 years",
#                                                 "Cord blood + 4 years + 8 years"))) %>%
#   ggplot(aes(x = age, y = abs(mean_diff))) +
# geom_boxplot() +  facet_wrap(~ method) +
#   scale_y_continuous("Epimutation mean methylation difference") +
#   scale_color_discrete(name = "Epimutation detection") +
#   scale_x_discrete("Age (in years)") +
#   theme_bw()
# 
# png("figures/INMAsab.replic.diff.coherent2.png")
# sab.replic.diff.coher.plot2
# dev.off()
# 

# ### Check batch in INMA0 ####
# batch_map <- colData(inma0)[, c("idnum", "Batch")]
# rownames(batch_map) <- batch_map$idnum
# quant.rep <- sab.rep.res2.meth %>% 
#   filter(sigDatasets2 == "Cord blood-4 years-8 years") %>%
#   filter(method == "quantile") %>%
#   mutate(batch = batch_map[as.character(idnum), "Batch"])
# 
# png("figures/INMAsab.replic.diff.quantile.png")
# quant.rep %>%
#   gather(age, mean_diff, c("mean_diff0", "mean_diff4", "mean_diff8")) %>%
#   mutate(age = gsub("mean_diff", "", age)) %>%
#   ggplot(aes(x = age, y = mean_diff, group = paste(epi_region_id, idnum), color = batch)) +
#   geom_line() + geom_point() +
#   scale_y_continuous("Epimutation mean methylation difference") +
#   scale_color_discrete(name = "Epimutation detection") +
#   scale_x_discrete("Age (in years)") +
#   theme_bw()
# dev.off()



# Correlation with gene expression ####
getGenesZ <- function(epi_df, row, gexp, gexpRange, window){
  
  exp_name <- epi_df[row, ][["SampleName"]]
  if (!exp_name %in% colnames(gexp)){
    return(  list(geneName = NA,
                  Zscore = NA,
                  outliers = NA,
                  distance = NA, 
                  rank = NA))
  }
  epi_range <- GRanges(gsub("_", ":", epi_df[row, ][["epi_region_id"]]))
  
  transRanges <- subsetByOverlaps(gexpRange, epi_range + window)
  if (length(transRanges) == 0){
    return(  list(geneName = NA,
                  Zscore = NA,
                  outliers = NA,
                  distance = NA, 
                  rank = NA))
  }
  
  selTranscripts <- names(transRanges)
  
  
  genesmat <- gexp[selTranscripts, , drop = FALSE]
  zmat <- apply(genesmat, 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
  zs <- zmat[exp_name, ]
  
  rankMat <- apply(genesmat, 1, rank)
  ranks <- rankMat[exp_name, ]
  
  
  isOutliersRow <- function(x){
    qs <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
    iq <- qs[2] - qs[1]
    out <- x < qs[1] - 1.5*iq | x > qs[2] + 1.5*iq
    out
  }
  outmat <- apply(genesmat, 1, isOutliersRow)
  outs <- outmat[exp_name, ]
  
  list(geneNames = selTranscripts,
       Zscores = zs,
       outliers = outs, 
       distances = distance(epi_range, transRanges),
       rank = ranks)
  
}



selectNearest <- function(x) {
  if (length(x$Zscore) == 1){
    data.frame(geneName = x$geneName, Zscore = x$Zscore, outliers = x$outliers,
               distance = x$Zscore, rank = x$rank) 
  }
  else {
    ind <- which.min(x$distance)
    data.frame(geneName =  x$geneName[ind], outliers = x$outliers[ind],
               Zscore = x$Zscore[ind], 
               distance = x$distance[ind], rank = x$rank[ind])
  }
}



isOutliers <- function(x, val){
  qs <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iq <- qs[2] - qs[1]
  out <- val < qs[1] - 1.5*iq | val > qs[2] + 1.5*iq
  out
}

getGenesTSS <- function(epi_df, row, gexp){
  
  gene <- unlist(epi_df[row, "tss_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  selTranscripts <- rownames(transcriptome_subcohort_f1)[fData(transcriptome_subcohort_f1)$GeneSymbolDB2 %in% gene]  
  if (length(selTranscripts) == 0){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  exp_name <- epi_df[row, ][["SampleName"]]
  
  if (!exp_name %in% colnames(gexp)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  gexp <- exprs(gexp)
  
  genesmat <- gexp[selTranscripts, , drop = FALSE]
  
  out <- getVals(genesmat, exp_name, selTranscripts)
}
getVals <- function(genesmat, exp_name, selTranscripts){
  zmat <- t(apply(genesmat, 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
  zs <- zmat[, exp_name]
  z <- zs[which.max(abs(zs))]
  
  genevec <- genesmat[which.max(abs(zs)), ]
  out <- isOutliers(genevec, genevec[exp_name])
  
  rank <- rank(genevec)[exp_name]
  
  list(z = z, out = out, rank = rank, gene = selTranscripts[which.max(abs(zs))])
  
}

getGeneseQTM <- function(epi_df, row, gexp){
  
  gene <- unlist(epi_df[row, "eqtm_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }

  exp_name <- epi_df[row, ][["SampleName"]]
  
  if (!exp_name %in% colnames(gexp)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  genesmat <- gexp[gene, , drop = FALSE]
  
  out <- getVals(genesmat, exp_name, gene)
}

## INMA4 ####
eqtm <- read_delim("data/eqtm.txt.gz", delim = "\t")

# inma4.exp.z <- t(apply(exprs(inma.expset), 1, function(x) (x - mean(x))/sd(x)))
inma4.exp.range <- makeGRangesFromDataFrame(fData(inma.expset))

res.inma4.filt <- res.inma4.df %>%
  filter(chromosome != 0 & method %in% methods)
res.inma4.filt$SampleName <- paste0("X04_", res.inma4.filt$idnum)

## Add epimutation magnitude
mean_diffs4 <- mclapply(seq_len(nrow(res.inma4.filt)), function(i) {
  getMeanDifference2(res.inma4.filt[i, ]$cpg_ids, res.inma4.filt[i, ]$sample, set = inma4)
}, mc.cores = 20)
res.inma4.filt$mean_diff <- unlist(mean_diffs4)

inma4.annot <- getAnnotation(inma4)

tss_genes4 <- mclapply(res.inma4.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  annot <- inma4.annot[cpgs, ]
  tss <- annot[grepl("TSS", annot$UCSC_RefGene_Group), , drop = FALSE]
  if (nrow(tss) == 0){
    return(NA)
  } else {
    gene <- strsplit(tss$UCSC_RefGene_Name, ";")
    unique(unlist(gene))
  }
}, mc.cores = 20)

res.inma4.filt$tss_genes <- tss_genes4
tss_list <- lapply(seq_len(nrow(res.inma4.filt)), getGenesTSS, 
                  epi_df = res.inma4.filt, gexp = inma.expset)
res.inma4.filt$tss_z <- sapply(tss_list, function(x) x$z)
res.inma4.filt$tss_out <- sapply(tss_list, function(x) x$out)
res.inma4.filt$tss_rank <- sapply(tss_list, function(x) x$rank)
res.inma4.filt$tss_gene <- sapply(tss_list, function(x) x$gene)


eqtm_genes4 <- mclapply(res.inma4.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  tab <- subset(eqtm, CpG %in% cpgs)
  if (nrow(tab) == 0){
    return(NA)
  } else {
    unique(tab$TC)  
    }
}, mc.cores = 20)

res.inma4.filt$eqtm_genes <- eqtm_genes4
eqtm_list <- lapply(seq_len(nrow(res.inma4.filt)), getGeneseQTM, 
                   epi_df = res.inma4.filt, gexp = exprs(inma.expset))
res.inma4.filt$eqtm_z <- sapply(eqtm_list, function(x) x$z)
res.inma4.filt$eqtm_out <- sapply(eqtm_list, function(x) x$out)
res.inma4.filt$eqtm_rank <- sapply(eqtm_list, function(x) x$rank)
res.inma4.filt$eqtm_gene <- sapply(eqtm_list, function(x) x$gene)

geneInfo4 <- mclapply(seq_len(nrow(res.inma4.filt)), getGenesZ, epi_df = res.inma4.filt,
                    gexp = exprs(inma.expset), gexpRange = inma4.exp.range, window = 250e3, mc.cores = 10)
# geneInfo4df <- Reduce(rbind, lapply(geneInfo4, selectTop))
# res.inma4.filt <- cbind(res.inma4.filt, geneInfo4df)

geneInfo4Near <- Reduce(rbind, lapply(geneInfo4, selectNearest))
colnames(geneInfo4Near) <- paste0(colnames(geneInfo4Near), "Near" )
res.inma4.filt$near_z <- geneInfo4Near$ZscoreNear
res.inma4.filt$near_out <- geneInfo4Near$outliersNear
res.inma4.filt$near_rank <- geneInfo4Near$rankNear
res.inma4.filt$near_gene <- geneInfo4Near$geneNameNear


res.inma4.filt.sumdf  <- res.inma4.filt %>% 
  select(method, mean_diff, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 3:8) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


inma4.gexp.plot <- res.inma4.filt.sumdf %>% 
  ggplot(aes(x = exp_type, y = value, color = method)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5, 
               binwidth = 0.2) +
  theme_bw() +
  facet_grid(measure ~ method, scales = "free") +
  scale_x_discrete(name = "Gene mapping")

png("figures/INMA4.genexp.png")
inma4.gexp.plot
dev.off()
# 
# inma4.gexp.methdiff.plot <- res.inma4.filt.sumdf %>%
#   filter(measure == "z") %>%
#   ggplot(aes(x = abs(mean_diff), y = abs(value), color = method)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   facet_grid(exp_type ~ method, scales = "free")
# 
# 
# png("figures/INMA4.genexp.methdiff.png")
# inma4.gexp.methdiff.plot
# dev.off()

meth.gexp.assoc4 <- lapply(methods, function(m){
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.inma4.filt.sumdf %>%
      filter(measure == "z" & method ==m & exp_type == me) 
    summary(lm(abs(value) ~ abs(mean_diff), tab))  
  })
})


# inma4.gexp.eqtm.plot <- res.inma4.filt %>% 
#   group_by(method, idnum) %>%
#   ggplot(aes(x = abs(mean_diff), y = abs(eqtm_z))) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   facet_grid(~ method)
# 
# png("figures/INMA4.genexp.eqtm.png")
# inma4.gexp.eqtm.plot
# dev.off()
# 
# res.inma4.filt$pvals <- pnorm(abs(res.inma4.filt$Zscore), lower = FALSE)
# inma4.gexp.pvals.plot <- res.inma4.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(Zscore)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
#   ggplot(aes(x = method, y = -log10(pvals))) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ n_cat)
# 
# png("figures/INMA4.genexp.pvals.png")
# inma4.gexp.pvals.plot
# dev.off()
# 
# 
# cot <- res.inma4.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(Zscore)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))
# 
# summary(lm(abs(Zscore) ~ n_cat, cot, subset = method == "beta"))  
# summary(lm(abs(Zscore) ~ n_cat, cot, subset = method == "mlm"))  
# 
# 
# cot <- res.inma4.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(ZscoreNear)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))
# 
# summary(lm(abs(ZscoreNear) ~ n_cat, cot, subset = method == "beta"))  
# summary(lm(abs(ZscoreNear) ~ n_cat, cot, subset = method == "mlm"))  
# 

## HELIX ####
helix.gexp <- transcriptome_subcohort_f1[, transcriptome_subcohort_f1$cohort != "KAN" & 
                                           transcriptome_subcohort_f1$h_ethnicity_3cat == "WhiteEur_WhiteOther"] 

# helix.exp.z <- t(apply(exprs(helix.gexp), 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
helix.exp.range <- makeGRangesFromDataFrame(fData(helix.gexp))

res.helix.filt <- res.helix.df %>%
  filter(chromosome != 0 & method %in% methods)

res.helix.filt$SampleName <- res.helix.filt$sample

## Add epimutation magnitude
mean_diffs8 <- mclapply(seq_len(nrow(res.helix.filt)), function(i) {
  getMeanDifference2(res.helix.filt[i, ]$cpg_ids, res.helix.filt[i, ]$sample, set = helix)
}, mc.cores = 20)
res.helix.filt$mean_diff <- unlist(mean_diffs8)

helix.annot <- getAnnotation(helix)

tss_genes8 <- mclapply(res.helix.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  annot <- helix.annot[cpgs, ]
  tss <- annot[grepl("TSS", annot$UCSC_RefGene_Group), , drop = FALSE]
  if (nrow(tss) == 0){
    return(NA)
  } else {
    gene <- strsplit(tss$UCSC_RefGene_Name, ";")
    unique(unlist(gene))
  }
}, mc.cores = 20)

res.helix.filt$tss_genes <- tss_genes8

tss_list8 <- lapply(seq_len(nrow(res.helix.filt)), getGenesTSS, 
                   epi_df = res.helix.filt, gexp = helix.gexp)
res.helix.filt$tss_z <- sapply(tss_list8, function(x) x$z)
res.helix.filt$tss_out <- sapply(tss_list8, function(x) x$out)
res.helix.filt$tss_rank <- sapply(tss_list8, function(x) x$rank)
res.helix.filt$tss_gene <- sapply(tss_list8, function(x) x$gene)



eqtm_genes8 <- mclapply(res.helix.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  tab <- subset(eqtm, CpG %in% cpgs)
  if (nrow(tab) == 0){
    return(NA)
  } else {
    unique(tab$TC)  
  }
}, mc.cores = 20)

res.helix.filt$eqtm_genes <- eqtm_genes8

eqtm_list8 <- lapply(seq_len(nrow(res.helix.filt)), getGeneseQTM, 
                    epi_df = res.helix.filt, gexp = exprs(helix.gexp))
res.helix.filt$eqtm_z <- sapply(eqtm_list8, function(x) x$z)
res.helix.filt$eqtm_out <- sapply(eqtm_list8, function(x) x$out)
res.helix.filt$eqtm_rank <- sapply(eqtm_list8, function(x) x$rank)
res.helix.filt$eqtm_gene <- sapply(eqtm_list8, function(x) x$gene)

geneInfo8 <- mclapply(seq_len(nrow(res.helix.filt)), getGenesZ, epi_df = res.helix.filt,
                    gexp = exprs(helix.gexp), gexpRange = helix.exp.range, window = 250e3, mc.cores = 10)

geneInfo8Near <- Reduce(rbind, lapply(geneInfo8, selectNearest))
res.helix.filt$near_z <- geneInfo8Near$Zscore
res.helix.filt$near_out <- geneInfo8Near$out
res.helix.filt$near_rank <- geneInfo8Near$rank
res.helix.filt$near_gene <- geneInfo8Near$geneName


res.helix.filt.sumdf  <- res.helix.filt %>% 
  select(method, mean_diff, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 3:8) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


helix.gexp.plot <- res.helix.filt.sumdf %>% 
  ggplot(aes(x = exp_type, y = value, color = method)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5,
               binwidth = 0.2) +
  theme_bw() +
  facet_grid(measure ~ method, scales = "free")  +
  scale_x_discrete(name = "Gene mapping")

png("figures/HELIX.genexp.png")
helix.gexp.plot
dev.off()
# 
# helix.gexp.methdiff.plot <- res.helix.filt.sumdf %>%
#   filter(measure == "z") %>%
#   ggplot(aes(x = abs(mean_diff), y = abs(value), color = method)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   facet_grid(exp_type ~ method, scales = "free")
# 
# 
# png("figures/HELIX.genexp.methdiff.png")
# helix.gexp.methdiff.plot
# dev.off()
# 
meth.gexp.assoc8 <- lapply(methods, function(m){
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.helix.filt.sumdf %>%
      filter(measure == "z" & method ==m & exp_type == me) 
    summary(lm(abs(value) ~ abs(mean_diff), tab))  
  })
})

res.comb.filt.sumdf <- rbind(res.inma4.filt %>% 
                               select(method, ends_with("out")) %>%
                               gather(Measure, value, 2:4) %>%
                               mutate(age = "4 years"),
                             res.helix.filt %>% 
                               select(method, ends_with("out")) %>%
                               gather(Measure, value, 2:4) %>%
                               mutate(age = "8 years")) %>%
  filter(!is.na(value)) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


all.prop.gexp.plot <- res.comb.filt.sumdf %>%
  group_by(method, exp_type, age) %>%
  summarize(p = mean(value)) %>%
  ggplot(aes(x = exp_type, y = p*100, fill = method)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(age ~ method) +
  scale_y_continuous("Epimutaitons with outlier expression (%)", limits = c(0, 50)) +
  scale_x_discrete(name = "Gene mapping")

png("figures/allINMA.genexp.prop.png", height = 350)
all.prop.gexp.plot
dev.off()

all.prop.gexp.tab <- lapply(c("4 years", "8 years"), function(a) {
  
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.comb.filt.sumdf %>%
      filter(age  == a & exp_type == me) 
    table(tab$method, tab$value)
  })
})
lapply(all.prop.gexp.tab, lapply, chisq.test)

all.prop.gexp.mod <- lapply(c("4 years", "8 years"), function(a) {
  
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.comb.filt.sumdf %>%
      filter(age  == a & exp_type == me) 
    summary(glm(value ~ method, tab, family = "binomial"))
  })
})


all.prop.gexp.mod2 <- lapply(c("4 years", "8 years"), function(a) {
  
  lapply(methods, function(me){
    tab <- res.comb.filt.sumdf %>%
      filter(age  == a & method == me) 
    summary(glm(value ~ exp_type, tab, family = "binomial"))
  })
})


# 
# geneInfoHelix <- lapply(seq_len(nrow(res.helix.filt)), getGenesZ, epi_df = data.frame(res.helix.filt),
#                     gexpZ = helix.exp.z, gexpRange = helix.exp.range, window = 250e3)
# 
# 
# geneInfohelixdf <- Reduce(rbind, lapply(geneInfoHelix, selectTop))
# res.helix.filt <- cbind(res.helix.filt, geneInfohelixdf)
# 
# geneInfohelixNear <- Reduce(rbind, lapply(geneInfoHelix, selectNearest))
# colnames(geneInfohelixNear) <- paste0(colnames(geneInfohelixNear), "Near" )
# res.helix.filt <- cbind(res.helix.filt, geneInfohelixNear)
# 
# helix.gexp.plot <- res.helix.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(Zscore)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
#   ggplot(aes(color = method, x = Zscore)) +
#   geom_density() +
#   theme_bw() +
#   facet_grid(n_cat ~ .)
# 
# png("figures/HELIX.genexp.png")
# helix.gexp.plot
# dev.off()
# 
# res.helix.filt$pvals <- pnorm(abs(res.helix.filt$Zscore), lower = FALSE)
# helix.gexp.pvals.plot <- res.helix.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(Zscore)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers")) %>%
#   ggplot(aes(x = method, y = -log10(pvals))) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(~ n_cat)
# 
# png("figures/HELIX.genexp.pvals.png")
# helix.gexp.pvals.plot
# dev.off()
# 
# 
# pac <- res.helix.filt %>% 
#   group_by(method, idnum) %>%
#   filter(!is.na(Zscore)) %>%
#   mutate(n = n(),
#          n_cat = ifelse(n > 5, "Samples with >5 outlier", "Samples with 1-5 outliers"))
# 
# summary(lm(abs(Zscore) ~ n_cat, pac, subset = method == "beta"))  
# summary(lm(abs(Zscore) ~ n_cat, pac, subset = method == "mlm"))  


## Export epimtuations with gene expresion ####
gexp_res_list <- list(`4 years` = subset(res.inma4.filt[, c(20, 16, 13, 3:5, 7, 25, 11, 28:29, 33, 26:27, 34, 30:32)], 
                                         !is.na(tss_z) | !is.na(eqtm_z) | !is.na(near_z)),
                      `8 years` = subset(res.helix.filt[, c(20, 16, 13, 3:5, 7, 23, 11, 28:29, 34, 25:26, 33, 30:32)], 
                                         !is.na(tss_z) | !is.na(eqtm_z) | !is.na(near_z)))

writexl::write_xlsx(gexp_res_list, path = "tables/INMA.Epimutations.geneExp.xlsx")

## Overlap with obesity genes

library(disgenet)
pass <- pass
disgenet_api_key <- get_disgenet_api_key(
  email = "carlos.ruiza@upf.edu", 
  password = "" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)




obe_genes <- read_delim("data/obesity_genes.txt", delim = "\t", col_names = FALSE)
helix.epi.gexp <- subset(res.helix.filt[, c(2, 20, 16, 13, 3:5, 7, 23, 11, 28:29, 34, 25:26, 33, 30:32)], 
       !is.na(tss_z) | !is.na(eqtm_z) | !is.na(near_z))
helix.epi.gexp <- subset(helix.epi.gexp, (abs(eqtm_z) > 2 & (eqtm_rank < 10 | eqtm_rank > 1010)) | 
                           (abs(tss_z) > 2 & (tss_rank < 10 | tss_rank > 1010)) | 
                           (abs(near_z) > 2 & (near_rank < 10 | near_rank > 1010)))
helix.epi.gexp$eqtm_symbol <- fData(helix.gexp)[helix.epi.gexp$eqtm_gene, "GeneSymbolDB2"]  
helix.epi.gexp$tss_symbol <- fData(helix.gexp)[helix.epi.gexp$tss_gene, "GeneSymbolDB2"]  
helix.epi.gexp$near_symbol <- fData(helix.gexp)[helix.epi.gexp$near_gene, "GeneSymbolDB2"]  

helix.epi.filt <- subset(helix.epi.gexp, (!is.na(eqtm_symbol) & eqtm_symbol != "") | 
                           (!is.na(tss_symbol) & tss_symbol != "") | 
                           (!is.na(near_symbol) & near_symbol != ""))
selgenes <- c(subset(helix.epi.filt, abs(eqtm_z) > 2)$eqtm_symbol, 
              subset(helix.epi.filt, abs(tss_z) > 2)$tss_symbol, 
              subset(helix.epi.filt, abs(near_z) > 2)$near_symbol)
selgenes <- unique(selgenes)

helix.epi.obe <- subset(helix.epi.gexp, sample %in% obe)
selgenes2 <- c(helix.epi.obeX1)$eqtm_symbol, 
              subset(helix.epi.obe, tss_symbol %in% obe_genes$X1)$tss_symbol, 
              subset(helix.epi.obe, near_symbol %in% obe_genes$X1)$near_symbol)

res_df <- Reduce(rbind, lapply(selgenes, function(x) {
  res <- gene2disease(gene = x)
  if (is.character(res)){
    return(NULL)
  }
  extract(res)
}))
diseases <- c("Autism Spectrum Disorders", "Obesity", "Intellectual Disability", 
              "MENTAL RETARDATION, AUTOSOMAL DOMINANT 39", "Neurodevelopmental Disorders")
subset(gq, disease_name %in% diseases)

subset(helix.epi.filt, eqtm_symbol == "MYT1L" | tss_symbol == "MYT1L" | near_symbol == "MYT1L" )

## Replicates ####
rep_epis <- subset(sab.rep.res2, sigDatasets == "Cord blood-4 years-8 years")

### INMA 4 ####
res.inma4.rep <- res.inma4.filt %>%
  mutate(epi_id = paste(method, idnum, epi_region_id)) %>%
  filter(epi_id %in% mutate(rep_epis, epi_id = paste(method, idnum, epi_region_id))$epi_id)

res.inma4.rep.sumdf  <- res.inma4.rep %>% 
  select(method, mean_diff, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 3:8) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


inma4.gexp.rep.plot <- res.inma4.rep.sumdf %>% 
  ggplot(aes(x = exp_type, y = value, color = method)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5) +
  theme_bw() +
  facet_grid(measure ~ method, scales = "free")

png("figures/INMA4.genexp.rep.png")
inma4.gexp.rep.plot
dev.off()

inma4.gexp.methdiff.rep.plot <- res.inma4.rep.sumdf %>%
  filter(measure == "z") %>%
  ggplot(aes(x = abs(mean_diff), y = abs(value), color = method)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(exp_type ~ method, scales = "free")


png("figures/INMA4.genexp.rep.methdiff.png")
inma4.gexp.methdiff.rep.plot
dev.off()

meth.gexp.rep.assoc4 <- lapply(methods, function(m){
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.inma4.rep.sumdf %>%
      filter(measure == "z" & method ==m & exp_type == me) 
    summary(lm(abs(value) ~ abs(mean_diff), tab))  
  })
})

### HELIX ####
res.helix.rep <- res.helix.filt %>%
  mutate(epi_id = paste(method, idnum, epi_region_id)) %>%
  filter(epi_id %in% mutate(rep_epis, epi_id = paste(method, idnum, epi_region_id))$epi_id)

res.helix.rep.sumdf  <- res.helix.rep %>% 
  select(method, mean_diff, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 3:8) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


helix.gexp.rep.plot <- res.helix.rep.sumdf %>% 
  ggplot(aes(x = exp_type, y = value, color = method)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5,
               binwidth = 0.2) +
  theme_bw() +
  facet_grid(measure ~ method, scales = "free")

png("figures/HELIX.genexp.rep.png")
helix.gexp.rep.plot
dev.off()

helix.gexp.methdiff.rep.plot <- res.helix.rep.sumdf %>%
  filter(measure == "z") %>%
  ggplot(aes(x = abs(mean_diff), y = abs(value), color = method)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(exp_type ~ method, scales = "free")


png("figures/HELIX.genexp.rep.methdiff.png")
helix.gexp.methdiff.rep.plot
dev.off()

meth.gexp.rep.assoc8 <- lapply(methods, function(m){
  lapply(c("eqtm", "tss", "near"), function(me){
    tab <- res.helix.rep.sumdf %>%
      filter(measure == "z" & method ==m & exp_type == me) 
    summary(lm(abs(value) ~ abs(mean_diff), tab))  
  })
})

res.comb.rep.sumdf <- rbind(mutate(res.inma4.rep.sumdf, age = "4 years", N = ncol(inma4)),
                             mutate(res.helix.rep.sumdf, age = "8 years", N = ncol(helix)))

allrep.prop.gexp.plot <- res.comb.rep.sumdf %>%
  filter(measure == "rank") %>%
  group_by(method, exp_type, age) %>%
  summarize(p = mean(value < N*0.05 | value > N*0.95, na.rm = TRUE)) %>%
  ggplot(aes(x = exp_type, y = p*100, fill = method, color = method)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(age ~ method) +
  scale_y_continuous("Percentage of samples in top 5% (%)")

png("figures/allINMA.rep.genexp.prop.png")
allrep.prop.gexp.plot
dev.off()


## Create table with persistant epimutations and gene expression  ####
rep_epis$chromosome <- sapply(strsplit(rep_epis$epi_region_id, "_"),`[`, 1)

combRanges <- unique(c(rowRanges(inma0), rowRanges(inma4), rowRanges(helix)))

rep_epis$start <- sapply(rep_epis$cpg_list, function(x) {
  min(start(combRanges[x]))
})
rep_epis$end <- sapply(rep_epis$cpg_list, function(x) {
  max(end(combRanges[x]))
})
mean_diffs0 <- mclapply(seq_len(nrow(rep_epis)), function(i) {
  sample <- colnames(inma0)[inma0$SampleID == rep_epis[i, ]$idnum]
  getMeanDifference2(rep_epis[i, ]$Newborn, sample, set = inma0)
}, mc.cores = 20)
rep_epis$mean_diff0 <- unlist(mean_diffs0)


pers_epi_df <- rep_epis %>% 
  mutate(epi_id = paste(method, idnum, epi_region_id)) %>%
  select(epi_id, method, idnum, epi_region_id, chromosome, start, end, mean_diff0) %>%
  left_join(res.inma4.filt %>% 
              mutate(epi_id = paste(method, idnum, epi_region_id)) %>%
              select(epi_id, starts_with(c("eqtm", "tss", "near")), mean_diff), 
            by = "epi_id") %>%
  left_join(res.helix.rep %>% 
              mutate(epi_id = paste(method, idnum, epi_region_id)) %>%
              select(epi_id, starts_with(c("eqtm", "tss", "near")), mean_diff), 
            by = "epi_id")

colnames(pers_epi_df) <- gsub("y", "8years", colnames(pers_epi_df))
colnames(pers_epi_df) <- gsub("x", "4years", colnames(pers_epi_df))
writexl::write_xlsx(pers_epi_df[, colnames(pers_epi_df)[c(3, 2, 4:8, 20, 32, 10:12, 14:19, 22:24, 26:31)]],
                    path = "tables/INMA.Epimutations.persistentEpis.xlsx")



# sab.rep.gexp <- sab.rep.res2 %>%
#   filter(method %in% methods)
# 
# sab.rep.gexp$SampleName <-  paste0("X04_", sab.rep.gexp$idnum)
# geneInfoSab4 <- lapply(seq_len(nrow(sab.rep.gexp)), getGenesZ, epi_df = data.frame(sab.rep.gexp),
#                    gexpZ = inma4.exp.z, gexpRange = inma4.exp.range, window = 250e3)
# 
# sab.rep.gexp$SampleName <-  colnames(helix_sab)[match(sab.rep.gexp$idnum, helix_sab$idnum)]
# geneInfoSab8 <- lapply(seq_len(nrow(sab.rep.gexp)), getGenesZ, epi_df = data.frame(sab.rep.gexp),
#                        gexpZ = helix.exp.z, gexpRange = helix.exp.range, window = 250e3)
# 
# geneInfoSab4top <- Reduce(rbind, lapply(geneInfoSab4, selectTop))
# colnames(geneInfoSab4top) <- paste0(colnames(geneInfoSab4top), "Top4" )
# sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab4top)
# 
# geneInfoSab4Near <- Reduce(rbind, lapply(geneInfoSab4, selectNearest))
# colnames(geneInfoSab4Near) <- paste0(colnames(geneInfoSab4Near), "Near4" )
# sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab4Near)
# 
# geneInfoSab8top <- Reduce(rbind, lapply(geneInfoSab8, selectTop))
# colnames(geneInfoSab8top) <- paste0(colnames(geneInfoSab8top), "Top8" )
# sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab8top)
# 
# geneInfoSab8Near <- Reduce(rbind, lapply(geneInfoSab8, selectNearest))
# colnames(geneInfoSab8Near) <- paste0(colnames(geneInfoSab8Near), "Near8" )
# sab.rep.gexp <- cbind(sab.rep.gexp, geneInfoSab8Near)
# 
# sab.gexp4.top.plot <- sab.rep.gexp %>% 
#   ggplot(aes(fill = sigDatasets2, x = ZscoreTop4)) +
#   geom_histogram(binwidth = 0.2) +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y") +
#   theme(legend.position = "none")
# 
# sab.gexp8.top.plot <- sab.rep.gexp %>% 
#   ggplot(aes(fill = sigDatasets2, x = ZscoreTop8)) +
#   geom_histogram(binwidth = 0.2) +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y")
# 
# png("figures/INMAsab.top.genexp.png", width = 1200)
# plot_grid(sab.gexp4.top.plot, sab.gexp8.top.plot, nrow = 1, rel_widths = c(3, 5))
# dev.off()
# 
# 
# 
# sab.gexp4.near.plot <- sab.rep.gexp %>% 
#   ggplot(aes(fill = sigDatasets2, x = ZscoreNear4)) +
#   geom_histogram(binwidth = 0.2) +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y") +
#   theme(legend.position = "none")
# 
# sab.gexp8.near.plot <- sab.rep.gexp %>% 
#   ggplot(aes(fill = sigDatasets2, x = ZscoreNear8)) +
#   geom_histogram(binwidth = 0.2) +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y")
# 
# png("figures/INMAsab.near.genexp.png", width = 1200)
# plot_grid(sab.gexp4.near.plot, sab.gexp8.near.plot, nrow = 1, rel_widths = c(3, 5))
# dev.off()
# 
# 
# sab.rep.gexp %>%
#   filter(!grepl("4 years", sigDatasets2)) %>%
#   pull(ZscoreTop4) %>%
#   abs() %>%
#   quantile(., 0.95, na.rm = TRUE)
# 
# 
# sab.rep.gexp %>%
#   filter(!grepl("8 years", sigDatasets2)) %>%
#   pull(ZscoreTop8) %>%
#   abs() %>%
#   quantile(., 0.95, na.rm = TRUE)
# 
# 
# sab.rep.gexp %>%
#   filter(!grepl("4 years", sigDatasets2)) %>%
#   pull(ZscoreNear4) %>%
#   abs() %>%
#   quantile(., 0.95, na.rm = TRUE)
# 
# 
# sab.rep.gexp %>%
#   filter(!grepl("8 years", sigDatasets2)) %>%
#   pull(ZscoreNear8) %>%
#   abs() %>%
#   quantile(., 0.95, na.rm = TRUE)
# 
# sab.gexp4.distance.signal.plot <- sab.rep.gexp %>% 
#   ggplot(aes(x = unlist(distance4), color = sigDatasets2, y = Zscore4)) +
#   geom_point() +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y")
# 
# png("figures/INMAsab.distance.signal.genexp.png")
# sab.gexp4.distance.signal.plot
# dev.off()
# 
# 
# 
# 
# png("figures/INMAsab.distanceDistr.signal.genexp.png")
# 
# sab.rep.gexp %>% 
#   ggplot(aes(y = unlist(distance4), x = sigDatasets2)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_grid(method ~ ., scales = "free_y")
# dev.off()
# 


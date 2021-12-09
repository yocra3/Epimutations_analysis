#'#################################################################################
#'#################################################################################
#' Analyze GSE51032 dataset
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)
library(tidyverse)
library(robustbase)

load("results/preprocess/GSE51032/GSE51032.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Select individuals with two main cancers  ####
gset$disease <- ifelse(is.na(gset$`cancer type (icd-10):ch1`), "Control",
                       ifelse(gset$`cancer type (icd-10):ch1` == "C18", "Colon", 
                              ifelse(gset$`cancer type (icd-10):ch1` == "C50", "Breast", NA)))
gset_filt <- gset[, !is.na(gset$disease)]

## Compute residuals of pcs  ####
beta <- meffil:::impute.matrix(getBeta(gset_filt), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(gset_filt), probe.range = 40000)
m <- getM(gset_filt)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(gset_filt) <- beta
save(gset_filt, file = "results/preprocess/GSE51032/GSE51032.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

case <- gset_filt[, gset_filt$disease != "Control"]
control <- gset_filt[, gset_filt$disease == "Control"]

## Run epimutations ####
methods <- c("beta", "quantile", "mlm")
names(methods) <- methods


res.gse51032.list <- lapply(methods, epimutations, case_samples = case, 
                            control_panel = control)
save(res.gse51032.list, file = "results/epimutations/GSE51032.epimutations.cases.residuals.Rdata")


## Process epimutations ####
nMethod <- sapply(res.gse51032.list, nrow)
nMethod[sapply(nMethod, is.null)] <- 0
res.gse51032.df <-  Reduce(rbind, res.gse51032.list) %>%
    mutate(method = rep(methods, unlist(nMethod))) %>%
    left_join(colData(gset_filt) %>% 
                data.frame() %>% 
                select(geo_accession, disease) %>% 
                mutate(sample = geo_accession))

res.gse51032.sum <- res.gse51032.df %>%
    group_by(method, sample, disease) %>%
    summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1", 
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


res.gse51032.sum %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ .) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Sex") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

epi_disease_tabs <- lapply(methods, function(m){
    tab <- res.gse51032.sum %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi")) 
    table(tab$cat, tab$disease)
  })
lapply(epi_disease_tabs, fisher.test)

disease.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ disease, data = res.gse51032.sum, family = "poisson", subset = method == x & n > 0))
})


## Get recurrent epimutations
recur.tumor.epi <- res.gse51032.df %>%
  filter(chromosome != 0) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  group_by(method, disease, epi_region_id, disease_n) %>%
  summarize(n = n(),
            freq = n()/disease_n) %>%
  ungroup() %>%
  distinct() 

## Colon tumor recurrent
plotDisease(gset_filt, candRegsGR["chr1_84326185"])


recur.disease.epi <- res.gse51032.df %>%
  filter(chromosome != 0) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  mutate(n_total = length(unique(sample))) %>%
  group_by(method, epi_region_id, n_total) %>%
  summarize(n = n(),
            freq = n()/n_total) %>%
  ungroup() %>%
  distinct() 

## Epimutations tumor
plotDisease(gset_filt, candRegsGR["chr3_9404422"])



res.gse51032.shared <- res.gse51032.df %>%
  filter(chromosome != 0) %>%
  group_by(sample, epi_region_id) %>%
  filter(length(unique(method)) == 2)


getMeanDifference <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}
magnitudes <- mclapply(seq_len(nrow(res.gse51032.shared)), 
                                        function(i) 
                                          getMeanDifference(res.gse51032.shared[i, ]$cpg_ids,
                                                            res.gse51032.shared[i, ]$sample, 
                                                            gset_filt), mc.cores = 10)
res.gse51032.shared$magnitude <- unlist(magnitudes)

annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

isTSS <- lapply(res.gse51032.shared$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]
  
  any(sapply(annot[cpgs, ]$UCSC_RefGene_Group, function(i) grepl("TSS", i)))
})
res.gse51032.shared$TSS <- unlist(isTSS)
subset(res.gse51032.shared, abs(magnitude) > 0.6 & TSS) %>%
  arrange(sample, epi_region_id, method) %>% data.frame() %>% head()

plotDisease(gset_filt, GRanges("chr11:64513056-64513809"))



plotDisease <- function(set, range){
  
  miniset <- subsetByOverlaps(set, range)
  
  df <- getBeta(miniset)
  
  df <- t(df) %>% data.frame()
  df$id <- colnames(miniset)
  
  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$disease <- miniset$disease
  df.gath$pos <- start(rowRanges(miniset)[df.gath$cpg])
  
  
  ggplot(df.gath, aes(x = pos, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation", limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}
arrange(recur.disease.epi, desc(freq))

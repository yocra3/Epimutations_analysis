#'#################################################################################
#'#################################################################################
#' Analyze GSE87650 dataset
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)
library(tidyverse)
library(robustbase)

load("results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Remove replicates  ####
gset_filt <- gset[, gset$description == "sample"]

## Compute residuals of pcs  ####
beta <- meffil:::impute.matrix(getBeta(gset_filt), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(gset_filt), probe.range = 40000)
m <- getM(gset_filt)
res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(gset_filt) <- beta
save(gset_filt, file = "results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

case <- gset_filt[, !gset_filt$disease %in% c("HS", "HL")]
control <- gset_filt[, gset_filt$disease %in% c("HS", "HL")]

## Run epimutations ####
methods <- c("beta", "quantile", "mlm")
names(methods) <- methods


res.gse87650.casecontrol.list <- lapply(methods, epimutations, case_samples = case, 
                            control_panel = control)
save(res.gse87650.casecontrol.list, file = "results/epimutations/GSE87650.epimutations.cases.residuals.Rdata")

res.gse87650.loo.list <- lapply(methods, epimutations_one_leave_out, methy = gset_filt, 
                                BPPARAM = MulticoreParam(2))
save(res.gse87650.loo.list, file = "results/epimutations/GSE87650.epimutations.loo.residuals.Rdata")

# Process epimutations ####
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

## Case control ####
### Preprocess data ####
nMethod <- sapply(res.gse87650.casecontrol.list, nrow)
nMethod[sapply(nMethod, is.null)] <- 0
res.gse87650.cc.df <-  Reduce(rbind, res.gse87650.casecontrol.list) %>%
    mutate(method = rep(methods, unlist(nMethod))) %>%
    left_join(colData(gset_filt) %>% 
                data.frame() %>% 
                select(Sample_Name, disease, age, Sex, smoking) %>% 
                mutate(sample = Sample_Name))

res.gse87650.cc.sum <- res.gse87650.cc.df %>%
    group_by(method, sample, disease, age, Sex, smoking) %>%
    summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1", 
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


res.gse87650.cc.sum %>%
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
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

### Differences between UC and CD ####
epi_disease_tabs <- lapply(methods, function(m){
    tab <- res.gse87650.cc.sum %>%
      filter(method %in% m) %>%
      mutate(cat = ifelse(n == 0, "No epi", "epi")) 
    table(tab$cat, tab$disease)
  })
lapply(epi_disease_tabs, fisher.test)

disease.mod <- lapply(methods, function(x) {
  summary(glmrob(n ~ disease, data = res.gse87650.cc.sum, family = "poisson", subset = method == x & n > 0))
})


### Get recurrent epimutations ####
recur.disease.epi <- res.gse87650.cc.df %>%
  filter(chromosome != 0) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  group_by(method, disease, epi_region_id, disease_n) %>%
  summarize(n = n(),
            freq = n()/disease_n) %>%
  ungroup() %>%
  distinct() 
recur.disease.epi %>% select(-n) %>% spread(method, freq) %>% arrange(desc(quantile))

candRegsGR <- epimutacions:::get_candRegsGR()

subset(res.gse87650.cc.df, epi_region_id == "chr1_50879560" & method == "quantile") %>% 
  select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>% data.frame() %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce()

## Crohn recurrent: cerca OTX1
plotDisease(gset_filt, GRanges("chr2:63279495-63286355")) 

## Crohn recurrent: cerca IRX1
plotDisease(gset_filt, GRanges("chr5:3592464-3592638")) 
plotDisease(gset_filt, GRanges("chr5:3599012-3602413"))

## Crohn recurrent: cerca NR2F2
plotDisease(gset_filt, GRanges("chr15:96884949-96888024")) 
plotDisease(gset_filt, GRanges("chr15:96890452-96890880")) 

## Crohn recurrent: cerca NR2F2
plotDisease(gset_filt, GRanges("chr7:96647021-96655889")) 
plotDisease(gset_filt, GRanges("chr7:96650096-96652115")) 
plotDisease(gset_filt, GRanges("chr7:96654782-96655889")) 

## Crohn recurrent: cerca TFAP2A
plotDisease(gset_filt, GRanges("chr6:10381196-10383147")) 
plotDisease(gset_filt, GRanges("chr6:10385320-10385903"))

plotDisease(gset_filt, candRegsGR["chr1_50879560"])

## Crohn recurrent: cerca DMRTA2
plotDisease(gset_filt, GRanges("chr1:50882910-50885352")) 
plotDisease(gset_filt, GRanges("chr1:50886393-50886782"))



recur.ibd.epi <- res.gse87650.cc.df %>%
  filter(chromosome != 0) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  mutate(n_total = length(unique(sample))) %>%
  group_by(method, epi_region_id, n_total) %>%
  summarize(n = n(),
            freq = n()/n_total) %>%
  ungroup() %>%
  distinct() 

recur.ibd.epi %>% select(-n) %>% spread(method, freq) %>% arrange(desc(quantile))
## UC has very few recurrent epimutations and we cannot find common epimutations between both

### Explore individual epimutations ####
res.gse87650.cc.shared <- res.gse87650.cc.df %>%
  filter(chromosome != 0) %>%
  group_by(sample, epi_region_id) %>%
  filter(length(unique(method)) == 3)


getMeanDifference <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}
magnitudes <- mclapply(seq_len(nrow(res.gse87650.cc.shared)), 
                                        function(i) 
                                          getMeanDifference(res.gse87650.cc.shared[i, ]$cpg_ids,
                                                            res.gse87650.cc.shared[i, ]$sample, 
                                                            gset_filt), mc.cores = 10)
res.gse87650.cc.shared$magnitude <- unlist(magnitudes)

annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

isTSS <- lapply(res.gse87650.cc.shared$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]
  
  any(sapply(annot[cpgs, ]$UCSC_RefGene_Group, function(i) grepl("TSS", i)))
})
res.gse87650.cc.shared$TSS <- unlist(isTSS)
subset(res.gse87650.cc.shared, abs(magnitude) > 0.5 & TSS) %>%
  select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(sample, epi_region_id, method) %>% data.frame() 



plotDisease(gset_filt, GRanges("chr12:118405665-118406158")) ## ksr2
plotDisease(gset_filt, GRanges("chr7:27160520-27160674")) ## HOXA3
plotDisease(gset_filt, GRanges("chr12:51717674-51718112")) ## BIN2
plotDisease(gset_filt, GRanges("chr11:18433554-18433745")) ## LDHC
plotDisease(gset_filt, GRanges("chr11:66511804-66512979")) ## C11orf80
plotDisease(gset_filt, GRanges("chr1:1851439-1851910")) ## TMEM52
plotDisease(gset_filt, GRanges("chr6:101846767-101846797")) ## GRIK2
plotDisease(gset_filt, GRanges("chr10:105978651-105978702")) ## WDR6
plotDisease(gset_filt, GRanges("chr11:60738971-60739178")) ## CD6


res.gse87650.cc.shared.tss <- res.gse87650.cc.shared[res.gse87650.cc.shared$TSS, ]
genesTSS <- lapply(res.gse87650.cc.shared.tss$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]
  
  unique(unlist(strsplit(annot[cpgs, ]$UCSC_RefGene_Name, ";")))
})

### Explore genes ####
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "carlos.ruiza@upf.edu", 
  password = "disgenet2806" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

res_df <- Reduce(rbind, lapply(unique(unlist(genesTSS)), function(x) {
  res <- gene2disease(gene = x)
  if (is.character(res)){
    return(NULL)
  }
  extract(res)
}))

disease_genes <- res_df[grep("Crohn|Colitis|Bowel", res_df$disease_name),]$gene_symbol
res.gse87650.cc.shared.tss$TSSgene <- genesTSS
filter(res.gse87650.cc.shared.tss, TSSgene %in% disease_genes & abs(magnitude) > 0.3) %>%
  select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(desc(abs(magnitude))) %>% data.frame()
plotDisease(gset_filt, GRanges("chr2:182321354-182321489")) ## ITGA4 - methylation PMID: 25902909
plotDisease(gset_filt, GRanges("chr11:1874017-1874320")) ## LSP1 - dudoso
plotDisease(gset_filt, GRanges("chr16:68676451-68676743")) ## CDH3 - dudoso
plotDisease(gset_filt, GRanges("chr6:133035150-133035266")) ## VNN1 - expression PMID: 23949622

## Leave-one-out ####
### Preprocess data ####
nMethod <- sapply(res.gse87650.loo.list, nrow)
nMethod[sapply(nMethod, is.null)] <- 0
res.gse87650.loo.df <-  Reduce(rbind, res.gse87650.loo.list) %>%
  mutate(method = rep(methods, unlist(nMethod))) %>%
  left_join(colData(gset_filt) %>% 
              data.frame() %>% 
              select(Sample_Name, disease, age, Sex, smoking) %>% 
              mutate(sample = Sample_Name)) %>%
  mutate(disease = ifelse(disease %in% c("HL", "HS"), "Healthy", disease),
         disease = factor(disease, levels = c("Healthy", "CD", "UC")))

res.gse87650.loo.sum <- res.gse87650.loo.df %>%
  group_by(method, sample, disease, age, Sex, smoking) %>%
  summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n == 1, "1", 
                               ifelse(n < 6, "2-5",
                                      ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


res.gse87650.loo.sum %>%
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
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")


res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm")),
         smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  group_by(method, smoking, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, smoking, disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = smoking, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ disease) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Smoking") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

res.gse87650.loo.sum %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, Sex, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, Sex, disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = Sex, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ disease) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Sex") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

res.gse87650.loo.sum %>%
  filter(n < 20) %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  ggplot(aes(x = age, y = n)) +
  geom_point() +
  theme_bw() +
  facet_grid(method ~ disease) +
  scale_y_continuous(name = "Number of epimutations") +
  scale_x_continuous(name = "Age")


### Factors influencing having an epimutation ####
epi_model_epi_risk <- lapply(methods, function(m){
  tab <- res.gse87650.loo.sum %>%
    filter(method %in% m) %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) 
  glm(out ~ disease + smoking + age + Sex, tab, family = "binomial")
})
lapply(epi_model_epi_risk, summary)
## No assocs

disease.mod.loo <- lapply(methods, function(x) {
  summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = method == x & n > 0))
})
## Higher number of epimutations in CD

disease.mod.loo.adj <- lapply(methods, function(m) {
  tab <- res.gse87650.loo.sum %>%
    filter(method %in% m) %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(n > 0) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current")))
  summary(glmrob(n ~ disease + smoking + age + Sex, data = tab, family = "poisson"))
})
## Higher number of epimutations for CD
disease.mod2.loo <- lapply(methods, function(x) {
  summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = method == x))
})
disease.mod2.loo.adj <- lapply(methods, function(m) {
  tab <- res.gse87650.loo.sum %>%
    filter(method %in% m) %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current")))
  summary(glmrob(n ~ disease + smoking + age + Sex, data = tab, family = "poisson"))
})  



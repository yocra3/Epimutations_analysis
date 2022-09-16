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
library(topGO)
library(pheatmap)
library(cowplot)

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
  df.gath$disease <- ifelse(df.gath$disease %in% c("HL", "HS"), "Healthy", df.gath$disease)
  df.gath$disease <- factor(df.gath$disease, levels = c("Healthy", "CD", "UC"))
  
  df.gath$Coordinates <- start(rowRanges(miniset)[df.gath$cpg])
  
  
  ggplot(df.gath, aes(x = Coordinates, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation level", limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = "Status", 
                       values = c("grey", "green", "blue")) 
    
}

## Case control ####
### Preprocess data ####
res.gse87650.cc.df <-  res.gse87650.casecontrol.list$quantile %>%
    left_join(colData(gset_filt) %>% 
                data.frame() %>% 
                dplyr::select(Sample_Name, disease, age, Sex, smoking) %>% 
                mutate(sample = Sample_Name))

res.gse87650.cc.sum <- res.gse87650.cc.df %>%
    group_by(sample, disease, age, Sex, smoking) %>%
    summarize(n = sum(start != 0 & cpg_n > 2)) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1", 
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))

res.gse87650.cc.sum %>%
  group_by(disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

extremeids.cc <- subset(res.gse87650.cc.sum, n > 20)$sample
extremeids.cc.cd <- subset(res.gse87650.cc.sum, n > 20 & disease == "CD")$sample
extremeids.cc.uc <- subset(res.gse87650.cc.sum, n > 20 & disease == "UC")$sample


### Get recurrent epimutations ####
recur.disease.epi <- res.gse87650.cc.df %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(disease, epi_region_id, disease_n) %>%
  summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct() 
arrange(recur.disease.epi, desc(freq))

gse87650.recur.epi <- filter(res.gse87650.cc.df, epi_region_id %in% subset(recur.disease.epi, n >= 3)$epi_region_id) %>%
  group_by(epi_region_id) %>%
  summarize(samples = paste(sample, collapse = ";"),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            cpgs = paste(unique(unlist(strsplit(cpg_ids, ","))), collapse = ";"),
            cpg_n = length(unique(unlist(strsplit(cpg_ids, ","))))) %>%
  mutate(size = end - start) %>%
  dplyr::select(epi_region_id, samples, chromosome, start, end, size, cpgs, cpg_n)

write.table(gse87650.recur.epi, file = "figures/GSE87650.recurrEpi.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")  

candRegsGR <- epimutacions:::get_candRegsGR()

subset(res.gse87650.cc.df, epi_region_id == "chr6_32793355" & cpg_n > 2) %>% 
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>% data.frame() %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce()

## Crohn recurrent: cerca NR2F2
### Full
plotDisease(gset_filt, GRanges("chr15:96884949-96890880")) 

### Regions
plotDisease(gset_filt, GRanges("chr15:96884949-96888024")) 
plotDisease(gset_filt, GRanges("chr15:96890452-96890880")) 

## Crohn recurrent: cerca DLX5
## Full 
plotDisease(gset_filt, GRanges("chr7:96650096-96655889")) 

## Regions
plotDisease(gset_filt, GRanges("chr7:96650096-96652115")) 
plotDisease(gset_filt, GRanges("chr7:96654782-96655889")) 

## Crohn recurrent: cerca TFAP2A
## Full
plotDisease(gset_filt, GRanges("chr6:10381196-10385903")) 

## Regions
plotDisease(gset_filt, GRanges("chr6:10381196-10383147")) 
plotDisease(gset_filt, GRanges("chr6:10385320-10385903"))


## Crohn recurrent: chr6_32793355 - in HLA region 


## Crohn recurrent: cerca OTX1
plotDisease(gset_filt, GRanges("chr2:63279495-63286355")) 

## Crohn recurrent: cerca IRX1
plotDisease(gset_filt, GRanges("chr5:3592464-3592638")) 
plotDisease(gset_filt, GRanges("chr5:3599012-3602413"))


plotDisease(gset_filt, candRegsGR["chr1_50879560"])

## Crohn recurrent: cerca DMRTA2
plotDisease(gset_filt, GRanges("chr1:50882910-50885352")) 
plotDisease(gset_filt, GRanges("chr1:50886393-50886782"))



recur.ibd.epi <- res.gse87650.cc.df %>%
  filter(chromosome != 0 & cpg_n > 2 & !sample %in% extremeids.cc) %>%
  mutate(n_total = length(unique(sample))) %>%
  group_by(epi_region_id, n_total) %>%
  summarize(n = n(),
            freq = n()/n_total) %>%
  ungroup() %>%
  distinct() 

arrange(recur.ibd.epi, desc(freq))
## UC has very few recurrent epimutations and we cannot find common epimutations between both

#### Explore individuals with recurrent epimutations ####
recu.regs.cd <- subset(recur.disease.epi, n >= 3 & disease == "CD")$epi_region_id
recu.inds.cd <- lapply(recu.regs.cd, function(reg){
  subset(res.gse87650.cc.df, epi_region_id == reg)$sample
})
table(sapply(recu.inds.cd, function(x) mean(x %in% extremeids.cc.cd)))
mean(sapply(recu.inds.cd, function(x) mean(x %in% extremeids.cc.cd)) == 1)

gset_filt$disease2 <- ifelse(gset_filt$disease %in% c("HL", "HS"), "Control", gset_filt$disease)
gset_filt$epi_disease <- ifelse(gset_filt$Sample_Name %in% extremeids.cc.cd, "Epi_CD",
                                ifelse(gset_filt$Sample_Name %in% extremeids.cc.uc, "Epi_UC", gset_filt$disease2))

#### PC whole methylome
pcs_resid <- meffil.methylation.pcs(getBeta(gset_filt), full.obj = TRUE,
                                    probe.range = 40000)
resid.pcs.vars <- pcs_resid$sdev^2/sum(pcs_resid$sdev^2)

pcs.resid.plot <- data.frame(pcs_resid$x, disease = gset_filt$epi_disease) %>%
  mutate(Disease = factor(disease, levels = c("Control", "CD", "Epi_CD", "UC", "Epi_UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Disease)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", round(resid.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(resid.pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(name = "Disease", 
                     values = c("grey", "green", "darkgreen", "cyan", "blue")) +
  ggtitle("Whole Methylation")

#### PC recurr epimutations
rec.cpgs <- subset(res.gse87650.cc.df, epi_region_id %in% recu.regs.cd)$cpg_ids
rec.cpgs <- unique(unlist(strsplit(rec.cpgs, ",")))
save(rec.cpgs, file =  "results/epimutations/GSE87650.recurrent.cpgs.Rdata")

pcs_rec <- meffil.methylation.pcs(getBeta(gset_filt[rec.cpgs, ]), full.obj = TRUE,
                                  probe.range = 40000)
rec.pcs.vars <- pcs_rec$sdev^2/sum(pcs_rec$sdev^2)

pcs.rec.plot <- data.frame(pcs_rec$x, disease = gset_filt$epi_disease) %>%
  mutate(Disease = factor(disease, levels = c("Control", "CD", "Epi_CD", "UC", "Epi_UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Disease)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", round(rec.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(rec.pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Disease", 
                     values = c("grey", "green", "darkgreen", "cyan", "blue")) +
  ggtitle("CpGs in recurrent epimutations")

png("figures/GSE87650.PCs.recurrentEpis.png", width = 700, height = 350)
plot_grid(pcs.resid.plot, pcs.rec.plot, labels = LETTERS[1:2], rel_widths = c(3, 4))
dev.off()


col_colors <- list(
  epi_disease = c("Control" = "black", "CD" = "green", "UC" = "blue",
              "Epi_CD" = "pink", "Epi_UC" = "brown"),
  Sex = c("F" = "purple", "M" = "lightblue")
)
pheatmap(getBeta(gset_filt[rec.cpgs, ]), scale = "none", 
         annotation_col  = data.frame(colData(gset_filt)[, c("epi_disease", "Sex"), drop = FALSE]),
         annotation_colors =  col_colors, 
         show_rownames = FALSE, show_colnames = FALSE)


### Explore individual epimutations ####
getMeanDifference <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}
res.gse87650.cc.reg.df <- subset(res.gse87650.cc.df, chromosome != 0 & cpg_n >= 3)

magnitudes <- mclapply(seq_len(nrow(res.gse87650.cc.reg.df)), 
                                        function(i) 
                                          getMeanDifference(res.gse87650.cc.reg.df[i, ]$cpg_ids,
                                                            res.gse87650.cc.reg.df[i, ]$sample, 
                                                            gset_filt), mc.cores = 10)
res.gse87650.cc.reg.df$magnitude <- unlist(magnitudes)

annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

isTSS <- lapply(res.gse87650.cc.reg.df$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]
  
  any(sapply(annot[cpgs, ]$UCSC_RefGene_Group, function(i) grepl("TSS", i)))
})
res.gse87650.cc.reg.df$TSS <- unlist(isTSS)
subset(res.gse87650.cc.reg.df, abs(magnitude) > 0.5 & TSS & !sample %in% recu.inds.cd & !epi_region_id %in% recu.regs.cd) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(sample, epi_region_id) %>% data.frame() 



plotDisease(gset_filt, GRanges("chr12:118405665-118406158")) ## ksr2
plotDisease(gset_filt, GRanges("chr7:27160520-27160674")) ## HOXA3
plotDisease(gset_filt, GRanges("chr12:51717674-51718112")) ## BIN2
plotDisease(gset_filt, GRanges("chr11:18433554-18433745")) ## LDHC
plotDisease(gset_filt, GRanges("chr11:66511804-66512979")) ## C11orf80
plotDisease(gset_filt, GRanges("chr1:1851439-1851910")) ## TMEM52
plotDisease(gset_filt, GRanges("chr6:101846767-101846797")) ## GRIK2
plotDisease(gset_filt, GRanges("chr10:105978651-105978702")) ## WDR6
plotDisease(gset_filt, GRanges("chr11:60738971-60739178")) ## CD6


res.gse87650.cc.tss <- res.gse87650.cc.reg.df[res.gse87650.cc.reg.df$TSS, ]
genesTSS <- lapply(res.gse87650.cc.tss$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]
  
  unique(unlist(strsplit(annot[cpgs, ]$UCSC_RefGene_Name, ";")))
})

### Explore genes ####
library(disgenet2r)
pass <- ""
disgenet_api_key <- get_disgenet_api_key(
  email = "carlos.ruiza@upf.edu", 
  password = pass )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

res_df <- Reduce(rbind, lapply(unique(unlist(genesTSS)), function(x) {
  res <- gene2disease(gene = x)
  if (is.character(res)){
    return(NULL)
  }
  extract(res)
}))

disease_genes <- res_df[grep("Crohn|Colitis|Bowel", res_df$disease_name),]$gene_symbol
res.gse87650.cc.tss$TSSgene <- genesTSS
filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & !sample %in% extremeids.cc.cd & !epi_region_id %in% recu.regs.cd ) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(Sample_Name) %>% data.frame()


prior.epi1 <- plot_epimutations(filter(res.gse87650.cc.tss, sample == "GSM2337477" & epi_region_id == "chr1_247578552"),
                  gset_filt) + ggtitle("NLRP3 TSS")  ## Gene NLRP3 PMID: 28523573

filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & cpg_n >= 3 & abs(magnitude) > 0.3) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(TSSgene) %>% data.frame() 

prior.epi3 <- plotDisease(gset_filt, GRanges("chr2:182321354-182321489")) +
  ggtitle("ITGA4 TSS") ## ITGA4 - methylation PMID: 25902909

prior.epi2 <- plot_epimutations(filter(res.gse87650.cc.tss, sample == "GSM2337551" & epi_region_id == "chr4_154709441"),
                  gset_filt) + ggtitle("SFRP2 TSS") ## SFRP2 PMID: 18716850

plotDisease(gset_filt, GRanges("chr11:1874017-1874320")) ## LSP1 - dudoso
plotDisease(gset_filt, GRanges("chr16:68676451-68676743")) ## CDH3 - dudoso

png("figures/GSE87650.priorEpimutations.png", width = 484, height = 600)
plot_grid(prior.epi1, prior.epi2, prior.epi3, ncol = 1, labels = LETTERS[1:3])
dev.off()

gse87650.prior.epi <- filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & cpg_n >= 3 & abs(magnitude) > 0.25) %>%
  group_by(epi_region_id) %>%
  summarize(samples = paste(sample, collapse = ";"),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            cpgs = paste(unique(unlist(strsplit(cpg_ids, ","))), collapse = ";"),
            cpg_n = length(unique(unlist(strsplit(cpg_ids, ",")))),
            magnitudes = paste(round(magnitude, 2), collapse = ";"),
            TSSgene = unique(unlist(TSSgene))) %>%
  mutate(size = end - start) %>%
  dplyr::select(epi_region_id, samples, chromosome, start, end, size, cpgs, cpg_n, magnitudes, TSSgene)

write.table(gse87650.prior.epi, file = "figures/GSE87650.prioirizedEpi.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")  

## Leave-one-out ####
### Preprocess data ####
res.gse87650.loo.df <-  res.gse87650.loo.list$quantile %>%
  left_join(colData(gset_filt) %>% 
              data.frame() %>% 
              dplyr::select(Sample_Name, disease, age, Sex, smoking) %>% 
              mutate(sample = Sample_Name)) %>%
  mutate(disease = ifelse(disease %in% c("HL", "HS"), "Healthy", disease),
         disease = factor(disease, levels = c("Healthy", "CD", "UC")))

res.gse87650.loo.sum <- res.gse87650.loo.df %>%
  group_by(sample, disease, age, Sex, smoking) %>%
  summarize(n = sum(!is.na(cpg_n) & cpg_n >= 3)) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n == 1, "1", 
                               ifelse(n < 6, "2-5",
                                      ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


disease.burden <- res.gse87650.loo.sum %>%
  group_by(disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/GSE87650.diseaseburden.png", width = 500, height = 300)
disease.burden
dev.off()

res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  group_by(smoking, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(smoking, disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = smoking, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid( ~ disease) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Smoking") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")
# 
# res.gse87650.loo.sum %>%
#   group_by(Sex, disease, n_cat) %>%
#   summarize(n = n()) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   complete( Sex, disease, n_cat, fill = list(n = 0, p = 0)) %>%
#   ggplot(aes(x = Sex, y = p*100, color = n_cat, fill = n_cat)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   facet_grid( ~ disease) +
#   scale_y_continuous(name = "Proportion of individuals") +
#   scale_x_discrete(name = "Sex") +
#   scale_color_discrete(name = "Epimutations per sample") +
#   scale_fill_discrete(name = "Epimutations per sample")

res.gse87650.loo.sum %>%
  mutate(epimut = ifelse(n_cat == 0, "No Epi", "Epi"),
         epimut = factor(epimut, levels = c("No Epi", "Epi"))) %>%
  ggplot(aes(x = epimut, y = age)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid( ~ disease) +
  scale_x_discrete(name = "Number of epimutations") +
  scale_y_continuous(name = "Age")


### Recurrent epimutations - for comparison with case-control ####
recur.disease.epi.loo <- res.gse87650.loo.df %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(disease, epi_region_id, disease_n) %>%
  summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct() 
arrange(recur.disease.epi.loo, desc(freq))



### Factors influencing having an epimutation ####

### Having epimutation vs not having an epimutation
#### Crude
epi_model_crude_risk <- res.gse87650.loo.sum %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
  glm(out ~ disease, ., family = "binomial") 
summary(epi_model_crude_risk)

#### Adjusted
epi_model_adj_risk <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  mutate(out = ifelse(n > 1, 1, n)) %>%
  glm(out ~ disease + smoking + age + Sex, ., family = "binomial")
summary(epi_model_adj_risk)

#### Stratify by sex
epi_model_sex_risk <- lapply(c("F", "M"), function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(Sex == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glm(out ~ disease + smoking + age, ., family = "binomial")
})
lapply(epi_model_sex_risk, summary)

### Other covariables
epi_model_covars <- lapply(disease, function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(disease == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glm(out ~ Sex + smoking + age, ., family = "binomial")
})
lapply(epi_model_covars, summary)



### Total number of epimutations (including 0)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson"))

#### Adjusted
disease.mod2.loo.adj <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod2.loo.adj)


### Total number of epimutations (excluding 0)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = n > 0))

#### Adjusted
disease.mod.loo.adj <- res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(n > 0) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod.loo.adj)

### Other covariates
epi_num_covars <- lapply(disease, function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(disease == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glmrob(n ~ smoking + age + Sex, data = ., family = "poisson")
})
lapply(epi_num_covars, summary)



### Total number of epimutations (excluding 0 and extreme)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = n > 0 & n < 21))

#### Adjusted
disease.mod.loo.adj2 <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  filter(n > 0 & n < 21) %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod.loo.adj2)


### Explore genes associated with epimutations ####
cpg.annot <- getAnnotation(gset_filt)

disease <- as.character(unique(res.gse87650.loo.df$disease))
names(disease) <- disease
#### Map to genes with Illumina annotation ####
genes_list <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & cpg_n > 2)$cpg_ids
  cpgs_list <- strsplit(cpgs, ",")
  gene_list <- lapply(cpgs_list, function(cpg){
    if (length(cpg) == 1 & is.na(cpg)){
      return(NA)
    }
    genes <- cpg.annot[cpg, "UCSC_RefGene_Name"]
    genes_vec <- unlist(strsplit(genes, ";"))
    genes_vec <- unique(genes_vec[!is.na(genes_vec)])
  })
  gene_list
})
genic_tab <- sapply(genes_list, function(x) table(sapply(x, length) > 0))

## Healthy vs CD
chisq.test(genic_tab[, 1:2])

## Healthy vs UC
chisq.test(genic_tab[, 2:3])

genic_tab %>%
  as_tibble() %>%
  mutate(Gene = c("Intergenic", "Gene")) %>%
  gather(Disease, N, 1:3) %>%
  group_by(Disease) %>%
  mutate(prop = N/sum(N),
         Disease = factor(Disease, levels = c("Healthy", "CD", "UC"))) %>%
  ggplot(aes(x = Disease, y = prop, fill = Gene)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_fill_discrete(name = "Epimutation location")

genes_vec <- lapply(genes_list, function(x){
  unique(unlist(x))
})


all_cpgs <- names(subsetByOverlaps(rowRanges(gset_filt), candRegsGR))
all_genes <- cpg.annot[all_cpgs, "UCSC_RefGene_Name"]
all_genes <- unlist(strsplit(all_genes, ";"))
all_genes <- unique(all_genes[!is.na(all_genes)])

getGOdata <- function(gene_vec){
  genesv <- factor(as.numeric(all_genes %in% gene_vec))
  names(genesv) <- all_genes
  
  Data <- new("topGOdata", 
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.org,
              nodeSize = 10,
              ID = "alias", mapping = "org.Hs.eg")
}
godata <- lapply(genes_vec, getGOdata)
gotests <- lapply(godata, runTest, algorithm = "weight01", statistic = "fisher")
  
finTab <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
                     topNodes = length(score(y)))
}, godata, gotests)

## Repeat tests without epimutations from individuals with extreme number
extremeids <- subset(res.gse87650.loo.sum, n > 20)$sample

genes_list2 <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & !sample %in% extremeids & cpg_n > 2)$cpg_ids
  cpgs_list <- strsplit(cpgs, ",")
  gene_list <- lapply(cpgs_list, function(cpg){
    if (length(cpg) == 1 & is.na(cpg)){
      return(NA)
    }
    genes <- cpg.annot[cpg, "UCSC_RefGene_Name"]
    genes_vec <- unlist(strsplit(genes, ";"))
    genes_vec <- unique(genes_vec[!is.na(genes_vec)])
  })
  gene_list
})
genic_tab2 <- sapply(genes_list2, function(x) table(sapply(x, length) > 0))

## Healthy vs CD
chisq.test(genic_tab2[, 1:2])

## Healthy vs UC
chisq.test(genic_tab2[, 2:3])

genic_tab2 %>%
  as_tibble() %>%
  mutate(Gene = c("Intergenic", "Gene")) %>%
  gather(Disease, N, 1:3) %>%
  group_by(Disease) %>%
  mutate(prop = N/sum(N),
         Disease = factor(Disease, levels = c("Healthy", "CD", "UC"))) %>%
  ggplot(aes(x = Disease, y = prop, fill = Gene)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_fill_discrete(name = "Epimutation location")

genes_vec2 <- lapply(genes_list2, function(x){
  unique(unlist(x))
})

godata2 <- lapply(genes_vec2, getGOdata)
gotests2 <- lapply(godata2, runTest, algorithm = "weight01", statistic = "fisher")

finTab2 <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata2, gotests2)
### Resultados no interesantes

#### Map to genes with eQTMs ####
eqtm <- read_delim("data/eqtm.txt.gz", delim = "\t")

all_tcs <- subset(eqtm, CpG %in% all_cpgs)$TC

getGOdata_eQTM <- function(cpg_vec){
  
  tcs <- subset(eqtm, CpG %in% cpg_vec)$TC
  
  genesv <- factor(as.numeric(all_tcs %in% tcs))
  names(genesv) <- all_tcs
  
  Data <- new("topGOdata", 
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.db,
              nodeSize = 10,
              affyLib = "hta20transcriptcluster.db")
}

cpgs_list <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & cpg_n > 2)$cpg_ids
  cpgs_list <- unlist(strsplit(cpgs, ","))
  cpgs_vec <- unique(cpgs_list[!is.na(cpgs_list)])
})

godata_eqtm <- lapply(cpgs_list, getGOdata_eQTM)
gotests_eqtm <- lapply(godata_eqtm, runTest, algorithm = "weight01", statistic = "fisher")

finTab_eqtm <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata_eqtm, gotests_eqtm)


cpgs_list_filt <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & !sample %in% extremeids & cpg_n > 2)$cpg_ids
  cpgs_list <- unlist(strsplit(cpgs, ","))
  cpgs_vec <- unique(cpgs_list[!is.na(cpgs_list)])
})

godata_eqtm_f <- lapply(cpgs_list_filt, getGOdata_eQTM)
gotests_eqtm_f <- lapply(godata_eqtm_f, runTest, algorithm = "weight01", statistic = "fisher")

finTab_eqtm_f <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata_eqtm_f, gotests_eqtm_f)




finTab_eqtm_df <- cbind(finTab_eqtm_f[[2]][, 1:2], Reduce(cbind, lapply(c(2, 1, 3), function(x){
  df <- finTab_eqtm_f[[x]]
  rownames(df) <- df$GO.ID
  df <- df[finTab_eqtm_f[[2]]$GO.ID, 3:6]
  colnames(df) <- paste0(names(finTab_eqtm_f)[x], colnames(df))
  df
}))) %>%
  filter((CDSignificant > 2 & as.numeric(CDw0) < 0.01) |
           (HealthySignificant > 2 & as.numeric(Healthyw0) < 0.01) |
           (UCSignificant > 2 & as.numeric(UCw0) < 0.01))
write.table(finTab_eqtm_df, file = "figures/GSE87650.loo.GOs.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")  
  
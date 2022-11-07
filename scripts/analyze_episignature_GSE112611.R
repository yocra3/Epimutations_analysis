#'#################################################################################
#'#################################################################################
#' Check episignature in GSE112611
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(tidyverse)
library(meffil)
library(epimutacions)

## Load data
load("data/GSE112611/GSE112611.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gset_GSE112611 <- gset

cellcounts_GSE112611 <- meffil.estimate.cell.counts.from.betas(getBeta(gset_GSE112611), cell.type.reference = "blood gse35069 complete",verbose=T)
cellcounts_GSE112611 <- cellcounts_GSE112611 %>%
  as.data.frame() %>%
  mutate(Disease = gset_GSE112611$`diagnosis:ch1`)

cellcounts_GSE112611 %>%
  gather(CellType, proportion, 1:7) %>%
  ggplot(aes(x = CellType, y = proportion, color = Disease)) +
  geom_boxplot() +
  theme_bw()

boxplot(NK ~ Disease, cellcounts_GSE112611)


load("results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gset_GSE87650 <- gset

gset_comb <- combineArrays (gset_GSE112611, gset_GSE87650)
gset_comb$Dataset <- ifelse(is.na(gset_comb$disease), "GSE112611", "GSE87650")
gset_comb$Disease <- ifelse(is.na(gset_comb$disease), gset_comb$`diagnosis:ch1`, gset_comb$disease)
gset_comb$Disease[gset_comb$Disease == "Crohn's disease"] <- "CD"
gset_comb$Disease[gset_comb$Disease %in% c("non-IBD control", "HS", "HL")] <- "Control"
gset_comb$Time <- ifelse(is.na(gset_comb$`baseline vs follow-up:ch1`), "BL", gset_comb$`baseline vs follow-up:ch1`)
gset_comb$Sex <- ifelse(is.na(gset_comb$Sex), gset_comb$`gender:ch1`, gset_comb$Sex)
gset_comb$Sex[gset_comb$Sex == "F"] <- "Female"
gset_comb$Sex[gset_comb$Sex == "M"] <- "Male"
gset_comb$Age <- ifelse(is.na(gset_comb$age), as.numeric(gset_comb$`age:ch1`), gset_comb$age)

## Check global patterns ####
pcs <- meffil.methylation.pcs(getBeta(gset_comb), full.obj = TRUE)

pcsdf <- data.frame(pcs$x[, 1:10]) %>%
  mutate(diagnosis = gset_comb$Disease,
          dataset = gset_comb$Dataset,
         Time = gset_comb$Time)
ggplot(pcsdf, aes(x = PC1, y = PC2, col = dataset)) + geom_point()

### Very strong batch effect -> apply ComBat
modcombat <- model.matrix(~ Disease + Sex + Age, data = colData(gset_comb))

## Run comBat
m <- getM(gset_comb)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = gset_comb$Dataset, mod = modcombat, par.prior=TRUE, prior.plots = FALSE)
beta <- ilogit2(combat_M)

gset_combat <- gset_comb
assay(gset_combat) <- beta
save(gset, file = "results/preprocess/GSE112611/GSE112611_GSE87650comb.normalizedComBat.GenomicRatioSet.Rdata")

pcs_comb <- meffil.methylation.pcs(getBeta(gset_combat), full.obj = TRUE)

pcsdf_comb <- data.frame(pcs_comb$x[, 1:10]) %>%
  mutate(diagnosis = gset_comb$Disease,
          dataset = gset_comb$Dataset,
          Sex = gset_comb$Sex,
         Time = gset_comb$Time)
ggplot(pcsdf_comb, aes(x = PC1, y = PC2, col = dataset)) + geom_point()
ggplot(pcsdf_comb, aes(x = PC1, y = PC2, col = Sex)) + geom_point()
ggplot(pcsdf_comb, aes(x = PC1, y = PC2, col = Time)) + geom_point()
ggplot(pcsdf_comb, aes(x = PC1, y = PC2, col = diagnosis)) + geom_point()


gsecomb.residuals.cc <- epimutations(gset_comb[, gset_comb$Disease != "Control"], gset_comb[, gset_comb$Disease == "Control"], method = "quantile")
save(gsecomb.residuals.cc, file = "results/epimutations/GSE112611_GSE87650.epimutations.cases.residuals.Rdata")

gset_comb$Sample_Name <- colnames(gset_comb)

res.comb.cc.df <-  gsecomb.residuals.cc %>%
    left_join(colData(gset_comb) %>%
                data.frame() %>%
                dplyr::select(Sample_Name, Disease, Age, Sex, Time, Dataset) %>%
                mutate(sample = Sample_Name))

#
recur.disease_comb.epi <- res.comb.cc.df %>%
  group_by(Disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(Disease, Dataset, epi_region_id, disease_n) %>%
  dplyr::summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
spread(recur.disease_comb.epi, Dataset, n ) %>%
  arrange(desc(GSE87650 )) %>%
  filter(!is.na(GSE112611))

## chr17_76875678 in TIMP-2 with CD -> PMID: 27034654

# save(gse112611.residuals.cc, file = "results/epimutations/GSE112611.epimutations.casecontrol.residuals.Rdata")




vars <- c("sex", "age", "diagnosis", "change", "diag_comb")
lm1 <- lapply(vars, function(x) summary(lm(formula(paste("PC1 ~ ", x)), pcsdf)))
lm2 <- lapply(vars, function(x) summary(lm(formula(paste("PC2 ~ ", x)), pcsdf)))

## Compare epimutations identified in both studies
load("results/epimutations/GSE112611.epimutations.casecontrol.residuals.Rdata")
load("results/epimutations/GSE87650.epimutations.cases.residuals.Rdata")

com.epis <- rbind(mutate(gse112611.residuals.cc , dataset = "GSE112611")[, -13],
                  mutate(res.gse87650.casecontrol.list$quantile , dataset = "GSE87650"))

rec.epis.com <- com.epis %>%
  group_by(epi_region_id, dataset) %>%
  dplyr::summarize(n = length(unique(sample))) %>%
  spread(dataset, n) %>%
  arrange(desc(GSE87650 ))


filter(rec.epis.com, !is.na(GSE87650) & !is.na(GSE112611) & GSE87650 > 1 & GSE112611 > 1)

## Get recurrent Cpgs from GSE87650 ####
load( "results/epimutations/GSE87650.recurrent.cpgs.Rdata")

com.cpgs <- intersect(rec.cpgs, rownames(gset))
pcs.rec <- meffil.methylation.pcs(getBeta(gset[com.cpgs, ]), full.obj = TRUE)


pcs.recdf <- data.frame(pcs.rec$x[, 1:10]) %>%
  mutate(sex = gset$`gender:ch1`,
         age = as.numeric(gset$`age:ch1`),
         diagnosis = gset$`diagnosis:ch1`,
         change = gset$`baseline vs follow-up:ch1`,
         diag_base = gset$`disease stage at diagnosis:ch1`,
         diag_follow = gset$`disease stage at at follow up:ch1`,
         diag_comb = ifelse(diag_base == "NA", diag_follow, diag_base))


ggplot(pcs.recdf, aes(x = PC1, y = PC2, col = sex)) + geom_point()
ggplot(pcs.recdf, aes(x = PC1, y = PC2, col = age)) + geom_point()
ggplot(pcs.recdf, aes(x = PC1, y = PC2, col = diagnosis)) + geom_point()
ggplot(pcs.recdf, aes(x = PC1, y = PC2, col = change)) + geom_point()
ggplot(pcs.recdf, aes(x = PC1, y = PC2, col = diag_comb)) + geom_point()


col_colors <- list(
  `diagnosis:ch1` = c("Crohn's disease" = "black", "non-IBD control" = "green")
)
pheatmap(getBeta(gset[com.cpgs, ]), scale = "none",
         annotation_col  = data.frame(colData(gset)[, "diagnosis:ch1", drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)


## Compute residuals ####
beta <- meffil:::impute.matrix(getBeta(gset), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
m <- getM(gset)
res <- residuals(lmFit(m, pcs$x[, seq_len(ndim)]), m)
beta <- ilogit2(res)

gset_resid <- gset
assay(gset_resid) <- beta

# Check global patterns ####
pcs_resid <- meffil.methylation.pcs(getBeta(gset_resid), full.obj = TRUE)

pcs.residdf <- data.frame(pcs_resid$x[, 1:10]) %>%
  mutate(sex = gset$`gender:ch1`,
         age = as.numeric(gset$`age:ch1`),
         diagnosis = gset$`diagnosis:ch1`,
         change = gset$`baseline vs follow-up:ch1`,
         diag_base = gset$`disease stage at diagnosis:ch1`,
         diag_follow = gset$`disease stage at at follow up:ch1`,
         diag_comb = ifelse(diag_base == "NA", diag_follow, diag_base))


ggplot(pcs.residdf, aes(x = PC1, y = PC2, col = sex)) + geom_point()
ggplot(pcs.residdf, aes(x = PC1, y = PC2, col = age)) + geom_point()
ggplot(pcs.residdf, aes(x = PC1, y = PC2, col = diagnosis)) + geom_point()
ggplot(pcs.residdf, aes(x = PC1, y = PC2, col = change)) + geom_point()
ggplot(pcs.residdf, aes(x = PC1, y = PC2, col = diag_comb)) + geom_point()

## Get recurrent Cpgs from GSE87650 ####
pcs.rec_resid <- meffil.methylation.pcs(getBeta(gset_resid[com.cpgs, ]), full.obj = TRUE)


pcs_resid.recdf <- data.frame(pcs.rec_resid$x[, 1:10]) %>%
  mutate(sex = gset$`gender:ch1`,
         age = as.numeric(gset$`age:ch1`),
         diagnosis = gset$`diagnosis:ch1`,
         change = gset$`baseline vs follow-up:ch1`,
         diag_base = gset$`disease stage at diagnosis:ch1`,
         diag_follow = gset$`disease stage at at follow up:ch1`,
         diag_comb = ifelse(diag_base == "NA", diag_follow, diag_base))


ggplot(pcs_resid.recdf, aes(x = PC1, y = PC2, col = sex)) + geom_point()
ggplot(pcs_resid.recdf, aes(x = PC1, y = PC2, col = age)) + geom_point()
ggplot(pcs_resid.recdf, aes(x = PC1, y = PC2, col = diagnosis)) + geom_point()
ggplot(pcs_resid.recdf, aes(x = PC1, y = PC2, col = change)) + geom_point()
ggplot(pcs_resid.recdf, aes(x = PC1, y = PC2, col = diag_comb)) + geom_point()


pheatmap(getBeta(gset_resid[com.cpgs, ]), scale = "none",
         annotation_col  = data.frame(colData(gset_resid)[, "diagnosis:ch1", drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)

#'#################################################################################
#'#################################################################################
#' Check episignature in GSE112611
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(tidyverse)
library(meffil)

## Load data 
load("data/GSE112611/GSE112611.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Check global patterns #### 
pcs <- meffil.methylation.pcs(getBeta(gset), full.obj = TRUE)

pcsdf <- data.frame(pcs$x[, 1:10]) %>%
  mutate(sex = gset$`gender:ch1`, 
         age = as.numeric(gset$`age:ch1`),
         diagnosis = gset$`diagnosis:ch1`,
         change = gset$`baseline vs follow-up:ch1`,
         diag_base = gset$`disease stage at diagnosis:ch1`,
         diag_follow = gset$`disease stage at at follow up:ch1`, 
         diag_comb = ifelse(diag_base == "NA", diag_follow, diag_base))


ggplot(pcsdf, aes(x = PC1, y = PC2, col = sex)) + geom_point()
ggplot(pcsdf, aes(x = PC1, y = PC2, col = age)) + geom_point()
ggplot(pcsdf, aes(x = PC1, y = PC2, col = diagnosis)) + geom_point()
ggplot(pcsdf, aes(x = PC1, y = PC2, col = change)) + geom_point()
ggplot(pcsdf, aes(x = PC1, y = PC2, col = diag_comb)) + geom_point()


vars <- c("sex", "age", "diagnosis", "change", "diag_comb")
lm1 <- lapply(vars, function(x) summary(lm(formula(paste("PC1 ~ ", x)), pcsdf)))
lm2 <- lapply(vars, function(x) summary(lm(formula(paste("PC2 ~ ", x)), pcsdf)))

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


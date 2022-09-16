#'#################################################################################
#'#################################################################################
#' Explore signatures on inversions in COVID dataset
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(meffil)
library(tidyverse)
library(xgboost)
library(caret)
library(cowplot)
library(readxl)

## Load data
load("data/all_scoreinvhap_final.RData")
load("results/preprocess/GSE168739/GSE168739.normalizedComBat.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

pheno_scourge <- read_xlsx("data/REV_PHENO_EUROPE.xlsx")
index <- read_xlsx("data/CIBERER_index_Carlos Ruiz.xlsx")
colnames(index) <- c("CIBERER", "Interno", "Redcap")

## Map methylation ids to genotypes
cD <- colData(gset)
cD <- left_join(as_tibble(cD), index, by = c("ID3" = "Redcap"))
cD$ID_geno <- gsub("-", "_", cD$CIBERER)
cD <- left_join(cD, all_inv, by = c("ID_geno" = "sample_ID"))
cD$Batch2 <- recode(cD$Batch, "Aurora Pujol (IDIBELL)" = "IDIBELL",
  "Francesc Vidal (Hospital Universitari Joan XXIII de Tarragona)" = "Joan XXIII",
  "Israel Fernandez (Hospital Sant Pau)" = "Sant Pau",
  "Pablo Lapunzina (Hospital Universitario La Paz)" = "La Paz",
  "Pere Soler (Hospital Universitari Vall d’Hebron)" = "Vall d'Hebron")

comb_pheno <- inner_join(cD, pheno_scourge, by = c("ID_geno" = "FID"))

pc <- meffil.methylation.pcs(getBeta(gset), full.obj = TRUE)


## Explore principal components
vars <- round((pc$sdev)**2/sum((pc$sdev)**2)*100, 2)
df.pc <- cbind(as_tibble(pc$x[, 1:5]), tibble(batch = cD$Batch2, disease = gset$Disease, sex = gset$Sex)) %>%
cbind(colData(gset)[, c("CD4T", "Bcell", "CD8T", "Eos", "Mono", "Neu", "NK")]) %>%
  as_tibble()



png("figures/GSE168739/combat_pcs_vs_batch.png", width = 1000)
ggplot(df.pc, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", vars[1], "%)")) +
  scale_y_continuous(name = paste0("PC2 (", vars[2], "%)"))

dev.off()

png("figures/GSE168739/combat_pcs_vs_batch_boxplot.png", width = 1000)
ggplot(df.pc, aes(x = batch, y = PC1)) +
  geom_boxplot() +
  theme_bw()
dev.off()

png("figures/GSE168739/combat_pcs_vs_disease.png", width = 1000)
ggplot(df.pc, aes(x = PC1, y = PC2, color = disease)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", vars[1], "%)")) +
  scale_y_continuous(name = paste0("PC2 (", vars[2], "%)"))

dev.off()


png("figures/GSE168739/combat_pcs_vs_disease_boxplot.png", width = 1000)
ggplot(df.pc, aes(x = disease, y = PC1)) +
    geom_point() +
    theme_bw()
dev.off()



summary(lm(PC1 ~ disease + batch + sex + CD4T + Bcell + CD8T + Eos + Mono + Neu + NK, df.pc))
summary(lm(PC2 ~ disease + batch + sex + CD4T + Bcell + CD8T + Eos + Mono + Neu + NK, df.pc))


## Inversions

invs <- scoreInvHap::inversionGR[c("inv8_001" , "inv17_007", "inv16_009"  )]

pcs <- lapply(seq_len(length(invs)), function(x){
  mini <- subsetByOverlaps(gset, invs[x])
  meffil.methylation.pcs(getBeta(mini))
})

names(pcs) <- c("inv8", "inv17", "inv16")

#Function for creating new model with the CpG sites available
create_model <- function(predict_inv,cpgs){
  #Create DMatrix with the CpGs provided and the labels
  predict_inv$train_mat <-
    predict_inv$data_ref %>%
    select(all_of(cpgs)) %>%
    as.matrix() %>%
    xgb.DMatrix(data = ., label = predict_inv$target)

  #Set parameters
  params <- list(booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss")

  #Create model
  predict_inv$model_new <- xgb.train(params = params, data = predict_inv$train_mat, nrounds = 300)
  return(predict_inv$model_new)
}

#Function for predicting the genotype
predict_geno <- function(predict_inv,model_feat="model_all",dataset,min.cpgs=0.7,inv){

  cat(paste0("INVERSION ",inv," PREDICTION\n\n"))

  ### 1. Select the CpGs intersected

  features <- predict_inv[[model_feat]]$feature_names
  cpgs <- intersect(features, colnames(dataset))
  dataset <- dataset[,cpgs]

  ### 2. Reuse the model if all the CpGs needed are avaibale or create a new model if not

  if(length(cpgs)==length(features)){
    model <- predict_inv[[model_feat]]
    cat("All the CpG sites needed are available in the dataset...\n\n")
    cat("Predicting the genotypes using the reference model...\n\n")
  }

  else {
    if (length(cpgs)>length(features)*min.cpgs){
      cat(paste0("Creating a new model with ",length(cpgs)," CpG sites...\n\n"))
      model <- create_model(predict_inv,cpgs)
    }
    else{
      cat(paste0("There are not enough CpG sites(",length(cpgs),"/",length(features),") for genotyping the inversion.\n\n"))
      return("Not enough Cpg sites.")
    }
  }

  ### 3. Predict the inversion genotypes

  preds <- predict(model, newdata = dataset)

  preds <- matrix(preds, nrow = 3, ncol = length(preds) / 3) %>%
    t() %>%
    data.frame() %>%
    mutate(max = max.col(., ties.method = "last"),
           label=str_replace_all(as.character(max), c("1" = "NN", "2" = "NI", "3" = "II")))
  colnames(preds)[5] <- paste0("label_",inv)
  rownames(preds) <- rownames(dataset)
  return(preds)
}


inv <- paste0("inv", c(17, 8, 16))
names(inv) <- inv

inv_pred <- lapply(inv, function(x){
  load(paste0("data/predict_", x, ".Rdata"))
  preds <- predict_geno(predict_inv, dataset =  t(getBeta(gset)), inv= x)
  out <- preds[, 5]
  names(out) <- rownames(preds)
  out
})

inv_pred_sel <- lapply(inv, function(x){
  load(paste0("data/predict_", x, ".Rdata"))
  preds <- predict_geno(predict_inv, model="model_important", t(getBeta(gset)), inv= x)
  out <- preds[, 5]
  names(out) <- rownames(preds)
  out
})


inv_pred_sel_df <- lapply(inv, function(x){
  load(paste0("data/predict_", x, ".Rdata"))
  preds <- predict_geno(predict_inv, model="model_important", t(getBeta(gset)), inv= x)
})


pc_plots <- lapply(inv, function(i){
  df <- data.frame(pcs[[i]])
  df$inv <- inv_pred[[i]][rownames(df)]

  ggplot(df, aes(x = PC1, y = PC2, color = inv)) +
    geom_point() +
    theme_bw() +
    ggtitle(i)

})
png("figures/GSE168739/inversion_allCpgs.png", width = 1500)
plot_grid(plotlist = pc_plots, nrow = 1)
dev.off()

pc_plots_sel <- lapply(inv, function(i){
  df <- data.frame(pcs[[i]])
  df$inv <- inv_pred_sel[[i]][rownames(df)]

  ggplot(df, aes(x = PC1, y = PC2, color = inv)) +
    geom_point() +
    theme_bw() +
  ggtitle(i)
})

png("figures/GSE168739/inversion_topCpgs.png", width = 1500)
plot_grid(plotlist = pc_plots_sel, nrow = 1)
dev.off()

invs_df <- data.frame(cD[, c("inv8_001", "inv17_007",  "inv16_009" )])
colnames(invs_df) <- c("inv8", "inv17", "inv16")


pc_plots_real <- lapply(inv, function(i){
  df <- data.frame(pcs[[i]])
  df$inv <- invs_df[, i]

  ggplot(df, aes(x = PC1, y = PC2, color = inv)) +
    geom_point() +
    theme_bw() +
  ggtitle(i)
})

png("figures/GSE168739/inversion_scoreInvHap.png", width = 1500)
plot_grid(plotlist = pc_plots_real, nrow = 1)
dev.off()


cD$inv17_pc <- pcs$inv17[, 2]
cD$inv17_num <- ifelse(cD$inv17_pc < -0.3, 2, ifelse(cD$inv17_pc < 0, 1, 0))


summary(glm(factor(Disease) ~ inv17_pc + age + Sex + CD4T + Bcell + CD8T + Eos + Mono + Neu + NK, family = "binomial", cD))
summary(glm(factor(Disease) ~ inv17_num + age + Sex + CD4T + Bcell + CD8T + Eos + Mono + Neu + NK, family = "binomial", cD))

cells <- c("CD4T", "Bcell", "CD8T", "Eos", "Mono", "Neu", "NK")
box_cells <- lapply(cells, function(i){
  ggplot(cD, aes_string(x = "Disease", y = i)) +
    geom_boxplot() +
    theme_bw() +
  ggtitle(i)
})

png("figures/GSE168739/cell_types_disease.png", width = 1500)
plot_grid(plotlist = box_cells, nrow = 2)
dev.off()

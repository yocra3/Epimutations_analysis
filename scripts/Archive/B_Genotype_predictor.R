#######################################
## GENOTYPE PREDICTOR FOR INVERSIONS ##
#######################################

library(xgboost)
library(tidyverse)
library(caret)
library(minfi)

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

#Test in HELIX data using model with all features
for (inv in c("8","16","17")){
  #Load model
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/predict_inv",inv,".Rdata"))
  
  #Load dataset
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy",inv,".Rdata"))
  methy <- eval(parse(text = paste0("methy",inv)))
  
  #Rows: IDs // Columns: CpGs
  df <- t(assay(methy))

  preds <- predict_geno(predict_inv, dataset=df, inv=inv)
  save(preds,file=paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/HELIX/preds",inv,"_all.Rdata"))
}

#Test in HELIX data using model with important features
for (inv in c("8","16","17")){
  #Load model
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/predict_inv",inv,".Rdata"))
  
  #Load dataset
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy",inv,".Rdata"))
  methy <- eval(parse(text = paste0("methy",inv)))
  
  #Rows: IDs // Columns: CpGs
  df <- t(assay(methy))
  
  preds <- predict_geno(predict_inv, model="model_important",df, inv=inv)
  save(preds,file=paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/HELIX/preds",inv,"_important.Rdata"))
}

#Test in TruDiagnostic using model with all features
for (inv in c("8","16","17")){
  #Load model
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/predict_inv",inv,".Rdata"))
  
  #Load dataset
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/TruDiagnostic/methy",inv,"_clean.Rdata"))
  methy <- eval(parse(text = paste0("methy",inv)))
  
  #Rows: IDs // Columns: CpGs
  df <- t(assay(methy))

  preds <- predict_geno(predict_inv, dataset=df, inv=inv)
  save(preds,file=paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/TruDiagnostic/preds",inv,"_all.Rdata"))
}

#Test in TruDiagnostic using model with important features
for (inv in c("8","16","17")){
  #Load model
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/predict_inv",inv,".Rdata"))
  
  #Load dataset
  load(paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/TruDiagnostic/methy",inv,"_clean.Rdata"))
  methy <- eval(parse(text = paste0("methy",inv)))
  
  #Rows: IDs // Columns: CpGs
  df <- t(assay(methy))
  
  preds <- predict_geno(predict_inv, model="model_important",df, inv=inv)
  save(preds,file=paste0("homews/PhD/4_Algorithms_and_functions/Predict_Inversion_Methylation/TruDiagnostic/preds",inv,"_important.Rdata"))
}

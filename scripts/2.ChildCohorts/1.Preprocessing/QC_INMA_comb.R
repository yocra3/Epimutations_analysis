#'#################################################################################
#'#################################################################################
#' Configuration for MeDALL QC of methylation data
#' We will run the normalization including all INMA samples (Esteller and MeDALL)
#' Specifications for some important input files:
#' - phenoPath: Rdata file with a pheno object. It should contain a column called sex, 
#' with sex defined as Female or Male.
#'#################################################################################

## Input files paths
idatsFold <- "data/IDATs_MeDALL"
sampleSheetPattern <- "*.csv"
phenoPath <- "results/preprocess/phenotypes/INMA_phenotypes.Rdata" ## Path to Rdata containing pheno object
genosPath <- "results/preprocess/methylation/INMA.methSNPs.raw"

## Discard samples for the project
discardSamples <- function(samplesheet){
  samplesheet <- samplesheet[!grepl("H2O", samplesheet$Sample_Name), ]
}

## Create common sample ID between IDATs and phenotypes
addSampID <- function(samplesheet) {
  
  ## Modify part of MeDALL
  medall <- subset(samplesheet, is.na(Scan_Date))
  medall$SampleID <- str_extract(medall$Sample_Name, "_[0-9]*_") ## Get internal number
  medall$SampleID <- gsub("_", "", medall$SampleID) ## Remove _
  medall$SampleID <- as.character(as.numeric(medall$SampleID)) ## Remove extra 0s
  medall$Batch <- "MeDALL"
  
  ## Modify part of Esteller
  esteller <- subset(samplesheet, is.na(TEM))
  esteller$SampleID <- as.character(as.numeric(substring(esteller$Sample_Name, 7, 10)))
  esteller$Batch <- "Esteller"
  esteller$Sample_Group <- "Age_0"
  
  ## Merge
  samplesheet <- rbind(esteller, medall)
  
  ## Create variable to correct in ComBat
  samplesheet$TechVar <- ifelse(is.na(samplesheet$TEM), samplesheet$Batch, samplesheet$TEM)
  samplesheet
}

## Create common sample ID between IDATs and genotypes
adaptSampID <- function(geno) {
  geno <- geno[, grep("^SAB", colnames(geno))] ## Select Sabadell samples
  samps <- colnames(geno)
  samps <- gsub("SAB", "", samps) ## Remove SAB prefix
  samps <- as.character(as.numeric(samps)) ## Remove extra 0s
  colnames(geno) <- samps
  geno
}

## Pipeline parameters
cores <- 16

## Variables from phenotype to check batch
batch_var <- c("Slide", "Array",  "TEM", "Sample_Well", "pred.sex", "Sample_Group", "Batch",
               "TechVar")
               
## QC parameters 
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)
pcs <- 15 

## Annotation
array <- "IlluminaHumanMethylation450k"
annotation <- "ilmn12.hg19"

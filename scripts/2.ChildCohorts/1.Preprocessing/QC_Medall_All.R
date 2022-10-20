#'#################################################################################
#'#################################################################################
#' Configuration for MeDALL QC of methylation data
#' We will run the normalization including only samples of Sabadell at 0 years.
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
  samplesheet$SampleID <- str_extract(samplesheet$Sample_Name, "_[0-9]*_") ## Get internal number
  samplesheet$SampleID <- gsub("_", "", samplesheet$SampleID) ## Remove _
  samplesheet$SampleID <- as.character(as.numeric(samplesheet$SampleID)) ## Remove extra 0s
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
batch_var <- c("Slide", "Array",  "TEM", "Sample_Well", "pred.sex", "Sample_Group")
               
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
pcs <- 12 

## Annotation
array <- "IlluminaHumanMethylation450k"
annotation <- "ilmn12.hg19"

#'#################################################################################
#'#################################################################################
#' Configuration for MeDALL QC of methylation data
#' We will run the normalization including only samples of Sabadell at 0 years.
#' Specifications for some important input files:
#' - phenoPath: Rdata file with a pheno object. It should contain a column called sex, 
#' with sex defined as Female or Male.
#'#################################################################################

## Input files paths
idatsFold <- "data/IDATs_Esteller"
sampleSheetPattern <- "*.csv"
phenoPath <- "results/preprocess/phenotypes/INMA_phenotypes.Rdata" ## Path to Rdata containing pheno object
genosPath <- "results/preprocess/methylation/INMA.methSNPs.raw"

## Discard samples for the project
discardSamples <- function(samplesheet){
  samplesheet 
}

## Create common sample ID between IDATs and phenotypes
addSampID <- function(samplesheet) {
  samplesheet$SampleID <- as.character(as.numeric(substring(samplesheet$Sample_Name, 7, 10)))
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
batch_var <- c("Slide", "Array",  "Scan_Date", "Sample_Well", "Sample_Plate", "Sex", 
               "sges", "preterm", "BW", "BL", "HC", "tippart",
               "breastfeeding",  "msmk", "meduc", "edadm", "m_not_eur")
               
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
pcs <- 20

## Annotation
array <- "IlluminaHumanMethylation450k"
annotation <- "ilmn12.hg19"

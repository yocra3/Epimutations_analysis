#'#################################################################################
#'#################################################################################
#' Configuration for MeDALL QC of methylation data
#' Specifications for some important input files:
#' - phenoPath: Rdata file with a pheno object. It should contain a column called sex, 
#' with sex defined as Female or Male.
#'#################################################################################

## Input files paths
idatsFold <- "/DATA/Methylation_INMA/450K_blood/DATA/IDAT/"
sampleSheetPattern <- "*.csv"
phenoPath <- "results/preprocess/phenotypes/INMA_phenotypes.Rdata" ## Path to Rdata containing pheno object
genosPath <- "results/preprocess/methylation/INMA.methSNPs.raw"

## Output name
out <- "MeDALL.normalizedRaw.GenomicRatioSet.Rdata"

## Create common sample ID
addSampID <- function(samplesheet) {
  samplesheet$SampleID <- str_extract(samplesheet$Sample_Name, "_[0-9]*_") 
  samplesheet
}

## Pipeline parameters
cores <- 1

## Discard samples
badSamples <- c()

## Variables from phenotype to check batch
batch_var <- c("Slide", "Array", "Sex", "GestAge", "Status", "pathGroup", "MolecularCause", "SampleBatch", "Rearrangements", "Del22q11", "N_genes")

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
pcs <- 5 

## Annotation
array <- "IlluminaHumanMethylation450k"
annotation <- "ilmn12.hg19"

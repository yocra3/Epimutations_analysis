#'#################################################################################
#'#################################################################################
#' Run methylation QC data
#' This code performs quality control to methylation data. 
#' Important decisions:
#' - Remove samples not belonging to the project.
#' - Remove samples based on QC
#' - Use values from meffil vignette in all parameters
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
configFile <- args[1]

## Load libraries ####
library(meffil)
library(minfi)
library(readxl)
library(ggplot2)
library(dplyr)

## Load parameters
source(configFile)

## Set number of cores
options(mc.cores = cores)


## Prepare sample sheet ####
### Load predefined sample sheet
samplesheet <- meffil.read.samplesheet(base = idatsFold, pattern = sampleSheetPattern)

### Load and adapt samples data
load(phenoPath)

## Merge both (Check in other datasets)
samplesheet <- addSampID(samplesheet)   
combSheet <- left_join(select(samplesheet, -Sex), pheno, by = "SampleID")

## Change sex to M / F
combSheet <- mutate(combSheet, Sex = substring(Sex, 1, 1))

## Discard some samples
combSheet <- subset(combSheet, !Sample_Name %in% badSamples)


## Generate QC report
### Load genotypes
genos <- meffil.extract.genotypes(genosPath)
genotypes <- genos[, match(combSheet$SampleID, colnames(genos))]

## Load methylation data and QC ####
qc.objects <- meffil.qc(combSheet, verbose = TRUE)
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
  genotypes = genotypes
)

outs <- c()
## Remove bad samples based on QC report and rerun QC
outlier <- qc.summary$bad.samples
round <- 1
save(qc.objects, file = paste0("qc.objects.round", round, ".Rdata"))
save(qc.summary, file = paste0("qcsummary.round", round, ".Rdata"))


while (!is.null(outlier)){
  outs <- c(outs, outlier)
  round <- round + 1
  qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
  save(qc.objects, file = paste0("qc.objects.round", round, ".Rdata"))
  
  qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
  save(qc.summary, file = paste0("qcsummary.round", round, ".Rdata"))
}
save(qc.objects, file = "qc.objects.clean.Rdata")
save(qc.summary, file = "qcsummary.clean.Rdata")

## Report filtered samples and probes
write.table(outs, file = "removed.samples.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = "removed.probes.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")

## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = "pc.fit.pdf", height = 6, width = 6)

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)
save(norm.objects, file = "norm.obj.pc.Rdata")

## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(pcs),
  batch.pcs = seq_len(pcs),
  batch.threshold = 0.01
)
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name, verbose = TRUE)
save(norm.beta, file = "norm.beta.Rdata")

## Check covariables
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 8e5)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = "norm.summary.Rdata")

## Create GenomicRatioSet
rownames(combSheet) <- combSheet$SampleID
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = combSheet[colnames(norm.beta), ],
                                      array = array,
                                      annotation = annotation)
save(gset, file = out)
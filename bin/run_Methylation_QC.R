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
library(ggplot2)
library(dplyr)
library(stringr)

## Load parameters
source(configFile)

## Set number of cores
options(mc.cores = cores)


## Prepare sample sheet ####
### Load predefined sample sheet
samplesheet <- meffil.read.samplesheet(base = idatsFold, pattern = sampleSheetPattern)

## Discard some samples
samplesheet <- discardSamples(samplesheet)

### Load and adapt samples data
load(phenoPath)

## Merge both (Check in other datasets)
samplesheet <- addSampID(samplesheet)   
combSheet <- left_join(select(samplesheet, -Sex), pheno, by = "SampleID")

## Change sex to M / F
combSheet <- mutate(combSheet, Sex = substring(Sex, 1, 1))

## Generate QC report
### Load genotypes
genos <- meffil.extract.genotypes(genosPath)
genos <- adaptSampID(genos)

## Map genotype IDs to IDAT IDs
comIDs <- intersect(combSheet$SampleID, colnames(genos))
combSheet.filt <- subset(combSheet, SampleID %in% comIDs)
genotypes <- genos[, combSheet.filt$SampleID]
colnames(genotypes) <- combSheet.filt$Sample_Name

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
save(qc.objects, file = paste0(outPrefix, ".qc.objects.round", round, ".Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.raw.html"))


while (nrow(outlier)> 0){
  outs <- rbind(outs, outlier)
  round <- round + 1
  qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
  save(qc.objects, file = paste0(outPrefix,"qc.objects.round", round, ".Rdata"))
  
  qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
  save(qc.summary, file = paste0(outPrefix, "qcsummary.round", round, ".Rdata"))
  outlier <- qc.summary$bad.samples
}
save(qc.objects, file = paste0(outPrefix, ".qc.objects.clean.Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.clean.html"))

## Report filtered samples and probes
write.table(outs, file = paste0(outPrefix, ".removed.samples.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = paste0(outPrefix, ".removed.probes.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")

## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = paste0(outPrefix, ".pc.fit.pdf"), height = 6, width = 6)

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)
save(norm.objects, file = paste0(outPrefix, ".norm.obj.pc.Rdata"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name, verbose = TRUE)

## Check covariables
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 40000)
## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(pcs),
  batch.pcs = seq_len(pcs),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = paste0(outPrefix, ".norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.normalization.html"))

## Create GenomicRatioSet
rownames(combSheet) <- combSheet$Sample_Name
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = combSheet[colnames(norm.beta), ],
                                      array = array,
                                      annotation = annotation)
save(gset, file = paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata"))




#'#################################################################################
#'#################################################################################
#' Preprocess CD datasets from GEO
#'#################################################################################
#'#################################################################################

## Uncompress files
tar -xvf data/GSE32148_RAW.tar  -C data/
gunzip -d data/GSE32148/*.idat.gz

mkdir data/GSE81961
tar -xvf data/GSE81961_RAW.tar  -C data/GSE81961/
gunzip -d data/GSE81961/*.idat.gz

tar -xvf data/GSE105798_RAW.tar  -C data/
mkdir data/GSE105798
mv data/*idat.gz data/GSE105798
mv data/GPL* data/GSE105798
gunzip -d data/GSE105798/*.idat.gz

# Preprocess DNA methylation ####
### Try to preprocess the different arrays simoultaneously
## Load libraries
library(GEOquery)
library(minfi)
library(meffil)
library(tidyverse)

options(timeout = 10000)
options(mc.cores = 5)

## Load phenotypes
gse81961.full <- getGEO("GSE81961", GSEMatrix = TRUE, getGPL = FALSE)

pheno  <- data.frame(Sample_Name = colnames(gse81961.full[[1]]),
                     Disease = gse81961.full[[1]]$`disease state:ch1`,
                     Sex =  gse81961.full[[1]]$`gender:ch1`, 
                     Age = as.numeric(gse81961.full[[1]]$`age (yr):ch1`), 
                     atnf = gse81961.full[[1]]$`atnf usage at time of phlebotomy:ch1`)

pheno$Disease <- ifelse(pheno$Disease %in% c("UC", "ulcerative colitis"), "UC",
  ifelse(pheno$Disease %in% c("Crohn", "CD", "Crohn's disease"), "CD", "Control") )
pheno$Sex[pheno$Sex == "Female"] <- "F"
pheno$Sex[pheno$Sex == "Male"] <- "M"
rownames(pheno) <- pheno$Sample_Name

## Create samplesheet
samplesheet <- meffil.create.samplesheet(path = "data/GSE81961/")
samplesheet$Sex <- pheno[samplesheet$Sample_Name, "Sex"]
samplesheet$Disease <- pheno[samplesheet$Sample_Name, "Disease"]
samplesheet$Age <- pheno[samplesheet$Sample_Name, "Age"]
samplesheet$atnf <- pheno[samplesheet$Sample_Name, "atnf"]

## Remove samples without phenotypes
samplesheet <- subset(samplesheet, !is.na(Sex))

qc.objects <- meffil.qc(samplesheet, verbose = TRUE,  cell.type.reference="blood gse35069 complete")
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1,
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
)

## Define function to select outliers
filterOutliers <- function(outlier){
  subset(outlier, issue %in% c("Methylated vs Unmethylated",
                               "X-Y ratio outlier",
                               "Low bead numbers",
                               "Detection p-value",
                               "Sex mismatch",
                               "Genotype mismatch")
  )
}

## Remove bad samples based on QC report and rerun QC
outlier <- qc.summary$bad.samples
# outlier <- filterOutliers(outlier)
round <- 1

outPrefix <- "results/preprocess/GSE81961/GSE81961"
save(qc.objects, file = paste0(outPrefix, ".qc.objects.round", round, ".Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.raw.html"))

## Report filtered samples and probes
write.table(qc.summary$bad.cpgs, file = paste0(outPrefix, ".removed.probes.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")


## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = paste0(outPrefix, ".pc.fit.pdf"), height = 6, width = 6)


pcs <- 4
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)

## Add predicted sex as sample sheet variable
for (i in seq_len(length(norm.objects))){
  norm.objects[[i]]$samplesheet$pred.sex <- norm.objects[[i]]$predicted.sex
}
save(norm.objects, file = paste0(outPrefix, ".norm.obj.pc.Rdata"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name,
                                      verbose = TRUE)
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 40000)

batch_var <- c("Slide", "sentrix_col",  "sentrix_row", "Sex", "Disease", "Age", "atnf")

norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(8),
  batch.pcs = seq_len(8),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = paste0(outPrefix, ".norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.normalization.html"))

rownames(samplesheet) <- samplesheet$Sample_Name
samplesheet.final <- samplesheet[colnames(norm.beta), ]
## Add predicted sex
samplesheet.final$pred.sex <- vapply(norm.objects, function(x) x$predicted.sex, character(1))

## Add cell counts
cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc <- data.frame(IID=row.names(cc),cc)
samplesheet.final <- cbind(samplesheet.final, cc[rownames(samplesheet.final), ])

## Save genomicratioset
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = samplesheet.final,
                                      array = "IlluminaHumanMethylation450k",
                                      annotation = "ilmn12.hg19")
save(gset, file = paste0(outPrefix, ".GenomicRatioSet.Rdata"))



## Load annotation ####
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

## Create Initial and final dataset with missings
dp <- meffil.load.detection.pvalues(qc.objects)
dp.f <- dp[rownames(gset), colnames(gset)]

beta <- assay(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

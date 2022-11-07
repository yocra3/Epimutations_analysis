#'#################################################################################
#'#################################################################################
#' Preprocess GSE32148
#'#################################################################################
#'#################################################################################

## Load libraries
library(GEOquery)
library(minfi)

# GSE32148 ####
gse32148.full <- getGEO("GSE32148", GSEMatrix = TRUE, getGPL = FALSE)

gset <- makeGenomicRatioSetFromMatrix(exprs(gse32148.full[[1]]), pData = pData(gse32148.full[[1]]))

outPrefix <- "results/preprocess/GSE32148/GSE32148"

gset$Disease <- gset$`disease state:ch1`
gset$Disease <- ifelse(gset$Disease %in% c("UC", "ulcerative colitis"), "UC",
                             ifelse(gset$Disease %in% c("Crohn", "CD", "Crohn's disease"), "CD", "Control") )

gset$Sex <- gset$`gender:ch1`
gset$Age <- as.numeric(gset$`age (y):ch1`)
gset$Twin <- gset$`twin:ch1`

save(gset, file = paste0(outPrefix, "withNA.GenomicRatioSet.Rdata"))


## Load annotation ####
grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".withNA.allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".withNA.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".withNA.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

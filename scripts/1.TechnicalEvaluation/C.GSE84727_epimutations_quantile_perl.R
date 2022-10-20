#'#################################################################################
#'#################################################################################
#' Prepare GSE84727 to run on original quantile-perl implementation
#'#################################################################################
#'#################################################################################

#'#################################################################################
# Prepare data in R
### Load libraries
library(minfi)
library(BiocParallel)
library(epimutacions)
library(readxl)

## Load data
load("data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Create file
annot <- data.frame(rowRanges(gset))
out <- cbind(rownames(gset), annot, getBeta(gset))
write.table(out, file = "results/barbosa_original/GSE84727.input.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
#'#################################################################################


#'#################################################################################
## Run detection in perl (script found in the same folder than this script)
perl barbosa_paper2.pl -f results/barbosa_original/GSE84727.input.txt \
-c 1 -s 2 -e 3 -a 6 -1 0.005 -9 0.995 -t 0.15 -p 3 \
-o results/barbosa_original/GSE84727.barbosa_result.txt 
#'#################################################################################
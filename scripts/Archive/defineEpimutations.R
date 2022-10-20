#'#################################################################################
#'#################################################################################
#' Define all possible epimutations for 450K
#'#################################################################################
#'#################################################################################

## Load libraries
#' Steps:
#' 1. Make a dataset with two samples. This dataset contains all the CpGs from the array
#' and one samples has all 0s and the other all 1s.
#' 2. Run bumphunter using the same configuration than in the package.
#' 3. Get the total number of regions computed. This is total 
#' number of effective tests for methods with p-value.

# Load libraries
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(rtracklayer)

data(Locations)

### Select CpGs (names starting by cg) in autosomic chromosomes
locs.450 <- subset(Locations, grepl("^cg", rownames(Locations)) & chr %in% paste0("chr", 1:22))
locs.450GR <- makeGRangesFromDataFrame(locs.450, start.field = "pos", end.field = "pos", strand = "*")
locs.450GR <- sort(locs.450GR)
mat <- matrix(0, nrow = length(locs.450GR), ncol = 2, 
              dimnames = list(names(locs.450GR), c("A", "B")))

## Set sample B to all 1
mat[, 2] <- 1

## Define model matrix
pheno <- data.frame(var = c(0, 1))
model <- model.matrix(~ var, pheno)

## Run bumphunter
bumps <- bumphunter(mat, design = model, pos = start(locs.450GR), 
                    chr = as.character(seqnames(locs.450GR)),
                    cutoff = 0.05, maxGap = 1000)$table
bumps.fil <- subset(bumps, L >= 3)
candRegsGR <- makeGRangesFromDataFrame(bumps.fil, keep.extra.columns = TRUE)
export.bed(candRegsGR, con = "results/candidateRegions/candidateRegions.450K.hg19.bed")

## Load ranges in hg38 remapped with NCBI
hg38 <- read.delim("results/candidateRegions/candidateRegions.450K.hg38_liftOver.tab")

## Remove regions with more than 1 mapping
feat_tab <- table(hg38$X.feat_name)
sel_regs <- names(feat_tab)[feat_tab == 1]
hg38.filt <- subset(hg38, X.feat_name %in% sel_regs)

## Select regions with coverage of 1
hg38.sel <- subset(hg38.filt, coverage == "1")

## Rename candRegs to map to hg38
names(candRegsGR) <- paste(seqnames(candRegsGR), start(candRegsGR), sep = "_")
save(candRegsGR, file = "results/candidateRegions/candidateRegions.450K.GRanges.Rdata")

candRegsGR38 <- makeGRangesFromDataFrame(hg38.sel, seqnames.field = "mapped_id",
                                         start.field = "mapped_start", 
                                         end.field = "mapped_stop")
names(candRegsGR38) <- hg38.sel$X.feat_name
mcols(candRegsGR38) <- mcols(candRegsGR)[names(candRegsGR38), ]
save(candRegsGR38, file = "results/candidateRegions/candidateRegions.450K_hg38.GRanges.Rdata")


## Overlap with CREs
CREs <- read.table("data/GRCh38-ccREs.bed")
colnames(CREs) <- c("chr", "start", "end", "var", "Accession", "Type")
CREsGR <- makeGRangesFromDataFrame(CREs, keep.extra.columns = TRUE)
names(CREsGR) <- CREsGR$Accession

creOver <- findOverlaps(candRegsGR38, CREsGR)
candRegsGR38$CRE <- sapply(seq_len(length(candRegsGR38)), function(i){
  paste(names(CREsGR)[to(creOver)[from(creOver) == i]], collapse = ";")
})

candRegsGR38$CRE_type <- sapply(seq_len(length(candRegsGR38)), function(i){
  paste(CREsGR$Type[to(creOver)[from(creOver) == i]], collapse = ";")
})
save(candRegsGR38, file = "results/candidateRegions/candidateRegions.450K_hg38.withCREs.GRanges.Rdata")

df <- data.frame(reg = names(candRegsGR38), CRE = candRegsGR38$CRE, 
                 CRE_type = candRegsGR38$CRE_type)
rownames(df) <- names(candRegsGR38)
candRegsGR$CRE <- df[names(candRegsGR), ]$CRE
candRegsGR$CRE_type <- df[names(candRegsGR), ]$CRE_type
candRegsGR$CRE[is.na(candRegsGR$CRE)] <- ""
candRegsGR$CRE_type[is.na(candRegsGR$CRE_type)] <- ""

save(candRegsGR, file = "results/candidateRegions/candidateRegions.450K.withCREs.GRanges.Rdata")

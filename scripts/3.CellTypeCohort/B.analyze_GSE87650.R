#'#################################################################################
#'#################################################################################
#' Analyze GSE87650 dataset
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)
library(tidyverse)
library(robustbase)
library(UpSetR)
library(cowplot)

load("results/preprocess/GSE87650/GSE87650.cell..autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")

## Remove replicates  ####
gset_filt <- gset[, gset$description == "sample"]

## Check PCAs ####
pcs <- meffil.methylation.pcs(getBeta(gset_filt), probe.range = 40000)
ggplot(data.frame(pcs, cells = gset_filt$`cell type:ch1`), aes(x = PC1, y = PC2, color = cells)) +
  geom_point()

pcsdf <- data.frame(pcs, cells = gset_filt$`cell type:ch1`)
pcsdf$shape <- ifelse(pcsdf$PC1 > -20 & pcsdf$cells == "monocytes", "outliers",
                      ifelse(pcsdf$PC1 < -20 & pcsdf$cells == "CD8", "outliers", "good"))
ggplot(pcsdf, aes(x = PC1, y = PC2, color = shape, shape = cells)) +
  geom_point()

## No differences between the CD datasets
gset_cds <- gset_filt[, gset_filt$`cell type:ch1` %in% c("CD4", "CD8")]
pcs2 <- meffil.methylation.pcs(getBeta(gset_cds), 
                               probe.range = 40000)
ggplot(data.frame(pcs2, cells = gset_cds$`cell type:ch1`), aes(x = PC1, y = PC2, color = cells)) +
  geom_point()

## Remove sample placed on another cluster 
gset_filt <- gset_filt[, pcsdf$shape == "good"]
save(gset_filt, file = "results/preprocess/GSE87650/GSE87650.autosomic.filterAnnotatedProbes.filterIndividuals.GenomicRatioSet.Rdata")

## Split gset by cell type
cells <- unique(gset_filt$`cell type:ch1`)
names(cells) <- cells
gset_list <- lapply(cells, function(x) gset_filt[, gset_filt$`cell type:ch1` == x])

## Compute residuals of pcs  ####
getResiduals <- function(grs){
  beta <- meffil:::impute.matrix(getBeta(grs), margin = 1)
  ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
  pcs <- meffil.methylation.pcs(getBeta(grs), probe.range = 40000)
  m <- getM(grs)
  res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
  beta <- ilogit2(res)
  assay(grs) <- beta
  grs
}
gset_res_list <- lapply(gset_list, getResiduals)

save(gset_res_list, file = "results/preprocess/GSE87650/GSE87650.cellSplit.autosomic.filterAnnotatedProbes.PCAresiduals.GenomicRatioSet.Rdata")

getPCs <- function(grs){
  beta <- meffil:::impute.matrix(getBeta(grs), margin = 1)
  ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
  ndim
}
sapply(gset_list, getPCs)

## Run epimutations ####
methods <- c("beta", "quantile", "mlm")
names(methods) <- methods


res.gse87650.cells.casecontrol.list <- mclapply(gset_res_list, function(set) {
  case <- set[, set$`simplified_diagnosis:ch1` != "HC"]
  control <- set[, set$`simplified_diagnosis:ch1` == "HC"]
  
  lapply(methods, epimutations, case_samples = case, 
         control_panel = control)
  }, mc.cores = 4)
save(res.gse87650.cells.casecontrol.list, file = "results/epimutations/GSE87650.cellsplit.epimutations.cases.residuals.Rdata")

res.gse87650.cells.loo.list <- mclapply(gset_res_list, function(set) {
  lapply(methods, epimutations_one_leave_out, methy = set, 
         BPPARAM = MulticoreParam(2))
}, mc.cores = 4)
save(res.gse87650.cells.loo.list, file = "results/epimutations/GSE87650.cellsplit.epimutations.loo.residuals.Rdata")

# Process epimutations ####
plotDisease <- function(set, range){
  
  miniset <- subsetByOverlaps(set, range)
  
  df <- getBeta(miniset)
  
  df <- t(df) %>% data.frame()
  df$id <- colnames(miniset)
  
  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$disease <- miniset$`simplified_diagnosis:ch1`
  df.gath$cell <- miniset$`cell type:ch1`
  df.gath$pos <- start(rowRanges(miniset)[df.gath$cpg])
  
  
  ggplot(df.gath, aes(x = pos, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation", limits = c(0, 1)) +
    theme_bw() +
    facet_wrap(~ cell) +
    theme(plot.title = element_text(hjust = 0.5))
}

## Leave-one-out ####
### Preprocess data ####
res.gse87650.cell.loo.list2 <- lapply(res.gse87650.cells.loo.list, function(x) x$quantile)

res.gse87650.cell.loo.df <-  Reduce(rbind, res.gse87650.cell.loo.list2) %>%
  mutate(cell = rep(names(res.gse87650.cell.loo.list2), 
                    sapply(res.gse87650.cell.loo.list2, nrow))) %>%
  left_join(colData(gset_filt) %>% 
              data.frame() %>% 
              mutate(sample = geo_accession,
                     disease = .data[["simplified_diagnosis.ch1"]],
                     id = .data[["patient_number.ch1"]],
                     age = as.numeric(.data[["age.at.sample.ch1"]]),
                     sex = .data[["Sex.ch1"]]) %>%
              select(sample, disease, age, sex, id)) 

### Explore epimutations replicability in different tissues ####
res.gse87650.cell.loo.filt <- res.gse87650.cell.loo.df %>%
  filter(chromosome != 0 & cpg_n >= 3) %>%
  mutate(reg_id = paste(epi_region_id, id)) %>%
  group_by(id, epi_region_id) %>%
  mutate(n_cells = length(unique(cell)))

upset.list <- lapply(cells, function(y) {
  subset(res.gse87650.cell.loo.filt, cell == y)$reg_id
})


res.gse87650.cell.loo.props <- res.gse87650.cell.loo.df %>%
  filter(chromosome != 0 & cpg_n >= 3) %>%
  group_by(id, epi_region_id) %>%
  summarize(cells = paste(sort(unique(cell)), collapse = ","))


## Plot some epimutations ####
getMeanDifference <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}
magnitudes <- mclapply(seq_len(nrow(res.gse87650.cell.loo.filt)), 
                       function(i) 
                         getMeanDifference(res.gse87650.cell.loo.filt[i, ]$cpg_ids,
                                           res.gse87650.cell.loo.filt[i, ]$sample, 
                                           gset_res_list[[res.gse87650.cell.loo.filt[i, ]$cell]]), 
                       mc.cores = 10)
res.gse87650.cell.loo.filt$magnitude <- unlist(magnitudes)
subset(res.gse87650.cell.loo.filt, n_cells == 4 & abs(magnitude) > 0.4) %>%
  arrange(reg_id) %>% select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%  data.frame()

res.gse87650.cell.loo.cpgs <- res.gse87650.cell.loo.filt %>%
  select(cell, id, cpg_ids, epi_region_id, reg_id, chromosome, start, end) %>%
  group_by(cell, id, epi_region_id) %>%
  summarize(cpg_ids = paste(unique(cpg_ids), collapse = ",")) %>%
  spread(cell, cpg_ids)
res.gse87650.cell.loo.cpgs$cpg_list <- lapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), function(i) {
  cpg_string <- paste(res.gse87650.cell.loo.cpgs[i, c("CD4", "CD8", "monocytes", "wh blood")], collapse = ",")
  cpg_list <- strsplit(cpg_string, ",")[[1]]
  unique(cpg_list[cpg_list != "NA"])
})

res.gse87650.cell.loo.cpgs$intersect_cpg_list <- lapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), function(i) {
  
  cpg_list <- as.list(res.gse87650.cell.loo.cpgs[i, c("CD4", "CD8", "monocytes", "wh blood")])
  cpg_list <- cpg_list[!sapply(cpg_list, function(x) length(x) == 1 & is.na(x))]
  out <- Reduce(intersect, lapply(cpg_list, function(x) strsplit(x, ",")[[1]]))
  out 
  })

## Remove regions without overlap between the cells
res.gse87650.cell.loo.cpgs <- subset(res.gse87650.cell.loo.cpgs, 
                                     lengths(res.gse87650.cell.loo.cpgs$intersect_cpg_list) >= 3)

getRep <- function(i, col, set) {
  if (!is.na(res.gse87650.cell.loo.cpgs[i, col])){
    return(0)
  } else{
    getMeanQuantile(res.gse87650.cell.loo.cpgs[i, ]$cpg_list[[1]], res.gse87650.cell.loo.cpgs[i, ]$id, set) 
  }
}
getMeanQuantile <- function(cpgs, idnum, set){
  samps <- colnames(set[, set$`patient_number:ch1` == idnum])
  cpgs <- cpgs[cpgs %in% rownames(set)]
  betas <- getBeta(set[cpgs, ])
  quant <- apply(betas, 1, function(x) {
    f <- ecdf(x)
    f(x[colnames(betas) == samps])
  })
  
  mean(quant)
}
res.gse87650.cell.loo.cpgs$rep_quant_CD4 <- sapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), getRep, 
                                    col = "CD4", set = gset_res_list$CD4)
res.gse87650.cell.loo.cpgs$rep_quant_CD8 <- sapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), getRep, 
                                                   col = "CD8", set = gset_res_list$CD8)
res.gse87650.cell.loo.cpgs$rep_quant_mono <- sapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), getRep, 
                                                   col = "monocytes", set = gset_res_list$monocytes)
res.gse87650.cell.loo.cpgs$rep_quant_whole <- sapply(seq_len(nrow(res.gse87650.cell.loo.cpgs)), getRep, 
                                                   col = "wh blood", set = gset_res_list$`wh blood`)

isSig <- function(x) x < 0.04 | x > 0.96


res.gse87650.cell.loo.rep <- res.gse87650.cell.loo.cpgs %>%
  ungroup() %>%
  filter(!is.na(rep_quant_CD4) & !is.na(rep_quant_CD8) & !is.na(rep_quant_mono) & !is.na(rep_quant_whole)) %>%
  mutate(cd4_st = ifelse(!is.na(CD4), "CD4", ""),
         cd8_st = ifelse(!is.na(CD8), "CD8", ""),
         mono_st = ifelse(!is.na(monocytes), "monocytes", ""),
         whole_st = ifelse(!is.na(`wh blood`), "wh blood", ""),
         sigDatasets = paste(cd4_st, cd8_st, mono_st, whole_st, sep = "-"),
         logcd4 = isSig(rep_quant_CD4),
         logcd8 = isSig(rep_quant_CD8),
         logmono = isSig(rep_quant_mono),
         logwhole = isSig(rep_quant_whole),
         cd4_st2 = ifelse(logcd4, "CD4", ""),
         cd8_st2 = ifelse(logcd8, "CD8", ""),
         mono_st2 = ifelse(logmono, "monocytes", ""),
         whole_st2 = ifelse(logwhole, "wh blood", ""),
         sigDatasets2 = paste(cd4_st2, cd8_st2, mono_st2, whole_st2, sep = "-"))


rep.plot <- res.gse87650.cell.loo.rep %>%
  select(id, epi_region_id, sigDatasets, sigDatasets2) %>%
  gather("Type", "CellTypes", 3:4) %>%
  mutate(Type = ifelse(Type == "sigDatasets", "Epimutations", "Epimutations + outliers"),
         simpDatasets = gsub("-", "", CellTypes),
         cells = ifelse(simpDatasets == "CD4CD8monocyteswh blood", "Whole blood + 3 cell types", 
                        ifelse(simpDatasets %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets, 
                               ifelse(grepl("wh", simpDatasets), "Whole blood + 1-2 cell types",  "2+ cell types"))),
         cells = ifelse(cells == "wh blood", "Whole blood", cells),
         cells = factor(cells, levels = c("CD4", "CD8", "monocytes",  "2+ cell types", "Whole blood", 
                                          "Whole blood + 1-2 cell types", 
                                          "Whole blood + 3 cell types"))) %>%
  group_by(cells, Type) %>%
  summarize(n = n()) %>%
  group_by(Type) %>%
  mutate(p = n/sum(n)) %>%
  ggplot(aes( x = Type, y = p*100, fill = cells)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Proportion of epimutations") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Epimutations replicability") +
  scale_fill_manual(name = "", values = c("pink", "darkolivegreen", "turquoise4", "brown", "grey80", "grey50", "grey20"))

save(rep.plot, file = "results/GSE87650.replicability_plot.Rdata")

## Figure 7
png("figures/GSE87650.replicateEpis.cells.png", width = 2400, height = 1400, res = 300)
rep.plot
dev.off()

rep.df <- res.gse87650.cell.loo.rep %>%
  select(id, epi_region_id, sigDatasets, sigDatasets2) %>%
  gather("Type", "CellTypes", 3:4) %>%
  mutate(Type = ifelse(Type == "sigDatasets", "Epimutations", "Epimutations + outliers"),
         simpDatasets = gsub("-", "", CellTypes),
         cells = ifelse(simpDatasets == "CD4CD8monocyteswh blood", "Whole blood + 3 cell types", 
                        ifelse(simpDatasets %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets, 
                               ifelse(grepl("wh", simpDatasets), "Whole blood + 1-2 cell types",  "2+ cell types"))),
         cells = ifelse(cells == "wh blood", "Whole blood", cells),
         cells = factor(cells, levels = c("CD4", "CD8", "monocytes",  "2+ cell types", "Whole blood", 
                                          "Whole blood + 1-2 cell types", 
                                          "Whole blood + 3 cell types"))) %>%
  group_by(cells, Type) %>%
  summarize(n = n()) %>%
  group_by(Type) %>%
  mutate(p = n/sum(n)) 

plotCase <- function(reg_name, res, samp_name){
  
  selRows <- subset(res, reg_id == reg_name)
  
  avcells <- unique(selRows$cell)
  
  if (length(avcells) != 4){
    
    df <- selRows[, c("sample", "chromosome", "start", "end", "cpg_ids", "cell"), drop = FALSE]
    newcells <- setdiff(cells, avcells)
    selRows <- rbind(df, data.frame(sample = sapply(gset_res_list[newcells], function(set){
      colnames(set)[set$`patient_number:ch1` == selRows$id[1]]
    }), 
    chromosome = df$chromosome[1], 
    start = min(df$start), 
    end = max(df$end), cpg_ids = df$cpg_ids[1], 
    cell = newcells))
  }
  
  plots <- lapply(cells, function(selcell) {
    p <- plot_epimutations(subset(selRows, cell == selcell), methy = gset_res_list[[selcell]]) +
    ggtitle(selcell) + scale_color_manual(labels = c("control", samp_name, "mean"), 
                                            values = c("black", "red", "darkblue"))
    if (selcell == "wh blood"){
      p <- p + ggtitle("Whole blood")
    }
    p
  })
  
  plots
}
ex1 <- plotCase("chr12_10095902 9005", res.gse87650.cell.loo.filt)
plot_grid(plotlist = ex1, ncol = 2 )

## Epimutation in 4 cell types
ex2 <- plotCase("chr20_62886664 8975", res.gse87650.cell.loo.filt, "Samp_8975")

## Sup Figure 25
png("figures/GSE87650.replicateEpicell.png", width = 4000, height = 2000, res = 300)
plot_grid(plotlist = ex2[c(3, 1, 2, 4)], ncol = 2 )
dev.off()



## Epimutation in monocytes + whole blood
subset(res.gse87650.cell.loo.rep, sigDatasets2  == "--monocytes-wh blood" & method == "quantile") %>% 
  data.frame()
ex3 <- plotCase("chr10_133793442 8918", "quantile", res.gse87650.cell.loo.filt)
plot_grid(plotlist = ex3, ncol = 2 )

## Epimutation in whole blood
subset(res.gse87650.cell.loo.rep, sigDatasets2  == "---wh blood" & method == "quantile") %>% 
  data.frame()
ex4 <- plotCase("chr10_133793442 8918", "quantile", res.gse87650.cell.loo.filt)
plot_grid(plotlist = ex3, ncol = 2 )


## Epimutation in CD8
subset(res.gse87650.cell.loo.rep, sigDatasets2  == "-CD8--") %>% 
  data.frame()
ex5 <- plotCase("chr11_59823993 9034", res.gse87650.cell.loo.filt, "Samp_9034")

## Sup Figure 26
png("figures/GSE87650.Epicell_specific.png", width = 4000, height = 2000, res = 300)
plot_grid(plotlist = ex5[c(3, 1, 2, 4)], ncol = 2 )
dev.off()


nowh_cpgs <- subset(res.gse87650.cell.loo.rep, !grepl("wh", sigDatasets2)) %>% 
  data.frame() %>% pull("cpg_list") %>% unlist() %>% unique()

pcs_nowh <- meffil.methylation.pcs(getBeta(gset_filt[cd8_cpgs, ]))
pheatmap::pheatmap(getBeta(gset_filt[nowh_cpgs, ]),
                           annotation_col  = data.frame(colData(gset_filt)[, "cell type:ch1", drop = FALSE]))


sp_epis <- subset(res.gse87650.cell.loo.rep, sigDatasets2 %in% c("-CD8--", "--monocytes-"))

diffMedian <-  function(i){
  cpgs <- unlist(sp_epis[i, ]$cpg_list)
  
  cells <- c("CD4", "CD8", "monocytes")[c(sp_epis$logcd4[i], sp_epis$logcd8[i], sp_epis$logmono[i])]
  whmed <- rowMedians(getBeta(gset_filt[cpgs, gset_filt$`cell type:ch1` == "wh blood"]))
  cellmed <- rowMedians(getBeta(gset_filt[cpgs, gset_filt$`cell type:ch1` %in% cells]))
  mean(whmed - cellmed)
}
sp_epis$diffMedians <- sapply(seq_len(nrow(sp_epis)), diffMedian)
mean(abs(sp_epis$diffMedians) > 0.2)

## Correlation with gene expression ####
### Functions ####
getGenesZ <- function(epi_df, row, gexp_list, gexpRange, window){
  
  gexp <- assay(gexp_list[[epi_df[row, "cell"][[1]]]])
  
  exp_name <- epi_df[row, ][["id"]]
  if (!exp_name %in% colnames(gexp)){
    return(  list(geneName = NA,
                  Zscore = NA,
                  outliers = NA,
                  distance = NA, 
                  rank = NA))
  }
  epi_range <- GRanges(gsub("_", ":", epi_df[row, ][["epi_region_id"]]))
  
  transRanges <- subsetByOverlaps(gexpRange, epi_range + window)
  if (length(transRanges) == 0){
    return(  list(geneName = NA,
                  Zscore = NA,
                  outliers = NA,
                  distance = NA, 
                  rank = NA))
  }
  
  selTranscripts <- names(transRanges)
  
  
  genesmat <- gexp[selTranscripts, , drop = FALSE]
  zmat <- apply(genesmat, 1, function(x) (x - mean(x))/sd(x))
  zs <- zmat[exp_name, ]
  
  rankMat <- apply(genesmat, 1, rank)
  ranks <- rankMat[exp_name, ]
  

  isOutliersRow <- function(x){
    qs <- quantile(x, c(0.25, 0.75))
    iq <- qs[2] - qs[1]
    out <- x < qs[1] - 1.5*iq | x > qs[2] + 1.5*iq
    out
  }
  outmat <- apply(genesmat, 1, isOutliersRow)
  outs <- outmat[exp_name, ]

  list(geneNames = selTranscripts,
       Zscores = zs,
       outliers = outs, 
       distances = distance(epi_range, transRanges),
       rank = ranks)
  
}
selectNearest <- function(x) {
  if (length(x$Zscore) == 1){
    data.frame(geneName = x$geneName, Zscore = x$Zscore, outliers = x$outliers,
               distance = x$Zscore, rank = x$rank) 
  }
  else {
    ind <- which.min(x$distance)
    data.frame(geneName =  x$geneName[ind], outliers = x$outliers[ind],
               Zscore = x$Zscore[ind], 
               distance = x$distance[ind], rank = x$rank[ind])
  }
}

isOutliers <- function(x, val){
  qs <- quantile(x, c(0.25, 0.75))
  iq <- qs[2] - qs[1]
  out <- val < qs[1] - 1.5*iq | val > qs[2] + 1.5*iq
  out
}

getGenesTSS <- function(epi_df, row, gexp_list){
  
  gene <- unlist(epi_df[row, "tss_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  selTranscripts <- rownames(gexp)[rowData(gexp)$Symbol %in% gene]
  if (length(selTranscripts) == 0){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  exp_name <- epi_df[row, ][["id"]]
  
  gexp <- assay(gexp_list[[epi_df[row, "cell"][[1]]]])
  
  if (!exp_name %in% colnames(gexp)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  genesmat <- gexp[selTranscripts, , drop = FALSE]
  
  out <- getVals(genesmat, exp_name, selTranscripts)
}

getVals <- function(genesmat, exp_name, selTranscripts){
  zmat <- t(apply(genesmat, 1, function(x) (x - mean(x))/sd(x)))
  zs <- zmat[, exp_name]
  z <- zs[which.max(abs(zs))]
  
  genevec <- genesmat[which.max(abs(zs)), ]
  out <- isOutliers(genevec, genevec[exp_name])
  
  rank <- rank(genevec)[exp_name]
  
  list(z = z, out = out, rank = rank, gene = selTranscripts[which.max(abs(zs))])
  
}

getGeneseQTM <- function(epi_df, row, gexp_list){
  
  gene <- unlist(epi_df[row, "eqtm_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  selTranscripts <- rownames(gexp)[rowData(gexp)$Symbol %in% gene]
  if (length(selTranscripts) == 0){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  exp_name <- epi_df[row, ][["id"]]
  
  gexp <- assay(gexp_list[[epi_df[row, "cell"][[1]]]])
  
  if (!exp_name %in% colnames(gexp)){
    return(   list(z = NA, out = NA, rank = NA, gene = NA))
  }
  genesmat <- gexp[selTranscripts, , drop = FALSE]
  
  out <- getVals(genesmat, exp_name, selTranscripts)
}


### Prepare data ####
load("results/preprocess/GSE87650/GSE87650.gene_exp.SummarizedExperiment.Rdata")
eqtm <- read_delim("data/eqtm.txt.gz", delim = "\t")

### Create z-score matrices
se$cell_type <- se$`cell type:ch1`

## Remove probes without mapping
gexp <- se[rowData(se)$Probe_Coordinates != "", se$`cell type:ch1` != "n4"]
rowData(gexp)$Coords_clean <- gsub("-[0-9]*:[0-9]*", "", rowData(gexp)$Probe_Coordinates, perl = TRUE)
rowData(gexp)$start <- sapply(strsplit(rowData(gexp)$Coords_clean , "-"), 
                              function(x) min(as.numeric(x)))
rowData(gexp)$end <- sapply(strsplit(rowData(gexp)$Coords_clean , "-"),
                            function(x) max(as.numeric(x)))
rowData(gexp)$Chromosome <- paste0("chr", rowData(gexp)$Chromosome) 
gexp_ranges <- makeGRangesFromDataFrame(rowData(gexp))

gexp$cell_type <- ifelse(gexp$cell_type == "CD14", "monocytes", 
                         ifelse(gexp$cell_type == "WB", "wh blood", gexp$cell_type))

gexp.list <- lapply(cells, function(cell) {
  minise <- gexp[, gexp$cell_type == cell]
  colnames(minise) <- minise$`samplenumber:ch1`
  minise
})

### Get TSS ####
cpg.annot <- getAnnotation(gset_res_list[[1]])

tss_genes <- mclapply(res.gse87650.cell.loo.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  annot <- cpg.annot[cpgs, ]
  tss <- annot[grepl("TSS", annot$UCSC_RefGene_Group), , drop = FALSE]
  if (nrow(tss) == 0){
    return(NA)
  } else {
    gene <- strsplit(tss$UCSC_RefGene_Name, ";")
    unique(unlist(gene))
  }
}, mc.cores = 5)

res.gse87650.cell.loo.filt$tss_genes <- tss_genes
tss_list <- lapply(seq_len(nrow(res.gse87650.cell.loo.filt)), getGenesTSS, 
                   epi_df = res.gse87650.cell.loo.filt, gexp_list = gexp.list)

res.gse87650.cell.loo.filt$tss_z <- sapply(tss_list, function(x) x$z)
res.gse87650.cell.loo.filt$tss_out <- sapply(tss_list, function(x) x$out)
res.gse87650.cell.loo.filt$tss_rank <- sapply(tss_list, function(x) x$rank)
res.gse87650.cell.loo.filt$tss_gene <- sapply(tss_list, function(x) x$gene)

### Get eQTM ####
eqtm_genes <- mclapply(res.gse87650.cell.loo.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  tab <- subset(eqtm, CpG %in% cpgs)
  if (nrow(tab) == 0){
    return(NA)
  } else {
    unique(tab$TC_gene )  
  }
}, mc.cores = 5)

res.gse87650.cell.loo.filt$eqtm_genes <- eqtm_genes
eqtm_list <- lapply(seq_len(nrow(res.gse87650.cell.loo.filt)), getGeneseQTM, 
                    epi_df = res.gse87650.cell.loo.filt, gexp_list = gexp.list)
res.gse87650.cell.loo.filt$eqtm_z <- sapply(eqtm_list, function(x) x$z)
res.gse87650.cell.loo.filt$eqtm_out <- sapply(eqtm_list, function(x) x$out)
res.gse87650.cell.loo.filt$eqtm_rank <- sapply(eqtm_list, function(x) x$rank)
res.gse87650.cell.loo.filt$eqtm_gene <- sapply(eqtm_list, function(x) x$gene)

### Close gene ####
geneInfo <- lapply(seq_len(nrow(res.gse87650.cell.loo.filt)), getGenesZ, 
                   epi_df = res.gse87650.cell.loo.filt,
                   gexp_list = gexp.list, gexpRange = gexp_ranges, window = 250e3)
geneInfoNear <- Reduce(rbind, lapply(geneInfo, selectNearest))
colnames(geneInfoNear) <- paste0(colnames(geneInfoNear), "Near" )
res.gse87650.cell.loo.filt$near_z <- geneInfoNear$Zscore
res.gse87650.cell.loo.filt$near_out <- geneInfoNear$outliers
res.gse87650.cell.loo.filt$near_rank <- geneInfoNear$rank
res.gse87650.cell.loo.filt$near_gene <- geneInfoNear$geneName

### Make plots ####
res.gse87650.cell.loo.out.sum <- res.gse87650.cell.loo.filt %>%
  ungroup() %>%
  select(cell, magnitude, reg_id, ends_with("out")) %>%
  gather(Measure, logical, 4:6) %>%
  filter(!is.na(logical)) %>%
  mutate(Measure = ifelse(Measure == "eqtm_out", "eQTM",
                          ifelse(Measure == "near_out", "Near gene", "TSS")),
         Measure = factor(Measure, levels = c("eQTM", "TSS", "Near gene")))
  
gse87650.gexp.outliers.plot <- res.gse87650.cell.loo.out.sum %>%
  mutate(cell = ifelse(cell == "wh blood", "Whole blood", cell)) %>%
  group_by(cell, Measure) %>%
  summarize(p = mean(logical)) %>%
  ggplot(aes(x = cell, y = p*100)) +
    geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5)) + 
  scale_y_continuous(name = "Proportion of outliers", limits = c(0, 50)) +
  scale_x_discrete(name = "Tissue") +
  facet_wrap(Measure ~ .) +
  scale_fill_discrete(name = "")

## Sup Figure 27
png("figures/GSE87650.genexp.outliers.png", height = 1200, width = 1920, res = 300)
gse87650.gexp.outliers.plot
dev.off()


res.gse87650.cell.loo.out.sum %>%
  mutate(Expression = ifelse(logical, "Outlier", "Not Outlier")) %>%
  ggplot(aes(x = Expression, y = abs(magnitude))) +
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(name = "Epimutations magnitude") +
  scale_x_discrete(name = "Gene status") +
  facet_grid(cell ~ Measure)
### Lack of pairs
  

res.gse87650.cell.loo.gexp.sum  <- res.gse87650.cell.loo.filt %>% 
  ungroup() %>%
  select(cell, magnitude, reg_id, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 4:9) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = ifelse(exp_type == "eqtm", "eQTM",
                                 ifelse(exp_type == "near", "Near gene", "TSS")),
         exp_type = factor(exp_type, levels = c("eQTM", "TSS", "Near gene")),
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")))


gse87650.gexp.plots <- res.gse87650.cell.loo.gexp.sum %>% 
  mutate(cell = ifelse(cell == "wh blood", "Whole blood", cell)) %>%
  ggplot(aes(x = exp_type, y = value, color = cell)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5, 
               binwidth = 0.2) +
  theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    facet_grid(measure ~ cell, scales = "free") +
  scale_x_discrete(name = "Gene mapping")

## Sup Figure 28
png("figures/GSE87650.genexp.scores.png", height = 1920, width = 2800, res = 300)
gse87650.gexp.plots
dev.off()
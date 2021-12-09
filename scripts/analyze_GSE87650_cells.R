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

## Remove sample places on another cluster 
gset_filt <- gset_filt[, pcsdf$shape == "good"]

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

## Case control ####
### Preprocess data ####
res.gse87650.cell.casecontrol.list2 <- lapply(res.gse87650.cells.casecontrol.list, function(x){
  nMethod <- sapply(x, nrow)
  nMethod[sapply(nMethod, is.null)] <- 0
  res.gse87650.cc.df <-  Reduce(rbind, x) %>%
    mutate(method = rep(methods, unlist(nMethod)))
})

res.gse87650.cell.cc.df <-  Reduce(rbind, res.gse87650.cell.casecontrol.list2) %>%
    mutate(cell = rep(names(res.gse87650.cell.casecontrol.list2), 
                      sapply(res.gse87650.cell.casecontrol.list2, nrow))) %>%
    left_join(colData(gset_filt) %>% 
                data.frame() %>% 
                mutate(sample = geo_accession,
                       disease = .data[["simplified_diagnosis.ch1"]],
                       id = .data[["patient_number.ch1"]],
                       age = as.numeric(.data[["age.at.sample.ch1"]]),
                       sex = .data[["Sex.ch1"]]) %>%
                select(sample, disease, age, sex, id)) 

res.gse87650.cell.cc.sum <- res.gse87650.cell.cc.df %>%
    group_by(method, sample, disease, age, sex, id, cell) %>%
    summarize(n = sum(start != 0 & (is.na(pvalue) | pvalue < 0.05/40408))) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1", 
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


res.gse87650.cell.cc.sum %>%
  mutate(method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, cell, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(method, disease, cell, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(method ~ cell) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")


### Plot recurrent epimutations identified in whole blood ####
## Crohn recurrent: cerca OTX1 -- revisar
plotDisease(gset_filt, GRanges("chr2:63279495-63286355")) 

## Crohn recurrent: cerca IRX1 - no signal
plotDisease(gset_filt, GRanges("chr5:3592464-3592638")) 
plotDisease(gset_filt, GRanges("chr5:3599012-3602413"))

## Crohn recurrent: cerca NR2F2 - no signal
plotDisease(gset_filt, GRanges("chr15:96884949-96888024")) 
plotDisease(gset_filt, GRanges("chr15:96890452-96890880")) 

## Crohn recurrent: cerca NR2F2 - no signal
plotDisease(gset_filt, GRanges("chr7:96647021-96655889")) 
plotDisease(gset_filt, GRanges("chr7:96650096-96652115")) 
plotDisease(gset_filt, GRanges("chr7:96654782-96655889")) 

## Crohn recurrent: cerca TFAP2A
plotDisease(gset_filt, GRanges("chr6:10381196-10383147")) 
plotDisease(gset_filt, GRanges("chr6:10385320-10385903"))

## Crohn recurrent: cerca DMRTA2 - no signal
plotDisease(gset_filt, GRanges("chr1:50882910-50885352")) 
plotDisease(gset_filt, GRanges("chr1:50886393-50886782"))

## individual epimtuations - no signal in cells 
plotDisease(gset_filt, GRanges("chr12:118405665-118406158")) ## ksr2
plotDisease(gset_filt, GRanges("chr7:27160520-27160674")) ## HOXA3
plotDisease(gset_filt, GRanges("chr12:51717674-51718112")) ## BIN2
plotDisease(gset_filt, GRanges("chr11:18433554-18433745")) ## LDHC
plotDisease(gset_filt, GRanges("chr11:66511804-66512979")) ## C11orf80
plotDisease(gset_filt, GRanges("chr1:1851439-1851910")) ## TMEM52
plotDisease(gset_filt, GRanges("chr6:101846767-101846797")) ## GRIK2
plotDisease(gset_filt, GRanges("chr10:105978651-105978702")) ## WDR6
plotDisease(gset_filt, GRanges("chr11:60738971-60739178")) ## CD6

plotDisease(gset_filt, GRanges("chr2:182321354-182321489")) ## ITGA4 - methylation PMID: 25902909
plotDisease(gset_filt, GRanges("chr11:1874017-1874320")) ## LSP1 - dudoso
plotDisease(gset_filt, GRanges("chr16:68676451-68676743")) ## CDH3 - dudoso
plotDisease(gset_filt, GRanges("chr6:133035150-133035266")) ## VNN1 - expression PMID: 23949622


## Leave-one-out ####
### Preprocess data ####
res.gse87650.cell.loo.list2 <- lapply(res.gse87650.cells.loo.list, function(x){
  nMethod <- sapply(x, nrow)
  nMethod[sapply(nMethod, is.null)] <- 0
  res.gse87650.cc.df <-  Reduce(rbind, x) %>%
    mutate(method = rep(methods, unlist(nMethod)))
})

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
  filter(chromosome != 0) %>%
  mutate(reg_id = paste(epi_region_id, id)) %>%
  group_by(id, method, epi_region_id) %>%
  mutate(n_cells = length(unique(cell)))

upset.list <- lapply(methods, function(x) {
  lapply(cells, function(y) {
  subset(res.gse87650.cell.loo.filt, method == x & cell == y)$reg_id
  })
})

x <- "beta"
png(paste0("figures/GSE87650.", x, ".cell_overlaps.png"), width = 700, height = 600)
upset(fromList(upset.list[[x]]), sets = cells, order.by = "freq",
      mainbar.y.label = "Common epimutations", 
      sets.x.label = "Epimutations per cell")
dev.off() 

x <- "quantile"
png(paste0("figures/GSE87650.", x, ".cell_overlaps.png"), width = 700, height = 600)
upset(fromList(upset.list[[x]]), sets = cells, order.by = "freq",
      mainbar.y.label = "Common epimutations", 
      sets.x.label = "Epimutations per cell")
dev.off() 

x <- "mlm"
png(paste0("figures/GSE87650.", x, ".cell_overlaps.png"), width = 700, height = 600)
upset(fromList(upset.list[[x]]), sets = cells, order.by = "freq",
      mainbar.y.label = "Common epimutations", 
      sets.x.label = "Epimutations per cell")
dev.off() 

res.gse87650.cell.loo.props <- res.gse87650.cell.loo.df %>%
  filter(chromosome != 0) %>%
  group_by(id, method, epi_region_id) %>%
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
subset(res.gse87650.cell.loo.filt, n_cells == 4 & abs(magnitude) > 0.4 & method == "quantile") %>%
  arrange(reg_id) %>% select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%  data.frame()

res.gse87650.cell.loo.cpgs <- res.gse87650.cell.loo.filt %>%
  select(cell, id, method, cpg_ids, epi_region_id, reg_id, chromosome, start, end) %>%
  group_by(cell, method, id, epi_region_id) %>%
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


rep.plot1 <- res.gse87650.cell.loo.rep %>%
  mutate(simpDatasets = gsub("-", "", sigDatasets),
         cells = ifelse(simpDatasets == "CD4CD8monocyteswh blood", "4 tissues", 
                        ifelse(simpDatasets %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets, "2-3 tissues")),
         cells = factor(cells, levels = c("CD4", "CD8", "monocytes", "wh blood", "2-3 tissues", "4 tissues")),
         method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, cells) %>%
  summarize(n = n()) %>%
  group_by(method) %>%
  mutate(p = n/sum(n)) %>%
  ggplot(aes( x = method, y = p*100, fill = cells)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Proportion of epimutations") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Detected epimutations")

rep.plot2 <- res.gse87650.cell.loo.rep %>%
  mutate(simpDatasets = gsub("-", "", sigDatasets2),
         cells = ifelse(simpDatasets == "CD4CD8monocyteswh blood", "4 tissues", 
                        ifelse(simpDatasets %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets, "2-3 tissues")),
         cells = factor(cells, levels = c("CD4", "CD8", "monocytes", "wh blood", "2-3 tissues", "4 tissues")),
         method = factor(method, levels = c("quantile", "beta", "mlm"))) %>%
  group_by(method, cells) %>%
  summarize(n = n()) %>%
  group_by(method) %>%
  mutate(p = n/sum(n)) %>%
  ggplot(aes( x = method, y = p*100, fill = cells)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Proportion of epimutations") +
  theme_bw()  +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Detected epimutations + outliers")

plot_grid(rep.plot1, rep.plot2)


plotCase <- function(reg_name, met, res){
  
  selRows <- subset(res, reg_id == reg_name & method == met)
  
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
    plot_epimutations(subset(selRows, cell == selcell), methy = gset_res_list[[selcell]]) +
    ggtitle(selcell)
  })
  
  plots
}
ex1 <- plotCase("chr12_10095902 9005", "quantile", res.gse87650.cell.loo.filt)
plot_grid(plotlist = ex1, ncol = 2 )

## Epimutation in 4 cell types
ex2 <- plotCase("chr20_62886664 8975", "quantile", res.gse87650.cell.loo.filt)
plot_grid(plotlist = ex2, ncol = 2 )

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

## Correlation with gene expression ####
### Functions ####
getGenesZ <- function(epi_df, row, gexpZ_list, gexpRange, window){
  
  gexpZ <- gexpZ_list[[epi_df[row, "cell"][[1]]]]
  
  exp_name <- epi_df[row, ][["id"]]
  if (!exp_name %in% colnames(gexpZ)){
    return(  list(geneName = NA,
                  Zscore = NA,
                  distance = NA, 
                  rank = NA))
  }
  epi_range <- GRanges(gsub("_", ":", epi_df[row, ][["epi_region_id"]]))
  
  transRanges <- subsetByOverlaps(gexpRange, epi_range + window)
  if (length(transRanges) == 0){
    return(  list(geneName = NA,
                  Zscore = NA,
                  distance = NA, 
                  rank = NA))
  }
  
  selTranscripts <- names(transRanges)
  
  zs <- gexpZ[selTranscripts, exp_name]
  rankMat <- apply(gexpZ[selTranscripts, , drop = FALSE], 1, rank)
  rownames(rankMat) <- colnames(gexpZ)
  ranks <- rankMat[exp_name, ]
  
  
  list(geneNames = selTranscripts,
       Zscores = zs,
       distances = distance(epi_range, transRanges),
       rank = ranks)
  
}
selectNearest <- function(x) {
  if (length(x$Zscore) == 1){
    data.frame(geneName = x$geneName, Zscore = x$Zscore, distance = x$Zscore, rank = x$rank) 
  }
  else {
    ind <- which.min(x$distance)
    data.frame(geneName =  x$geneName[ind], Zscore = x$Zscore[ind], 
               distance = x$distance[ind], rank = x$rank[ind])
  }
}

getGenesTSS <- function(epi_df, row, gexpZ_list){
  
  gene <- unlist(epi_df[row, "tss_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  selTranscripts <- rownames(se)[rowData(se)$Symbol %in% gene]
  if (length(selTranscripts) == 0){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  exp_name <- epi_df[row, ][["id"]]
  
  gexpZ <- gexpZ_list[[epi_df[row, "cell"][[1]]]]
  
  if (!exp_name %in% colnames(gexpZ)){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  zs <- gexpZ[selTranscripts, exp_name]
  z <- zs[which.max(abs(zs))]
  
  
  rankMat <- apply(gexpZ[selTranscripts, , drop = FALSE], 1, rank)
  rownames(rankMat) <- colnames(gexpZ)
  
  ranks <- rankMat[exp_name, ]
  rank <- ranks[which.max(abs(zs))]
  
  list(z = z, rank = rank, gene = selTranscripts[which.max(abs(zs))])
}

getGeneseQTM <- function(epi_df, row, gexpZ_list){
  
  gene <- unlist(epi_df[row, "eqtm_genes"])
  if (length(gene) == 1 && is.na(gene)){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  selTranscripts <- rownames(se)[rowData(se)$Symbol %in% gene]
  if (length(selTranscripts) == 0){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  exp_name <- epi_df[row, ][["id"]]
  
  gexpZ <- gexpZ_list[[epi_df[row, "cell"][[1]]]]
  
  if (!exp_name %in% colnames(gexpZ)){
    return(   list(z = NA, rank = NA, gene = NA))
  }
  zs <- gexpZ[selTranscripts, exp_name]
  z <- zs[which.max(abs(zs))]
  
  rankMat <- apply(gexpZ[selTranscripts, , drop = FALSE], 1, rank)
  rownames(rankMat) <- colnames(gexpZ)
  
  ranks <- rankMat[exp_name, ]
  rank <- ranks[which.max(abs(zs))]
  
  list(z = z, rank = rank, gene = selTranscripts[which.max(abs(zs))])
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
gexp.z.list <- lapply(cells, function(cell) {
  minise <- gexp[, gexp$cell_type == cell]
  mat <- t(apply(assay(minise), 1, function(x) (x - mean(x))/sd(x)))
  colnames(mat) <- minise$`samplenumber:ch1`
  mat
})
names(gexp.z.list) <- cells

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
                   epi_df = res.gse87650.cell.loo.filt, gexpZ_list = gexp.z.list)
res.gse87650.cell.loo.filt$tss_z <- sapply(tss_list, function(x) x$z)
res.gse87650.cell.loo.filt$tss_rank <- sapply(tss_list, function(x) x$rank)
res.gse87650.cell.loo.filt$tss_gene <- sapply(tss_list, function(x) x$gene)

### Get eQTM ####
eqtm_genes <- mclapply(res.gse87650.cell.loo.filt$cpg_ids, function(x) {
  
  cpgs <- strsplit(x, ",")[[1]]
  tab <- subset(eqtm, CpG %in% cpgs)
  if (nrow(tab) == 0){
    return(NA)
  } else {
    unique(tab$CpG_gene)  
  }
}, mc.cores = 5)

res.gse87650.cell.loo.filt$eqtm_genes <- eqtm_genes
eqtm_list <- lapply(seq_len(nrow(res.gse87650.cell.loo.filt)), getGeneseQTM, 
                    epi_df = res.gse87650.cell.loo.filt, gexpZ_list = gexp.z.list)
res.gse87650.cell.loo.filt$eqtm_z <- sapply(eqtm_list, function(x) x$z)
res.gse87650.cell.loo.filt$eqtm_rank <- sapply(eqtm_list, function(x) x$rank)
res.gse87650.cell.loo.filt$eqtm_gene <- sapply(eqtm_list, function(x) x$gene)

### Close gene ####
geneInfo <- lapply(seq_len(nrow(res.gse87650.cell.loo.filt)), getGenesZ, 
                   epi_df = res.gse87650.cell.loo.filt,
                   gexpZ_list = gexp.z.list, gexpRange = gexp_ranges, window = 250e3)
geneInfoNear <- Reduce(rbind, lapply(geneInfo, selectNearest))
colnames(geneInfoNear) <- paste0(colnames(geneInfoNear), "Near" )
res.gse87650.cell.loo.filt$near_z <- geneInfoNear$Zscore
res.gse87650.cell.loo.filt$near_rank <- geneInfoNear$rank
res.gse87650.cell.loo.filt$near_gene <- geneInfoNear$geneName

### Make plots ####
res.gse87650.cell.loo.gexp.sum  <- res.gse87650.cell.loo.filt %>% 
  ungroup() %>%
  select(method, cell, magnitude, reg_id, ends_with(c("z", "rank")), -sz) %>%
  gather(Measure, value, 5:10) %>%
  mutate(exp_type = sapply(strsplit(Measure, "_"), `[`, 1),
         exp_type = factor(exp_type, levels = c("eqtm", "tss", "near")), 
         measure = sapply(strsplit(Measure, "_"), `[`, 2),
         measure = factor(measure, levels = c("z", "rank")), 
         method = factor(method, levels = c("quantile", "beta", "mlm")))


gse87650.gexp.plots <- lapply(cells, function(x) {
  res.gse87650.cell.loo.gexp.sum %>% 
    filter(cell == x) %>%
  ggplot(aes(x = exp_type, y = value, color = method)) +
  geom_violin() +
  geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5, 
               binwidth = 0.2) +
  theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(x) +
    facet_grid(measure ~ method, scales = "free")
})
plot_grid(plotlist = gse87650.gexp.plots)

# png("figures/INMA4.genexp.png")
# inma4.gexp.plot
# dev.off()

gse87650.gexp.methdiff.plots <-  lapply(cells, function(x) {
  res.gse87650.cell.loo.gexp.sum %>% 
    filter(cell == x) %>%
      filter(measure == "z") %>%
  ggplot(aes(x = abs(magnitude), y = abs(value), color = method)) +
  geom_point() +
  geom_smooth(method = "lm") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(x) +
    theme_bw() +
  facet_grid(exp_type ~ method, scales = "free")
})
plot_grid(plotlist = gse87650.gexp.methdiff.plots)

# png("figures/INMA4.genexp.methdiff.png")
# inma4.gexp.methdiff.plot
# dev.off()

meth.gexp.assocs <- lapply(methods, function(m){
  lapply(cells, function(ce) {
    lapply(c("eqtm", "tss", "near"), function(me){
      tab <- res.gse87650.cell.loo.gexp.sum %>%
        filter(measure == "z" & method ==m & exp_type == me & cell == ce) 
      summary(lm(abs(value) ~ abs(magnitude), tab))  
    })
  })
})

res.gse87650.cell.loo.gexp.sum2 <- res.gse87650.cell.loo.gexp.sum %>% 
  mutate(comb_id = paste(reg_id, method)) %>%
  left_join(res.gse87650.cell.loo.rep %>%
              mutate(comb_id = paste(epi_region_id, id, method)) %>%
              select(comb_id, sigDatasets, sigDatasets2), by = "comb_id") %>%
  mutate(simpDatasets = gsub("-", "", sigDatasets),
       cells = ifelse(simpDatasets == "CD4CD8monocyteswh blood", "4 tissues", 
                      ifelse(simpDatasets %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets, "2-3 tissues")),
       cells = factor(cells, levels = c("CD4", "CD8", "monocytes", "wh blood", "2-3 tissues", "4 tissues")),
       simpDatasets2 = gsub("-", "", sigDatasets2),
       cells2 = ifelse(simpDatasets2 == "CD4CD8monocyteswh blood", "4 tissues", 
                      ifelse(simpDatasets2 %in% c("CD4", "CD8", "monocytes", "wh blood"), simpDatasets2, "2-3 tissues")),
       cells2 = factor(cells2, levels = c("CD4", "CD8", "monocytes", "wh blood", "2-3 tissues", "4 tissues")),
       method = factor(method, levels = c("quantile", "beta", "mlm"))) 
  


gse87650.gexp.plots2 <- lapply(cells, function(x) {
  res.gse87650.cell.loo.gexp.sum2 %>% 
    filter(cell == x & exp_type == "eqtm" & !is.na(cells)) %>%
    ggplot(aes(x = cells, y = value, color = cells)) +
    geom_violin() +
    geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 0.1, stackratio = .5, 
                 binwidth = 0.2) +
    theme_bw() +
    theme( plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(x) +
    facet_grid(measure ~ method, scales = "free")
})
plot_grid(plotlist = gse87650.gexp.plots2)

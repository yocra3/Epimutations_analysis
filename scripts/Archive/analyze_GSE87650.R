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
# library(topGO)
library(pheatmap)
library(cowplot)
library(ExperimentHub)
library(limma)

load("results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Remove replicates  ####
gset_filt <- gset[, gset$description == "sample"]

## Compute residuals of pcs  ####
# beta <- meffil:::impute.matrix(getBeta(gset_filt), margin = 1)
# ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
pcs <- meffil.methylation.pcs(getBeta(gset_filt), probe.range = 40000, full.obj = TRUE)
# m <- getM(gset_filt)
# res <- residuals(lmFit(m, pcs[, seq_len(ndim)]), m)
# beta <- ilogit2(res)
# assay(gset_filt) <- beta
# save(gset_filt, file = "results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")

gset_filt$PC1 <- pcs$x[, 1]
gset_filt$PC2 <- pcs$x[, 2]

top_clust <- gset_filt[, colnames(gset_filt)[pcs$x[, 2] > 0]]
low_clust <- gset_filt[, colnames(gset_filt)[pcs$x[, 2] < 0]]

pcs_top <- meffil.methylation.pcs(getBeta(top_clust), probe.range = 40000)
pcs_low <- meffil.methylation.pcs(getBeta(low_clust), probe.range = 40000)

gse87650.top.cc <- epimutations(top_clust[,! top_clust$disease %in% c("HS", "HL")], top_clust[, top_clust$disease %in% c("HS", "HL")], method = "quantile")
gse87650.low.cc <- epimutations(low_clust[,! low_clust$disease %in% c("HS", "HL")], low_clust[, low_clust$disease %in% c("HS", "HL")], method = "quantile")
gse87650.strat.cc <- rbind(gse87650.top.cc, gse87650.low.cc)
save(gse87650.strat.cc, file = "results/epimutations/GSE87650.epimutations_quantile.cases.stratified.Rdata")

#
# case <- gset_filt[, !gset_filt$disease %in% c("HS", "HL")]
# control <- gset_filt[, gset_filt$disease %in% c("HS", "HL")]
#
# ## Run epimutations ####
# methods <- c("beta", "quantile", "mlm")
# names(methods) <- methods
#
# gse87650.cc <- epimutations(case, control, method = "quantile")
# save(gse87650.cc, file = "results/epimutations/GSE87650.epimutations_quantile.cases.residuals.Rdata")
#
#
# res.gse87650.casecontrol.list <- lapply(methods, epimutations, case_samples = case,
#                             control_panel = control)
# save(res.gse87650.casecontrol.list, file = "results/epimutations/GSE87650.epimutations.cases.residuals.Rdata")
#
# res.gse87650.loo.list <- lapply(methods, epimutations_one_leave_out, methy = gset_filt,
#                                 BPPARAM = MulticoreParam(2))
# save(res.gse87650.loo.list, file = "results/epimutations/GSE87650.epimutations.loo.residuals.Rdata")

# Process epimutations ####
plotDisease <- function(set, range){

  miniset <- subsetByOverlaps(set, range)

  df <- getBeta(miniset)

  df <- t(df) %>% data.frame()
  df$id <- colnames(miniset)

  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$disease <- miniset$disease
  df.gath$disease <- ifelse(df.gath$disease %in% c("HL", "HS"), "Healthy", df.gath$disease)
  df.gath$disease <- factor(df.gath$disease, levels = c("Healthy", "CD", "UC"))

  df.gath$Coordinates <- start(rowRanges(miniset)[df.gath$cpg])


  ggplot(df.gath, aes(x = Coordinates, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation level", limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = "Status",
                       values = c("grey", "green", "blue"))

}

## Case control ####
### Preprocess data ####
res.gse87650.cc.df <-  gse87650.strat.cc %>%
    left_join(colData(gset_filt) %>%
                data.frame() %>%
                dplyr::select(Sample_Name, disease, age, Sex, smoking) %>%
                mutate(sample = Sample_Name))

res.gse87650.cc.sum <- res.gse87650.cc.df %>%
    group_by(sample, disease, age, Sex, smoking) %>%
    summarize(n = sum(start != 0 & cpg_n > 2)) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1",
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))

res.gse87650.cc.sum %>%
  group_by(disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

extremeids.cc <- subset(res.gse87650.cc.sum, n > 20)$sample
extremeids.cc.cd <- subset(res.gse87650.cc.sum, n > 45 & disease == "CD")$sample
extremeids.cc.uc <- subset(res.gse87650.cc.sum, n > 20 & disease == "UC")$sample


### Get recurrent epimutations ####
recur.disease.epi <- res.gse87650.cc.df %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(disease, epi_region_id, disease_n) %>%
  dplyr::summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
arrange(recur.disease.epi, desc(freq))

gse87650.recur.epi <- filter(res.gse87650.cc.df, epi_region_id %in% subset(recur.disease.epi, n >= 3)$epi_region_id) %>%
  group_by(epi_region_id) %>%
  summarize(samples = paste(sample, collapse = ";"),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            cpgs = paste(unique(unlist(strsplit(cpg_ids, ","))), collapse = ";"),
            cpg_n = length(unique(unlist(strsplit(cpg_ids, ","))))) %>%
  mutate(size = end - start) %>%
  dplyr::select(epi_region_id, samples, chromosome, start, end, size, cpgs, cpg_n)

write.table(gse87650.recur.epi, file = "figures/GSE87650.recurrEpi.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

candRegsGR <- epimutacions:::get_candRegsGR()

subset(res.gse87650.cc.df, epi_region_id == "chr6_32793355" & cpg_n > 2) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>% data.frame() %>%
  makeGRangesFromDataFrame() %>% GenomicRanges::reduce()

## Crohn recurrent: cerca NR2F2
### Full
plotDisease(gset_filt, GRanges("chr15:96884949-96890880"))

### Regions
plotDisease(gset_filt, GRanges("chr15:96884949-96888024"))
plotDisease(gset_filt, GRanges("chr15:96890452-96890880"))

## Crohn recurrent: cerca DLX5
## Full
plotDisease(gset_filt, GRanges("chr7:96650096-96655889"))

## Regions
plotDisease(gset_filt, GRanges("chr7:96650096-96652115"))
plotDisease(gset_filt, GRanges("chr7:96654782-96655889"))

## Crohn recurrent: cerca TFAP2A
## Full
plotDisease(gset_filt, GRanges("chr6:10381196-10385903"))

## Regions
plotDisease(gset_filt, GRanges("chr6:10381196-10383147"))
plotDisease(gset_filt, GRanges("chr6:10385320-10385903"))


## Crohn recurrent: chr6_32793355 - in HLA region


## Crohn recurrent: cerca OTX1
plotDisease(gset_filt, GRanges("chr2:63279495-63286355"))

## Crohn recurrent: cerca IRX1
plotDisease(gset_filt, GRanges("chr5:3592464-3592638"))
plotDisease(gset_filt, GRanges("chr5:3599012-3602413"))


plotDisease(gset_filt, candRegsGR["chr1_50879560"])

## Crohn recurrent: cerca DMRTA2
plotDisease(gset_filt, GRanges("chr1:50882910-50885352"))
plotDisease(gset_filt, GRanges("chr1:50886393-50886782"))



recur.ibd.epi <- res.gse87650.cc.df %>%
  filter(chromosome != 0 & cpg_n > 2 & !sample %in% extremeids.cc) %>%
  mutate(n_total = length(unique(sample))) %>%
  group_by(epi_region_id, n_total) %>%
  summarize(n = n(),
            freq = n()/n_total) %>%
  ungroup() %>%
  distinct()

arrange(recur.ibd.epi, desc(freq))
## UC has very few recurrent epimutations and we cannot find common epimutations between both

#### Explore individuals with recurrent epimutations ####
recu.regs.cd <- subset(recur.disease.epi, n >= 3 & disease == "CD")$epi_region_id
recu.inds.cd <- lapply(recu.regs.cd, function(reg){
  subset(res.gse87650.cc.df, epi_region_id == reg)$sample
})
table(sapply(recu.inds.cd, function(x) mean(x %in% extremeids.cc.cd)))
mean(sapply(recu.inds.cd, function(x) mean(x %in% extremeids.cc.cd)) == 1)

gset_filt$disease2 <- ifelse(gset_filt$disease %in% c("HL", "HS"), "Control", gset_filt$disease)
gset_filt$epi_disease <- ifelse(gset_filt$Sample_Name %in% extremeids.cc.cd, "CD_Epi",
                                ifelse(gset_filt$Sample_Name %in% extremeids.cc.uc, "UC_Epi", gset_filt$disease2))
gset_filt$epi_disease <- factor(gset_filt$epi_disease , levels = c("Control", "CD", "CD_Epi", "UC", "UC_Epi"))
#### PC whole methylome
# pcs_resid <- meffil.methylation.pcs(getBeta(gset_filt), full.obj = TRUE,
#                                     probe.range = 40000)
pcs.vars <- pcs$sdev^2/sum(pcs$sdev^2)

pcs.resid.plot <- data.frame(pcs$x, disease = gset_filt$epi_disease) %>%
  mutate(Disease = factor(disease, levels = c("Control", "CD", "Epi_CD", "UC", "Epi_UC"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Disease)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", round(pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(name = "Disease",
                     values = c("grey", "green", "darkgreen", "cyan", "blue")) +
  ggtitle("Whole Methylation")

#### PC recurr epimutations
rec.cpgs <- subset(res.gse87650.cc.df, epi_region_id %in% recu.regs.cd)$cpg_ids
rec.cpgs <- unique(unlist(strsplit(rec.cpgs, ",")))
save(rec.cpgs, file =  "results/epimutations/GSE87650.recurrent.cpgs.Rdata")

pcs_rec <- meffil.methylation.pcs(getBeta(gset_filt[rec.cpgs, ]), full.obj = TRUE,
                                  probe.range = 40000)
rec.pcs.vars <- pcs_rec$sdev^2/sum(pcs_rec$sdev^2)

pcs.rec.plot <- data.frame(pcs_rec$x, disease = gset_filt$epi_disease) %>%
  mutate(Disease = factor(disease, levels = c("Control", "CD", "CD_Epi", "UC", "UC_Epi"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = Disease)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1 (", round(rec.pcs.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(rec.pcs.vars[2]*100, 1), "%)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Disease",
                     values = c("grey", "green", "darkgreen", "cyan", "blue")) +
  ggtitle("CpGs in recurrent epimutations")

png("figures/GSE87650.PCs.recurrentEpis.png", width = 700, height = 350)
plot_grid(pcs.resid.plot, pcs.rec.plot, labels = LETTERS[1:2], rel_widths = c(3, 4))
dev.off()


col_colors <- list(
  epi_disease = c("Control" = "black", "CD" = "green", "UC" = "blue",
              "Epi_CD" = "pink", "Epi_UC" = "brown"),
  Sex = c("F" = "purple", "M" = "lightblue")
)
pheatmap(getBeta(gset_filt[rec.cpgs, ]), scale = "none",
         annotation_col  = data.frame(colData(gset_filt)[, c("epi_disease", "Sex"), drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)


## Explore association with cell types
colData(gset_filt)[, c("Bcell", "CD4T", "CD8T",  "Mono", "Neu", "NK", "epi_disease")] %>%
  as.data.frame() %>%
  gather(CellType, Proportion, 1:6) %>%
  ggplot(aes(x = CellType, y = Proportion, color = epi_disease)) +
  geom_boxplot() +
  theme_bw()

## Explore differences in methylation
gset_cdepi <- gset_filt[, gset_filt$epi_disease %in% c("Control", "CD_Epi")]
gset_cdepi$epi_disease <- droplevels(gset_cdepi$epi_disease )
model_cdepi <- model.matrix(~ epi_disease + Sex + age  + PC1 + PC2, colData(gset_cdepi) )
lm_cdepi <- lmFit(getBeta(gset_cdepi), model_cdepi) %>% eBayes()
res_cdepi <- topTable(lm_cdepi, coef = 2, n= Inf)

cpgs_cdepi <- rownames(subset(res_cdepi, adj.P.Val < 0.05))

gset_cd <- gset_filt[, gset_filt$epi_disease %in% c("Control", "CD")]
gset_cd$epi_disease <- droplevels(gset_cd$epi_disease )
model_cd <- model.matrix(~ epi_disease + Sex + age  + PC1 + PC2, colData(gset_cd) )
lm_cd <- lmFit(getBeta(gset_cd), model_cd) %>% eBayes()
res_cd <- topTable(lm_cd, coef = 2, n= Inf)
cpgs_cd <- rownames(subset(res_cd, adj.P.Val < 0.05))

com_cpgs <- intersect(cpgs_cdepi, cpgs_cd)
plot(res_cd[com_cpgs, "logFC"], res_cdepi[com_cpgs, "logFC"])
cor(res_cd[com_cpgs, "logFC"], res_cdepi[com_cpgs, "logFC"])

mean(rec.cpgs %in% cpgs_cdepi)
mean(rec.cpgs %in% cpgs_cd)
mean(rec.cpgs %in% com_cpgs)

library(FlowSorted.Blood.450k)
ref_gset <- preprocessFunnorm(FlowSorted.Blood.450k) %>% mapToGenome() %>% dropMethylationLoci()
ref_gset <- dropLociWithSnps(ref_gset[!seqnames(ref_gset) %in% c("chrX", "chrY"), ])
ref_gset_cells <- ref_gset[, !ref_gset$CellType %in% c("PBMC")]


## Compute epimutations
ref_epi <- epimutations(ref_gset_cells[, ref_gset_cells$CellType != "WBC"],ref_gset_cells[, ref_gset_cells$CellType == "WBC"], method = "quantile")
ref_epi.df <-  ref_epi %>%
    left_join(colData(ref_gset_cells) %>%
                data.frame() %>%
                dplyr::select(Sample_Name, CellType, CellTypeLong ) %>%
                mutate(sample = Sample_Name))

#
recur.ref_cell.epi <- ref_epi.df  %>%
  group_by(CellType) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(CellType, epi_region_id) %>%
  dplyr::summarize(n = length(unique(sample))) %>%
  ungroup() %>%
  spread(CellType, n, fill = 0)




rec.regions.cd <- subset(recur.disease.epi, disease == "CD" & n >= 3)$epi_region_id

subset(recur.ref_cell.epi, epi_region_id %in% rec.regions.cd & pmax(Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK) > 4)


ref_gset_cells$CellType <- relevel(factor(ref_gset_cells$CellType), "WBC")

pc_cells <- meffil.methylation.pcs(getBeta(ref_gset_cells), full.obj = TRUE,
                                  probe.range = 40000)
pc_cells$x %>%
  data.frame() %>%
  mutate(Cell = ref_gset_cells$CellTypeLong) %>%
  ggplot(aes(x = PC1, y = PC2, color = Cell)) +
  geom_point() +
  theme_bw()

#
pc_cells_rec <- meffil.methylation.pcs(getBeta(ref_gset_cells[intersect(rownames(ref_gset_cells), rec.cpgs), ]), full.obj = TRUE,
                                  probe.range = 40000)
pc_cells_rec$x %>%
  data.frame() %>%
  mutate(Cell = ref_gset_cells$CellTypeLong) %>%
  ggplot(aes(x = PC1, y = PC2, color = Cell)) +
  geom_point() +
  theme_bw()




model_ref <- model.matrix(~ CellType, colData(ref_gset_cells) )
lm_ref <- lmFit(getBeta(ref_gset_cells), model_ref) %>% eBayes()
res_ref <- topTable(lm_ref, coef = 2:9, n= Inf)


res_cdepi_sig <- topTable(lm_cdepi, coef = 2, n= Inf) %>%
  mutate(CpG = rownames(.)) %>%
  as_tibble() %>%
  mutate(padj = p.adjust(P.Value, "bonferroni")) %>%
  subset(padj < 0.05)

res_ref_rec <- res_ref %>%
  mutate(padj = p.adjust(P.Value, "bonferroni")) %>%
  subset(padj < 0.05) %>%
  mutate(CpG = rownames(.)) %>%
  as_tibble()

res_ref_cd <- inner_join(res_ref_rec, res_cdepi_sig, by = "CpG") %>%
  select("CpG", "logFC", starts_with("Cell"))
cor(res_ref_cd$logFC, res_ref_cd[,3:10])

makeTab <- function(coef){
  res_ref <- topTable(lm_ref, coef = coef, n= Inf)

  res_ref_rec <- res_ref %>%
    mutate(padj = p.adjust(P.Value, "bonferroni")) %>%
    subset(padj < 0.05) %>%
    mutate(CpG = rownames(.)) %>%
    as_tibble()

  res_ref_rec
}
cellTabs <- lapply(2:9, makeTab)
names(cellTabs) <- levels(ref_gset_cells$CellType)[-1]

pcs_neu <- meffil.methylation.pcs(getBeta(gset_filt[intersect(rownames(gset_filt), cellTabs$Neu$CpG), ]), full.obj = TRUE,
                                  probe.range = 40000)

cell_cpgs_plot <- lapply(names(cellTabs), function(cell){
  pcs_neu <- meffil.methylation.pcs(getBeta(gset_filt[intersect(rownames(gset_filt), cellTabs[[cell]]$CpG), ]), full.obj = TRUE,
                                    probe.range = 40000)
  neu.pcs.vars <- pcs_neu$sdev^2/sum(pcs_neu$sdev^2)

  pcs.neu.plot <- data.frame(pcs_neu$x, disease = gset_filt$epi_disease) %>%
    mutate(Disease = factor(disease, levels = c("Control", "CD", "CD_Epi", "UC", "UC_Epi"))) %>%
    ggplot(aes(x = PC1, y = PC2, color = Disease)) +
    geom_point() +
    theme_bw() +
    scale_x_continuous(name = paste0("PC1 (", round(neu.pcs.vars[1]*100, 1), "%)")) +
    scale_y_continuous(name = paste0("PC2 (", round(neu.pcs.vars[2]*100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = "Disease",
                       values = c("grey", "green", "darkgreen", "cyan", "blue")) +
    ggtitle(paste0("Top CpGs in ", cell))
})
plot_grid(plotlist = cell_cpgs_plot, ncol = 4, nrow = 2)


cell_cors <- sapply(cellTabs, function(x) {
  comb <- inner_join(x, res_cdepi_sig, by = "CpG")
  cor(comb$logFC.x, comb$logFC.y)
})

neu_comb <- inner_join(cellTabs$Neu, res_cdepi_sig, by = "CpG", suffix = c(".Neu", ".CD"))

ggplot(neu_comb, aes(x = logFC.CD, y = logFC.Neu)) +
  geom_point() +
  theme_bw()

neu_cpgs <- subset(neu_comb , sign(logFC.Neu ) == sign(logFC.CD))$CpG
disc_cpgs <- subset(neu_comb , sign(logFC.Neu ) != sign(logFC.CD))$CpG

mean(disc_cpgs %in% rec.cpgs)
mean(neu_cpgs %in% rec.cpgs)

mean(disc_cpgs %in% cpgs_cd)
mean(disc_cpgs %in% com_cpgs)
mean(neu_cpgs %in% cpgs_cd)
mean(neu_cpgs %in% com_cpgs)


eh <- ExperimentHub()
candRegsGR <- eh[["EH6692"]]

neu_GR <- sort(rowRanges(gset_filt[neu_cpgs,]))
neu_epi <- candRegsGR[to(findOverlaps(neu_GR,candRegsGR ))]
neu_epi <- neu_epi[]



## Region with recurrent epimutation matching pattern of neutrophils
a <- subset(res_ref_cd, CellTypeNeu < -0.1)

plotRegion <- function(set, range, var){

  miniset <- subsetByOverlaps(set, range)

  df <- getBeta(miniset)

  df <- t(df) %>% data.frame()
  df$id <- colnames(miniset)

  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$disease <- var
  df.gath$Coordinates <- start(rowRanges(miniset)[df.gath$cpg])


  ggplot(df.gath, aes(x = Coordinates, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation level", limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

}

gset_cd <- gset_filt[, gset_filt$disease != "UC"]
plot_grid(plotRegion(gset_cd, candRegsGR["chr17_79880258"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr17_79880258"],ref_gset_cells$CellTypeLong), ncol = 2)

subset(recur.ref_cell.epi, epi_region_id %in% rec.regions.cd & pmax(Bcell, CD4T, CD8T, Eos, Gran, Mono, Neu, NK) > 4)

## chr1_25253237 -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr1_25253237"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr1_25253237"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr11_67204799 -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr11_67204799"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr11_67204799"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr11_94274867 -> map with neutrophils
plot_grid(plotRegion(gset_cd, candRegsGR["chr11_94274867"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr11_94274867"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr16_29673465 -> map with CD19+ Bcells
plot_grid(plotRegion(gset_cd, candRegsGR["chr16_29673465"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr16_29673465"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr21_46340776  -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr21_46340776"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr21_46340776"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr22_24822802   -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr22_24822802"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr22_24822802"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr6_31695027    -> región demasiado grande
plot_grid(plotRegion(gset_cd, candRegsGR["chr6_31695027"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr6_31695027"],ref_gset_cells$CellTypeLong), ncol = 2)

## chr6_32793355   -> región demasiado grande
plot_grid(plotRegion(gset_cd, candRegsGR["chr6_32793355"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr6_32793355"],ref_gset_cells$CellTypeLong), ncol = 2)

## Download a second
library(GEOquery)
geo_ref_cells <- getGEO("GSE88824")
geo_ref_cells_gset <- makeGenomicRatioSetFromMatrix(exprs(geo_ref_cells[[1]]), pData = pData(geo_ref_cells[[1]]))

geo_ref_cells_gset <- dropLociWithSnps(geo_ref_cells_gset[!seqnames(geo_ref_cells_gset) %in% c("chrX", "chrY"), ]) %>% dropMethylationLoci()
geo_ref_cells_gset$CellType <- geo_ref_cells_gset$`cell type:ch1`
geo_cells_filt <- geo_ref_cells_gset[, geo_ref_cells_gset$`disease state:ch1` == "Control"]

pc_geo <- meffil.methylation.pcs(getBeta(geo_cells_filt), probe.range = 40000, full.obj = TRUE)
pc_geo$x %>%
  data.frame() %>%
  mutate(Cell = geo_cells_filt$CellType) %>%
  ggplot(aes(x = PC1, y = PC2, color = Cell)) +
  geom_point() +
  theme_bw()




## Compute epimutations
cell_epi <- epimutations(geo_cells_filt[, !geo_cells_filt$CellType %in% c("WBC", "WholeBlood")],geo_cells_filt[, geo_cells_filt$CellType %in% c("WBC", "WholeBlood")], method = "quantile")
cell_epi.df <-  cell_epi %>%
    left_join(colData(geo_cells_filt) %>%
                data.frame() %>%
                dplyr::select(geo_accession , CellType ) %>%
                mutate(sample = geo_accession ))

#
recur.cell_epi.epi <- cell_epi.df  %>%
  group_by(CellType) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(CellType, epi_region_id) %>%
  dplyr::summarize(n = length(unique(sample))) %>%
  ungroup() %>%
  spread(CellType, n, fill = 0)

subset(recur.cell_epi.epi, epi_region_id %in% rec.regions.cd & pmax(CD19B, CD4T, CD8T, Monocyte, Neutrophil, NKcell) > 5)


## chr1_25253237 -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr1_25253237"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr1_25253237"],ref_gset_cells$CellType), ncol = 2)

## chr11_67204799 -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr11_67204799"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr11_67204799"],ref_gset_cells$CellType), ncol = 2)

## chr11_94274867 -> map with neutrophils
plot_grid(plotRegion(gset_cd, candRegsGR["chr11_94274867"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr11_94274867"],ref_gset_cells$CellType), ncol = 2)

## chr21_46340776  -> No equivalence
plot_grid(plotRegion(gset_cd, candRegsGR["chr21_46340776"],gset_cd$epi_disease), plotRegion(ref_gset_cells, candRegsGR["chr21_46340776"],ref_gset_cells$CellType), ncol = 2)


### Explore individual epimutations ####
getMeanDifference <- function(cpglist, samp, set){
  cpgs <- strsplit(cpglist, ",")[[1]]
  betas <- getBeta(set[cpgs, ])
  means <- rowMedians(betas[, colnames(betas) != samp, drop = FALSE], na.rm = TRUE)
  diff  <- betas[, samp] - means
  mean(diff, na.rm = TRUE)
}
res.gse87650.cc.reg.df <- subset(res.gse87650.cc.df, chromosome != 0 & cpg_n >= 3)

magnitudes <- mclapply(seq_len(nrow(res.gse87650.cc.reg.df)),
                                        function(i)
                                          getMeanDifference(res.gse87650.cc.reg.df[i, ]$cpg_ids,
                                                            res.gse87650.cc.reg.df[i, ]$sample,
                                                            gset_filt), mc.cores = 10)
res.gse87650.cc.reg.df$magnitude <- unlist(magnitudes)

annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

isTSS <- lapply(res.gse87650.cc.reg.df$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]

  any(sapply(annot[cpgs, ]$UCSC_RefGene_Group, function(i) grepl("TSS", i)))
})
res.gse87650.cc.reg.df$TSS <- unlist(isTSS)
subset(res.gse87650.cc.reg.df, abs(magnitude) > 0.5 & TSS & !sample %in% recu.inds.cd & !epi_region_id %in% recu.regs.cd) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(sample, epi_region_id) %>% data.frame()



plotDisease(gset_filt, GRanges("chr12:118405665-118406158")) ## ksr2
plotDisease(gset_filt, GRanges("chr7:27160520-27160674")) ## HOXA3
plotDisease(gset_filt, GRanges("chr12:51717674-51718112")) ## BIN2
plotDisease(gset_filt, GRanges("chr11:18433554-18433745")) ## LDHC
plotDisease(gset_filt, GRanges("chr11:66511804-66512979")) ## C11orf80
plotDisease(gset_filt, GRanges("chr1:1851439-1851910")) ## TMEM52
plotDisease(gset_filt, GRanges("chr6:101846767-101846797")) ## GRIK2
plotDisease(gset_filt, GRanges("chr10:105978651-105978702")) ## WDR6
plotDisease(gset_filt, GRanges("chr11:60738971-60739178")) ## CD6


res.gse87650.cc.tss <- res.gse87650.cc.reg.df[res.gse87650.cc.reg.df$TSS, ]
genesTSS <- lapply(res.gse87650.cc.tss$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]

  unique(unlist(strsplit(annot[cpgs, ]$UCSC_RefGene_Name, ";")))
})

### Explore genes ####
library(disgenet2r)
pass <- ""
disgenet_api_key <- get_disgenet_api_key(
  email = "carlos.ruiza@upf.edu",
  password = pass )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

res_df <- Reduce(rbind, lapply(unique(unlist(genesTSS)), function(x) {
  res <- gene2disease(gene = x)
  if (is.character(res)){
    return(NULL)
  }
  extract(res)
}))

disease_genes <- res_df[grep("Crohn|Colitis|Bowel", res_df$disease_name),]$gene_symbol
res.gse87650.cc.tss$TSSgene <- genesTSS
filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & !sample %in% extremeids.cc.cd & !epi_region_id %in% recu.regs.cd ) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(Sample_Name) %>% data.frame()


prior.epi1 <- plot_epimutations(filter(res.gse87650.cc.tss, sample == "GSM2337477" & epi_region_id == "chr1_247578552"),
                  gset_filt) + ggtitle("NLRP3 TSS")  ## Gene NLRP3 PMID: 28523573

filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & cpg_n >= 3 & abs(magnitude) > 0.3) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(TSSgene) %>% data.frame()

prior.epi3 <- plotDisease(gset_filt, GRanges("chr2:182321354-182321489")) +
  ggtitle("ITGA4 TSS") ## ITGA4 - methylation PMID: 25902909

prior.epi2 <- plot_epimutations(filter(res.gse87650.cc.tss, sample == "GSM2337551" & epi_region_id == "chr4_154709441"),
                  gset_filt) + ggtitle("SFRP2 TSS") ## SFRP2 PMID: 18716850

plotDisease(gset_filt, GRanges("chr11:1874017-1874320")) ## LSP1 - dudoso
plotDisease(gset_filt, GRanges("chr16:68676451-68676743")) ## CDH3 - dudoso

png("figures/GSE87650.priorEpimutations.png", width = 484, height = 600)
plot_grid(prior.epi1, prior.epi2, prior.epi3, ncol = 1, labels = LETTERS[1:3])
dev.off()

gse87650.prior.epi <- filter(res.gse87650.cc.tss, TSSgene %in% disease_genes & cpg_n >= 3 & abs(magnitude) > 0.25) %>%
  group_by(epi_region_id) %>%
  summarize(samples = paste(sample, collapse = ";"),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            cpgs = paste(unique(unlist(strsplit(cpg_ids, ","))), collapse = ";"),
            cpg_n = length(unique(unlist(strsplit(cpg_ids, ",")))),
            magnitudes = paste(round(magnitude, 2), collapse = ";"),
            TSSgene = unique(unlist(TSSgene))) %>%
  mutate(size = end - start) %>%
  dplyr::select(epi_region_id, samples, chromosome, start, end, size, cpgs, cpg_n, magnitudes, TSSgene)

write.table(gse87650.prior.epi, file = "figures/GSE87650.prioirizedEpi.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

## Leave-one-out ####
### Preprocess data ####
res.gse87650.loo.df <-  res.gse87650.loo.list$quantile %>%
  left_join(colData(gset_filt) %>%
              data.frame() %>%
              dplyr::select(Sample_Name, disease, age, Sex, smoking) %>%
              mutate(sample = Sample_Name)) %>%
  mutate(disease = ifelse(disease %in% c("HL", "HS"), "Healthy", disease),
         disease = factor(disease, levels = c("Healthy", "CD", "UC")))

res.gse87650.loo.sum <- res.gse87650.loo.df %>%
  group_by(sample, disease, age, Sex, smoking) %>%
  summarize(n = sum(!is.na(cpg_n) & cpg_n >= 3)) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n == 1, "1",
                               ifelse(n < 6, "2-5",
                                      ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


disease.burden <- res.gse87650.loo.sum %>%
  group_by(disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = disease, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

png("figures/GSE87650.diseaseburden.png", width = 500, height = 300)
disease.burden
dev.off()

res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  group_by(smoking, disease, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(smoking, disease, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = smoking, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid( ~ disease) +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Smoking") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")
#
# res.gse87650.loo.sum %>%
#   group_by(Sex, disease, n_cat) %>%
#   summarize(n = n()) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   complete( Sex, disease, n_cat, fill = list(n = 0, p = 0)) %>%
#   ggplot(aes(x = Sex, y = p*100, color = n_cat, fill = n_cat)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   facet_grid( ~ disease) +
#   scale_y_continuous(name = "Proportion of individuals") +
#   scale_x_discrete(name = "Sex") +
#   scale_color_discrete(name = "Epimutations per sample") +
#   scale_fill_discrete(name = "Epimutations per sample")

res.gse87650.loo.sum %>%
  mutate(epimut = ifelse(n_cat == 0, "No Epi", "Epi"),
         epimut = factor(epimut, levels = c("No Epi", "Epi"))) %>%
  ggplot(aes(x = epimut, y = age)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid( ~ disease) +
  scale_x_discrete(name = "Number of epimutations") +
  scale_y_continuous(name = "Age")


### Recurrent epimutations - for comparison with case-control ####
recur.disease.epi.loo <- res.gse87650.loo.df %>%
  group_by(disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(disease, epi_region_id, disease_n) %>%
  summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
arrange(recur.disease.epi.loo, desc(freq))



### Factors influencing having an epimutation ####

### Having epimutation vs not having an epimutation
#### Crude
epi_model_crude_risk <- res.gse87650.loo.sum %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
  glm(out ~ disease, ., family = "binomial")
summary(epi_model_crude_risk)

#### Adjusted
epi_model_adj_risk <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  mutate(out = ifelse(n > 1, 1, n)) %>%
  glm(out ~ disease + smoking + age + Sex, ., family = "binomial")
summary(epi_model_adj_risk)

#### Stratify by sex
epi_model_sex_risk <- lapply(c("F", "M"), function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(Sex == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glm(out ~ disease + smoking + age, ., family = "binomial")
})
lapply(epi_model_sex_risk, summary)

### Other covariables
epi_model_covars <- lapply(disease, function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(disease == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glm(out ~ Sex + smoking + age, ., family = "binomial")
})
lapply(epi_model_covars, summary)



### Total number of epimutations (including 0)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson"))

#### Adjusted
disease.mod2.loo.adj <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod2.loo.adj)


### Total number of epimutations (excluding 0)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = n > 0))

#### Adjusted
disease.mod.loo.adj <- res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(n > 0) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod.loo.adj)

### Other covariates
epi_num_covars <- lapply(disease, function(x) {
  res.gse87650.loo.sum %>%
    filter(!is.na(smoking) & smoking != "Don't know") %>%
    filter(disease == x) %>%
    mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
    glmrob(n ~ smoking + age + Sex, data = ., family = "poisson")
})
lapply(epi_num_covars, summary)



### Total number of epimutations (excluding 0 and extreme)
#### Crude
summary(glmrob(n ~ disease, data = res.gse87650.loo.sum, family = "poisson", subset = n > 0 & n < 21))

#### Adjusted
disease.mod.loo.adj2 <- res.gse87650.loo.sum %>%
  filter(!is.na(smoking) & smoking != "Don't know") %>%
  filter(n > 0 & n < 21) %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Ex", "Current"))) %>%
  glmrob(n ~ disease + smoking + age + Sex, data = ., family = "poisson")
summary(disease.mod.loo.adj2)


### Explore genes associated with epimutations ####
cpg.annot <- getAnnotation(gset_filt)

disease <- as.character(unique(res.gse87650.loo.df$disease))
names(disease) <- disease
#### Map to genes with Illumina annotation ####
genes_list <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & cpg_n > 2)$cpg_ids
  cpgs_list <- strsplit(cpgs, ",")
  gene_list <- lapply(cpgs_list, function(cpg){
    if (length(cpg) == 1 & is.na(cpg)){
      return(NA)
    }
    genes <- cpg.annot[cpg, "UCSC_RefGene_Name"]
    genes_vec <- unlist(strsplit(genes, ";"))
    genes_vec <- unique(genes_vec[!is.na(genes_vec)])
  })
  gene_list
})
genic_tab <- sapply(genes_list, function(x) table(sapply(x, length) > 0))

## Healthy vs CD
chisq.test(genic_tab[, 1:2])

## Healthy vs UC
chisq.test(genic_tab[, 2:3])

genic_tab %>%
  as_tibble() %>%
  mutate(Gene = c("Intergenic", "Gene")) %>%
  gather(Disease, N, 1:3) %>%
  group_by(Disease) %>%
  mutate(prop = N/sum(N),
         Disease = factor(Disease, levels = c("Healthy", "CD", "UC"))) %>%
  ggplot(aes(x = Disease, y = prop, fill = Gene)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_fill_discrete(name = "Epimutation location")

genes_vec <- lapply(genes_list, function(x){
  unique(unlist(x))
})


all_cpgs <- names(subsetByOverlaps(rowRanges(gset_filt), candRegsGR))
all_genes <- cpg.annot[all_cpgs, "UCSC_RefGene_Name"]
all_genes <- unlist(strsplit(all_genes, ";"))
all_genes <- unique(all_genes[!is.na(all_genes)])

getGOdata <- function(gene_vec){
  genesv <- factor(as.numeric(all_genes %in% gene_vec))
  names(genesv) <- all_genes

  Data <- new("topGOdata",
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.org,
              nodeSize = 10,
              ID = "alias", mapping = "org.Hs.eg")
}
godata <- lapply(genes_vec, getGOdata)
gotests <- lapply(godata, runTest, algorithm = "weight01", statistic = "fisher")

finTab <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
                     topNodes = length(score(y)))
}, godata, gotests)

## Repeat tests without epimutations from individuals with extreme number
extremeids <- subset(res.gse87650.loo.sum, n > 20)$sample

genes_list2 <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & !sample %in% extremeids & cpg_n > 2)$cpg_ids
  cpgs_list <- strsplit(cpgs, ",")
  gene_list <- lapply(cpgs_list, function(cpg){
    if (length(cpg) == 1 & is.na(cpg)){
      return(NA)
    }
    genes <- cpg.annot[cpg, "UCSC_RefGene_Name"]
    genes_vec <- unlist(strsplit(genes, ";"))
    genes_vec <- unique(genes_vec[!is.na(genes_vec)])
  })
  gene_list
})
genic_tab2 <- sapply(genes_list2, function(x) table(sapply(x, length) > 0))

## Healthy vs CD
chisq.test(genic_tab2[, 1:2])

## Healthy vs UC
chisq.test(genic_tab2[, 2:3])

genic_tab2 %>%
  as_tibble() %>%
  mutate(Gene = c("Intergenic", "Gene")) %>%
  gather(Disease, N, 1:3) %>%
  group_by(Disease) %>%
  mutate(prop = N/sum(N),
         Disease = factor(Disease, levels = c("Healthy", "CD", "UC"))) %>%
  ggplot(aes(x = Disease, y = prop, fill = Gene)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_fill_discrete(name = "Epimutation location")

genes_vec2 <- lapply(genes_list2, function(x){
  unique(unlist(x))
})

godata2 <- lapply(genes_vec2, getGOdata)
gotests2 <- lapply(godata2, runTest, algorithm = "weight01", statistic = "fisher")

finTab2 <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata2, gotests2)
### Resultados no interesantes

#### Map to genes with eQTMs ####
eqtm <- read_delim("data/eqtm.txt.gz", delim = "\t")

all_tcs <- subset(eqtm, CpG %in% all_cpgs)$TC

getGOdata_eQTM <- function(cpg_vec){

  tcs <- subset(eqtm, CpG %in% cpg_vec)$TC

  genesv <- factor(as.numeric(all_tcs %in% tcs))
  names(genesv) <- all_tcs

  Data <- new("topGOdata",
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.db,
              nodeSize = 10,
              affyLib = "hta20transcriptcluster.db")
}

cpgs_list <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & cpg_n > 2)$cpg_ids
  cpgs_list <- unlist(strsplit(cpgs, ","))
  cpgs_vec <- unique(cpgs_list[!is.na(cpgs_list)])
})

godata_eqtm <- lapply(cpgs_list, getGOdata_eQTM)
gotests_eqtm <- lapply(godata_eqtm, runTest, algorithm = "weight01", statistic = "fisher")

finTab_eqtm <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata_eqtm, gotests_eqtm)


cpgs_list_filt <- lapply(disease, function(x){
  cpgs <- filter(res.gse87650.loo.df, chromosome != 0 & disease == x & !sample %in% extremeids & cpg_n > 2)$cpg_ids
  cpgs_list <- unlist(strsplit(cpgs, ","))
  cpgs_vec <- unique(cpgs_list[!is.na(cpgs_list)])
})

godata_eqtm_f <- lapply(cpgs_list_filt, getGOdata_eQTM)
gotests_eqtm_f <- lapply(godata_eqtm_f, runTest, algorithm = "weight01", statistic = "fisher")

finTab_eqtm_f <- Map(function(x, y) {
  GenTable(x, w0 = y, orderBy = "w0",
           topNodes = length(score(y)))
}, godata_eqtm_f, gotests_eqtm_f)




finTab_eqtm_df <- cbind(finTab_eqtm_f[[2]][, 1:2], Reduce(cbind, lapply(c(2, 1, 3), function(x){
  df <- finTab_eqtm_f[[x]]
  rownames(df) <- df$GO.ID
  df <- df[finTab_eqtm_f[[2]]$GO.ID, 3:6]
  colnames(df) <- paste0(names(finTab_eqtm_f)[x], colnames(df))
  df
}))) %>%
  filter((CDSignificant > 2 & as.numeric(CDw0) < 0.01) |
           (HealthySignificant > 2 & as.numeric(Healthyw0) < 0.01) |
           (UCSignificant > 2 & as.numeric(UCw0) < 0.01))
write.table(finTab_eqtm_df, file = "figures/GSE87650.loo.GOs.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

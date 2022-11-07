#'#################################################################################
#'#################################################################################
#' Combine datasets of CD from GEO
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(meffil)
library(tidyverse)
library(epimutacions)
library(ExperimentHub)
library(cowplot)


## Load data ####
### GSE32148  ####
load("results/preprocess/GSE32148/GSE32148.withNA.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata")
gse32148 <- gset

### GSE81961 ####
load("results/preprocess/GSE81961/GSE81961.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gse81961 <- gset

### GSE87650 ####
load("results/preprocess/GSE87650/GSE87650.wholeblood.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
gse87650 <- gset
gse87650$Disease <- ifelse(gse87650$disease %in% c("UC", "ulcerative colitis"), "UC",
                           ifelse(gse87650$disease %in% c("Crohn", "CD", "Crohn's disease"), "CD", "Control")) 
### GSE112611 ####
load("results/preprocess/GSE112611/GSE112611.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
gse112611 <- gset

# Epimutations
gse32148.cc <- epimutations(gse32148[, gse32148$Disease != "Control"], gse32148[, gse32148$Disease == "Control"], method = "quantile")
save(gse32148.cc, file = "results/epimutations/GSE32148.epimutations.cases.residuals.Rdata")

gse81961.cc <- epimutations(gse81961[, gse81961$Disease != "Control"], gse81961[, gse81961$Disease == "Control"], method = "quantile")
save(gse81961.cc, file = "results/epimutations/GSE81961.epimutations.cases.residuals.Rdata")


load("results/epimutations/GSE87650.epimutations_quantile.cases.stratified.Rdata")
load("results/epimutations/GSE112611.epimutations.casecontrol.residuals.Rdata")

com.epis <- rbind(mutate(gse32148.cc, dataset = "GSE32148"),
                  mutate(gse87650.strat.cc , dataset = "GSE87650"),
                  mutate(gse81961.cc , dataset = "GSE81961"),
                  mutate(gse112611.residuals.cc , dataset = "GSE112611")
)
#
rec.epis.com <- com.epis %>%
  filter(!is.na(epi_region_id)) %>%
  group_by(epi_region_id, dataset) %>%
  dplyr::summarize(n = length(unique(sample))) %>%
  spread(dataset, n,fill = 0) %>%
  arrange(desc(GSE32148 )) %>%
  mutate(N_datasets = (GSE112611 > 0) + (GSE32148 > 0) + (GSE81961 > 0) + (GSE87650 > 0),
         N_epis = GSE112611 + GSE32148 + GSE81961 + GSE87650 )



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
candRegsGR <- epimutacions:::get_candRegsGR()


plotComb <- function(GR){
  plot_grid(plotRegion(gse112611, GR, gse112611$`diagnosis:ch1`) + ggtitle("GSE112611"),
            plotRegion(gse32148, GR, gse32148$Disease) + ggtitle("GSE32148"),
            plotRegion(gse81961, GR, gse81961$Disease) + ggtitle("GSE81961"),
            plotRegion(gse87650, GR, gse87650$Disease)  + ggtitle("GSE87650"))

}
subset(rec.epis.com, N_datasets > 1) %>% arrange(desc(N_datasets), desc(N_epis))

plotComb(candRegsGR["chr2_3641709"]) ## OK
plotComb(GRanges("chr2:3642586-3642847")) ## Demasiada variabilidad - potencial SNP

plotComb("chr6_32060681") ## Muy grande
plotComb(GRanges("chr6:32063394-32064582")) ## Demasiada variabilidad

plotComb(candRegsGR["chr11_2919689"]) ## OK
plotComb(GRanges("chr11:2919763-2925847")) ## Región muy grande

plotComb(candRegsGR["chr4_1577932"])
plotComb(GRanges("chr4:1580377-1582248")) ## Demasiada variabilidad

plotComb(candRegsGR["chr17_6898738"])
plotComb(GRanges("chr17:6898738-6899758")) ## Demasiada variabilidad


plotComb(candRegsGR["chr2_3704363"]) ## Falso positivo

plotComb(candRegsGR["chr6_31275148"]) #
plotComb(GRanges("chr6:31275791-31276667")) ## No claro

plotComb(candRegsGR["chr6_32793355"]) ## Muy grande

plotComb(candRegsGR["chr1_2119531"]) ## no clara
plotComb(GRanges("chr1:2120978-2121521")) ## Mucha variabilidad

plotComb(candRegsGR["chr20_36146069"]) ## no clara
plotComb(GRanges("chr20:36147327-36149452")) ## no clara

plotComb(candRegsGR["chr1_68514311"]) ## no clara
plotComb(GRanges("chr1:68515788-68517454")) ## no clara

plotComb(candRegsGR["chr5_126625777"]) ## 
plotComb(GRanges("chr5:126625742-126626298")) ## Presente en 1 control / MEGF10 gene asociado a Crohn en un GWAS

plotComb(candRegsGR["chr6_33279563"]) ## HLA


plotComb(candRegsGR["chr7_94284258"]) ## 
plotComb(GRanges("chr7:94284486-94286484")) ## Mucha variabilidad

plotComb(candRegsGR["chr16_66612955"]) ## OK
plotComb(GRanges("chr16:66612955-66613334")) ## Epimutación en algunos controls - CMTM1 (imprinting?)


plotComb(candRegsGR["chr2_190044504"]) ## Epimutación en controles

plotComb(candRegsGR["chr11_2720229"]) ## 
plotComb(GRanges("chr11:2721207-2721409")) ## Mucha variabilidad


plotComb(candRegsGR["chr13_100085073"]) ## OK - Intergenic

plotComb(candRegsGR["chr16_88767124"]) ## OK - RNF600

plotComb(candRegsGR["chr17_4688569"]) ## 
plotComb(GRanges("chr17:4688569-4689893")) ## Mucha variabilidad

plotComb(candRegsGR["chr2_63273436"]) ## No solapamiento

plotComb(candRegsGR["chr6_13274151"]) ## OK - PHACTR1. Identified in ASD but with higher difference in these datasets

plotComb(candRegsGR["chr6_30874989"]) ## HLA
plotComb(candRegsGR["chr6_31618987"]) ## HLA
plotComb(candRegsGR["chr6_32128101"]) ## HLA

plotComb(candRegsGR["chr7_27180888"]) ## 
plotComb(GRanges("chr7:27184030-27184441")) ## False positive

plotComb(candRegsGR["chr14_105830248"]) ## Epimutacion en controles

plotComb(candRegsGR["chr13_113584082"]) ## MCFL2 candidato

plotComb(candRegsGR["chr6_25042090"]) ## OK - RIPOR2 

plotComb(candRegsGR["chr12_3861567"]) ## Epimutacion en controles

plotComb(candRegsGR["chr17_76875678"]) ## Epimutacion en controles

plotComb(candRegsGR["chr2_225265371"]) ## No claro

plotComb(candRegsGR["chr6_31686784"]) ## HLA

plotComb(candRegsGR["chr7_27159883"]) ## No overlap

plotComb(candRegsGR["chr8_144891810"]) ## Epimutacion en controles

plotComb(candRegsGR["chr13_105791890"]) ## Epimutacion en controles

plotComb(candRegsGR["chr17_37024020"]) ## 
plotComb(GRanges("chr17:37024020-37024625")) ## Epimutación en controles

plotComb(candRegsGR["chr19_46317840"]) ## 
plotComb(GRanges("chr19:46318514-46319153")) ## RSPH6A/SYMPK

plotComb(candRegsGR["chr19_9731681"]) ## ZNF561

plotComb(candRegsGR["chr6_27635558"]) ## Epimutación en controles

plotComb(candRegsGR["chr6_31650735"]) ## HLA

plotComb(candRegsGR["chr8_101224915"]) ## Epimutación en controles

plotComb(candRegsGR["chr8_125462982"]) ## Epimutación en controles

plotComb(candRegsGR["chr1_161007634"]) ## 
plotComb(GRanges("chr1:161008297-161008977")) ##F11R


plotComb(candRegsGR["chr1_68511835"]) ## 
plotComb(GRanges("chr1:68512539-68513063")) ## Mucha variabilidad

plotComb(candRegsGR["chr10_95461526"]) ## Epimutación en controles

plotComb(candRegsGR["chr11_128553855"]) ## No overlap

plotComb(candRegsGR["chr11_73356316"]) ##
plotComb(GRanges("chr11:73357095-73357396")) ## Mucha variabilidad

plotComb(candRegsGR["chr11_89223981"]) ##
plotComb(GRanges("chr11:89224684-89225014")) ## Mucha variabilidad

plotComb(candRegsGR["chr12_53297384"]) ##
plotComb(GRanges("chr12:53298383-53299310")) ## KRT8

plotComb(candRegsGR["chr15_39871808"]) ##
plotComb(GRanges("chr15:39871808-39872186")) ## Epimutación en controls

plotComb(candRegsGR["chr17_6817409"]) ##Epimutación en controls

plotComb(candRegsGR["chr17_80797881"]) ##
plotComb(GRanges("chr17:80797948-80798706")) ## No se ve claro

plotComb(candRegsGR["chr18_30349111"]) ##
plotComb(GRanges("chr18:30352975-30353699")) ## No se ve claro

plotComb(candRegsGR["chr19_52390810"]) ##
plotComb(GRanges("chr19:52390810-52391605")) ## No se ve claro

plotComb(candRegsGR["chr2_106014950"]) ##
plotComb(GRanges("chr2:106014950-106015767")) ## No se ve claro

plotComb(candRegsGR["chr20_19866743"]) ##
plotComb(GRanges("chr20:19866743-19867136")) ## No se ve claro

plotComb(candRegsGR["chr22_45608345"]) ##
plotComb(GRanges("chr22:45608345-45608686")) ## No se ve claro

plotComb(candRegsGR["chr5_150324958"]) ##
plotComb(GRanges("chr5:150326174-150326497")) ## No se ve claro

plotComb(candRegsGR["chr6_144328421"]) ##
plotComb(GRanges("chr6:144329052-144329829")) ## No se ve claro

plotComb(candRegsGR["chr6_29648161"]) ## HLA
plotComb(candRegsGR["chr6_31544694"]) ## HLA
plotComb(candRegsGR["chr6_32906460"]) ## HLA

plotComb(candRegsGR["chr7_130129946"]) ## 
plotComb(GRanges("chr7:130131676-130132453")) ## No se ve claro

plotComb(candRegsGR["chr7_150019477"]) ## No overlap

plotComb(candRegsGR["chr1_1288586"]) ## Epimutacion en controles

plotComb(candRegsGR["chr11_64115989"]) ## Epimutacion en controles
plotComb(candRegsGR["chr12_7780736"]) ## Epimutacion en controles

plotComb(candRegsGR["chr19_29217858"]) ## Epimutacion en controles
plotComb(candRegsGR["chr19_58877856"]) ## Epimutacion en controles

plotComb(candRegsGR["chr21_37442104"]) ## Epimutacion en controles
plotComb(candRegsGR["chr22_24890330"]) ## 
plotComb(GRanges("chr22:24890690-24890831")) ## No claro

plotComb(candRegsGR["chr3_113160071"]) ## 
plotComb(GRanges("chr3:113160437-113160554")) ## No claro


plotComb(candRegsGR["chr4_118005533"]) ## 
plotComb(GRanges("chr4:118006429-118006812")) ## TRAM1L1

plotComb(candRegsGR["chr5_16178784"]) ## No claro

plotComb(candRegsGR["chr1_92945668"]) ## No claro
plotComb(GRanges("chr1:92945668-92947035")) ## Epimutación en controles

plotComb(candRegsGR["chr1_9555113"]) ## 
plotComb(GRanges("chr1:9555342-9555557")) ## Intergenic

plotComb(candRegsGR["chr11_4628823"]) ## 
plotComb(GRanges("chr11:4629195-4629597")) ## TRIM68

plotComb(candRegsGR["chr11_61061699"]) ## Small

plotComb(candRegsGR["chr13_47370386"]) ## Small

plotComb(candRegsGR["chr14_101289470"]) ## 
plotComb(GRanges("chr14:101290867-101292392")) ## NO claro

plotComb(candRegsGR["chr14_20903410"]) ## small
plotComb(candRegsGR["chr14_76445988"]) ## small

plotComb(candRegsGR["chr15_25068564"]) ## No claro

plotComb(candRegsGR["chr15_65197576"]) ## No claro

plotComb(candRegsGR["chr16_1840421"]) ## 
plotComb(GRanges("chr16:1844927-1845460")) ## Epimutacion en controles

plotComb(candRegsGR["chr16_742426"]) ## 
plotComb(GRanges("chr16:745663-746132")) ## Epimutacion en controles

plotComb(candRegsGR["chr17_20798895"]) ## No claro

plotComb(candRegsGR["chr17_7036834"]) ## small

plotComb(candRegsGR["chr19_2041905"]) ## small

plotComb(candRegsGR["chr19_46525566"]) ## small

plotComb(candRegsGR["chr19_49000743"]) ## 
plotComb(GRanges("chr19:49000743-49002477")) ## Posible (mucha variabilidad) -LMKT3

plotComb(candRegsGR["chr19_50014987"]) ## 
plotComb(GRanges("chr19:50015523-50015780")) ## No claro

plotComb(candRegsGR["chr2_119602212"]) ## small

plotComb(candRegsGR["chr2_177021702"]) ## small
plotComb(candRegsGR["chr22_31317764"]) ## small

plotComb(candRegsGR["chr3_122295656"]) ## small
plotComb(candRegsGR["chr3_147120248"]) ## small

plotComb(candRegsGR["chr3_182816738"]) ## 
plotComb(GRanges("chr3:182816738-182817584")) ## No claro

plotComb(candRegsGR["chr4_188916496"]) ## 
plotComb(GRanges("chr4:188916496-188916875")) ## No claro

plotComb(candRegsGR["chr4_89299173"]) ## No claro

plotComb(candRegsGR["chr5_126205009"]) ## Small

plotComb(candRegsGR["chr5_135364060"]) ## 
plotComb(GRanges("chr5:135364552-135364827")) ## No claro


plotComb(GRanges("chr6:26224013-26224668")) ## OK - H3C6











new_rec <- subset(rec.epis.com, GSE87650 >= 2)$epi_region_id
new_rec <- new_rec[!is.na(new_rec)]

new_rec_cpgs <- subset(com.epis, epi_region_id %in% new_rec)$cpg_ids
new_rec_cpgs <- unique(unlist(strsplit(new_rec_cpgs, ",")))


new_rec_cpgs <- new_rec_cpgs[new_rec_cpgs %in% rownames(geoCD_comb)]


load( "results/epimutations/GSE87650.recurrent.cpgs.Rdata")
gse105798.rec_cpgs <- geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE105798"]
gse105798_pc_rec <- meffil.methylation.pcs(getBeta(gse105798.rec_cpgs), full.obj = TRUE)
pcsdf_gse105798_rec <- data.frame(gse105798_pc_rec$x[, 1:10]) %>%
  mutate(Disease = gse105798.rec_cpgs$`disease state:ch1` )

png("figures/gse105798_reccpgs.pc.png")
ggplot(pcsdf_gse105798_rec, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()

gse105798_pc <- meffil.methylation.pcs(getBeta(geoCD_comb[, geoCD_comb$Dataset == "GSE105798"]), full.obj = TRUE)
pcsdf_gse105798 <- data.frame(gse105798_pc$x[, 1:10]) %>%
  mutate(Disease = gse105798.rec_cpgs$`disease state:ch1` )

png("figures/gse105798_pc.png")
ggplot(pcsdf_gse105798, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()

gse32148.rec_cpgs <- geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE32148"]
gse32148_pc_rec <- meffil.methylation.pcs(getBeta(gse32148.rec_cpgs), full.obj = TRUE)
pcsdf_gse32148_rec <- data.frame(gse32148_pc_rec$x[, 1:10]) %>%
  mutate(Disease = gse32148.rec_cpgs$`disease state:ch1` )

png("figures/gse32148_reccpgs.pc.png")
ggplot(pcsdf_gse32148_rec, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()

gse32148_pc <- meffil.methylation.pcs(getBeta(geoCD_comb[, geoCD_comb$Dataset == "GSE32148"]), full.obj = TRUE)
pcsdf_gse32148 <- data.frame(gse32148_pc$x[, 1:10]) %>%
  mutate(Disease = gse32148.rec_cpgs$`disease state:ch1` )

png("figures/gse32148_pc.png")
ggplot(pcsdf_gse32148, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()


gse81961.rec_cpgs <- geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE81961"]
gse81961_pc_rec <- meffil.methylation.pcs(getBeta(gse81961.rec_cpgs), full.obj = TRUE)
pcsdf_gse81961_rec <- data.frame(gse81961_pc_rec$x[, 1:10]) %>%
  mutate(Disease = gse81961.rec_cpgs$`disease state:ch1` )

png("figures/gse81961_reccpgs.pc.png")
ggplot(pcsdf_gse81961_rec, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()

gse81961_pc <- meffil.methylation.pcs(getBeta(geoCD_comb[, geoCD_comb$Dataset == "GSE81961"]), full.obj = TRUE)
pcsdf_gse81961 <- data.frame(gse81961_pc$x[, 1:10]) %>%
  mutate(Disease = gse81961.rec_cpgs$`disease state:ch1` )

png("figures/gse81961_pc.png")
ggplot(pcsdf_gse81961, aes(x = PC1, y = PC2, col = Disease)) + geom_point()
dev.off()


col_colors <- list(
  `disease state:ch1` = c("Crohn's disease" = "black", "non-IBD control" = "green")
)
png("figures/gse81961_reccpgs.heatmap.png")
pheatmap(getBeta( geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE81961"]), scale = "none",
         annotation_col  = data.frame(colData(geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE81961"])[, "disease state:ch1", drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

png("figures/gse105798_reccpgs.heatmap.png")
pheatmap(getBeta( geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE105798"]), scale = "none",
         annotation_col  = data.frame(colData(geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE105798"])[, "disease state:ch1", drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

png("figures/gse32148_reccpgs.heatmap.png")
pheatmap(getBeta( geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE32148"]), scale = "none",
         annotation_col  = data.frame(colData(geoCD_comb[new_rec_cpgs, geoCD_comb$Dataset == "GSE32148"])[, "disease state:ch1", drop = FALSE]),
         annotation_colors =  col_colors,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

plotRegion <- function(cpgs, set){
  mat <- getBeta(set[cpgs, ]) %>% data.frame() %>%
    mutate(cpg = cpgs,
            position = start(rowRanges(set[cpgs, ]))) %>%
    gather(Sample_Name, methylation, seq_len(ncol(set))) %>%
    left_join(dplyr::select(as.data.frame(colData(set)), Sample_Name, Disease), by = "Sample_Name")

  ggplot(mat, aes(x = position, y = methylation, group = Sample_Name, color = Disease, fill = Disease)) +
    geom_line() +
    geom_point() +
    theme_bw()
}
plotRegion(cpgs2, geoCD_comb[, geoCD_comb$Dataset == "GSE87650"])
plotRegion(cpgs2, geoCD_comb[, geoCD_comb$Dataset == "GSE32148"]) + ggtitle("GSE32148")
plotRegion(cpgs2, geoCD_comb[, geoCD_comb$Dataset == "GSE81961"]) + ggtitle("GSE81961")

## Combine arrays
geoCD_comb <- Reduce(combineArrays, geo_list)

## Check global patterns ####
pcs <- meffil.methylation.pcs(getBeta(geoCD_comb), full.obj = TRUE)

pcsdf <- data.frame(pcs$x[, 1:10]) %>%
  mutate(dataset = geoCD_comb$Dataset)
ggplot(pcsdf, aes(x = PC1, y = PC2, col = dataset)) + geom_point()

geoCD_comb$Disease <- ifelse(is.na( geoCD_comb$`disease state:ch1`),  geoCD_comb$disease,  geoCD_comb$`disease state:ch1`)
geoCD_comb$Disease <- ifelse(geoCD_comb$Disease %in% c("UC", "ulcerative colitis"), "UC",
  ifelse(geoCD_comb$Disease %in% c("Crohn", "CD", "Crohn's disease"), "CD", "Control") )
geoCD_comb$Sex <- ifelse(is.na( geoCD_comb$`gender:ch1`),  geoCD_comb$Sex,  geoCD_comb$`gender:ch1`)
geoCD_comb$Sex[geoCD_comb$Sex == "F"] <- "Female"
geoCD_comb$Sex[geoCD_comb$Sex == "M"] <- "Male"
geoCD_comb$Sample_Name <- colnames(geoCD_comb)

dir.create("results/preprocess/CD_comb/")
save(geoCD_comb, file = "results/preprocess/CD_comb/CD_comb.GenomicRatioSet.Rdata")


cellcounts_list <- lapply(geo_list, function(x) meffil.estimate.cell.counts.from.betas(getBeta(x), cell.type.reference = "blood gse35069 complete",verbose=T))
cellcounts <- Reduce(rbind, cellcounts_list) %>% data.frame()
cellcounts$Sample_Name <- rownames(cellcounts)
pheno <- left_join(cellcounts, data.frame(colData(geoCD_comb))[, c("Sample_Name", "Disease", "Dataset")], by = "Sample_Name")
pheno$Disease2 <- ifelse(pheno$Sample_Name %in% extreme.cd & pheno$Disease == "CD", "CD+", pheno$Disease)
rownames(pheno) <- pheno$Sample_Name


ggplot(pheno, aes(x = Dataset, y = NK, color = Disease2)) +
  geom_boxplot() +
  theme_bw()





reference.object <- meffil:::get.cell.type.reference("blood gse35069 complete")
reference.object$beta[intersect(new_rec_cpgs, rownames(reference.object$beta)), ]

#
col_colors <- list(
  Disease2 = c("CD" = "black", "CD+" = "blue", "Control" = "white", "UC" = "green")
)
pheatmap(getBeta( geoCD_comb[intersect(rownames(geoCD_comb),rownames(reference.object$beta)), geoCD_comb$Dataset == "GSE87650"]), scale = "none",
         annotation_col  = pheno[pheno$Dataset == "GSE87650", c("Disease", "Disease2"), drop = FALSE],
         annotation_colors =  col_colors,

          show_rownames = FALSE, show_colnames = FALSE)

pc_cell <- meffil.methylation.pcs(getBeta(geoCD_comb[intersect(rownames(geoCD_comb),rownames(reference.object$beta)), geoCD_comb$Dataset == "GSE87650"]), full.obj = TRUE)
cot <- data.frame(pc_cell$x[, 1:2]) %>%
  mutate(Sample_Name = geoCD_comb[, geoCD_comb$Dataset == "GSE87650"]$Sample_Name) %>%
  left_join(pheno, by = "Sample_Name")

ggplot(cot, aes(x = PC1, y = PC2, col = Disease2)) + geom_point() + theme_bw()

pc_pac <- meffil.methylation.pcs(getBeta(geoCD_comb[, geoCD_comb$Dataset == "GSE87650"]), full.obj = TRUE)
pac <- data.frame(pc_pac$x[, 1:2]) %>%
  mutate(Sample_Name = geoCD_comb[, geoCD_comb$Dataset == "GSE87650"]$Sample_Name) %>%
  left_join(pheno, by = "Sample_Name")
ggplot(pac, aes(x = PC1, y = PC2, col = Disease2)) + geom_point() + theme_bw()

### Very strong batch effect -> apply ComBat
modcombat <- model.matrix(~ Disease + Sex, data = colData(geoCD_comb))

mat_NA <- is.na(getBeta(geoCD_comb))

assay(geoCD_comb) <- meffil:::impute.matrix(getBeta(geoCD_comb))

## Run comBat
m <- getM(geoCD_comb)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = geoCD_comb$Dataset, mod = modcombat, par.prior=TRUE, prior.plots = FALSE)
beta <- ilogit2(combat_M)
beta[mat_NA] <- NA

assay(geoCD_comb) <- beta
save(geoCD_comb, file = "results/preprocess/CD_comb/CD_comb.normalizedComBat.GenomicRatioSet.Rdata")


## Check global patterns ####
pcs_combat <- meffil.methylation.pcs(getBeta(geoCD_comb), full.obj = TRUE)

pcsdf_combat <- data.frame(pcs_combat$x[, 1:10]) %>%
  mutate(dataset = geoCD_comb$Dataset,
        sex = geoCD_comb$Sex,
        disease = geoCD_comb$Disease)
ggplot(pcsdf_combat, aes(x = PC1, y = PC2, col = dataset)) + geom_point()
ggplot(pcsdf_combat, aes(x = PC1, y = PC2, col = disease)) + geom_point()
ggplot(pcsdf_combat, aes(x = PC1, y = PC2, col = sex)) + geom_point()


## Compute residuals of pcs  ####
beta <- meffil:::impute.matrix(getBeta(geoCD_comb), margin = 1)
ndim <- isva::EstDimRMT(beta, FALSE)$dim + 1
m <- getM(geoCD_comb)
res <- residuals(lmFit(m, pcs$x[, seq_len(ndim)]), m)
beta <- ilogit2(res)
assay(geoCD_comb) <- beta
save(geoCD_comb, file = "results/preprocess/CD_comb/CD_comb.withNA.PCAresiduals.GenomicRatioSet.Rdata")

cdcomb.residuals.cc <- epimutations(geoCD_comb[, geoCD_comb$Disease == "CD"], geoCD_comb[, geoCD_comb$Disease == "Control"], method = "quantile")


pcs_resid <- meffil.methylation.pcs(getBeta(geoCD_comb), full.obj = TRUE)

pcsdf_resid <- data.frame(pcs_resid$x[, 1:10]) %>%
  mutate(dataset = geoCD_comb$Dataset,
        sex = geoCD_comb$Sex,
        disease = geoCD_comb$Disease)

ggplot(pcsdf_resid, aes(x = PC1, y = PC2, col = dataset)) + geom_point()
ggplot(pcsdf_resid, aes(x = PC1, y = PC2, col = disease)) + geom_point()
ggplot(pcsdf_resid, aes(x = PC1, y = PC2, col = sex)) + geom_point()




res.comb.cc.df <-  cdcomb.residuals.cc %>%
    left_join(colData(geoCD_comb) %>%
                data.frame() %>%
                dplyr::select(Sample_Name, Disease, Sex, Dataset) %>%
                mutate(sample = Sample_Name))

#
recur.disease_comb.epi <- res.comb.cc.df %>%
  group_by(Disease) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(Disease, Dataset, epi_region_id, disease_n) %>%
  dplyr::summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
spread(recur.disease_comb.epi, Dataset, n ) %>%
  arrange(desc(GSE87650    ), desc(GSE81961 ))

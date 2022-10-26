#'#################################################################################
#'#################################################################################
#' Produce plots and tables of epimutations in GSE84727
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(epimutacions)
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(UpSetR)
library(cowplot)

# Load data ####
load("results/epimutations/GSE84727.epimutations.allSamples.Rdata")
barbosa_ori <- read_delim("results/barbosa_original/GSE84727.barbosa_result.txt", 
                          delim = "\t", 
                          col_types = cols_only("rownames(gset)" = "c",
                                                "seqnames" = "c",
                                                "start" = "d",
                                                "end" = "d",
                                                "width" = "d",
                                                "strand" = "c",
                                                "Sign_individuals_t0.15_n3_w1k_BOTH" = "c",
                                                "Sign_direction_t0.15_n3_w1k_BOTH" = "c",
                                                "Sign_window_t0.15_n3_w1k_BOTH" = "c"))

# Combine epimutacions and perl results ####
GSE84727.comb <- Reduce(rbind, res.GSE84727.list)
GSE84727.comb$method <- rep(names(res.GSE84727.list), sapply(res.GSE84727.list, nrow))

barbosa_comb$outlier_direction <- barbosa_comb$dir_l
barbosa_comb$sample <- barbosa_comb$samp_l
barbosa_comb$method <- "quantile-perl"

## Add lines for samples without epimutations
no_samps <- setdiff(unique(GSE84727.comb$sample), unique(barbosa_comb$sample))
barbosa_comb <- rbind(barbosa_comb, data.frame(chromosome = 0,
                                               start = 0, end = 0, length = NA, 
                                               sz = NA, samp = NA,
                                               cpg_ids = NA, outlier_score = NA, 
                                               outlier_significance = NA,
                                              adj_pvalue = NA, 
                                              outlier_direction = NA,
                                              samp_l = no_samps,
                                              dir_l = NA,
                                              sample = no_samps,
                                              method = "quantile-perl"))

sel_cols <- c("sample", "chromosome", "start", "end", "cpg_ids", "method")
GSE84727.comb <- rbind(GSE84727.comb[, sel_cols], barbosa_comb[, sel_cols])
GSE84727.comb$method <- factor(GSE84727.comb$method, 
                               levels = c("quantile-perl", "quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd"))

levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="quantile"] <- "quantile-R"
levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="isoforest"] <- "iForest"
levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="mahdistmcd"] <- "mah-dist"

# Plots ####
p_samp_dist <- GSE84727.comb %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n < 6, "1-5",
                               ifelse(n < 10, "6-10", "10+"))),
         n_cat = factor(n_cat, levels = c("0", "1-5", "6-10", "10+"))) %>%
  count(method, n_cat) %>% 
  complete(method, n_cat, fill = list(n = 0)) %>%
  distinct() %>%
  ggplot(aes(fill = n_cat, color = n_cat, y = n, x=method)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Epimutations per sample") +
  scale_color_discrete(name = "Epimutations per sample") +
  theme_bw()

png("figures/GSE84727_epimut_samp_dist_stacked.png", width = 800, height = 400)
p_samp_dist
dev.off()

## Overlap between methods ####
### Add predefined epimutations regions
GSE84727.epi <- subset(GSE84727.comb, chromosome != 0)
rstGR <- GenomicRanges::makeGRangesFromDataFrame(GSE84727.epi)
ensembldb::seqlevelsStyle(rstGR) <- "UCSC" ## Ensure chromosomes have the same format
over <- GenomicRanges::findOverlaps(rstGR, candRegsGR)
GSE84727.epi$epi_region_id <- NA
GSE84727.epi$epi_region_id[S4Vectors::from(over)] <- names(candRegsGR[S4Vectors::to(over)])

### Compare overlap based on predefined regions and sample
GSE84727.epi$reg_id <- paste(GSE84727.epi$sample , GSE84727.epi$epi_region_id )
methods <- levels(GSE84727.epi$method)
names(methods) <- methods
upset.list <- lapply(methods, function(x) subset(GSE84727.epi, method == x)$reg_id)

png("figures/GSE84727_epimut_method_overlaps.png", width = 1000, height = 700)
upset(fromList(upset.list), sets = methods, order.by = "freq",
      mainbar.y.label = "Common epimutations", 
      sets.x.label = "Epimutations per method", 
      text.scale = c(1.7, 1.7, 1.2, 1.2, 2, 1.4))
dev.off()

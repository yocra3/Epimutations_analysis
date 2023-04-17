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
library(ExperimentHub)

# Load data ####
load("results/epimutations/GSE84727.epimutations.allSamples.Rdata")
load("results/epimutations/GSE84727.epimutations.allSamples.ramr.Rdata")
barbosa_ori <- read_delim("results/barbosa_original/GSE84727.barbosa_result_short.txt", 
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

get_regions <- function(flag_df, window_sz = 1000, N = 3, pref = "R") {
  
  if(nrow(flag_df) < N) {
    return(data.frame(chr=NA, pos=NA, region=NA))
  }
  
  # Order the input first by chromosome and then by position
  flag_df <- flag_df[with(flag_df, order(chr, pos)), ]
  
  x <- do.call(rbind, lapply(unique(flag_df$chr), function(chr) {
    # Subset by chromosome
    red_df <- flag_df[flag_df$chr == chr, ]
    if(nrow(red_df) < N) {
      return(data.frame(chr=NA, pos=NA, region=NA))
    }
    
    # Get the position of the previous and next probe for each probe in the
    # data.frame. The first and last position get its own position minus/plus
    # the window size to be sure to include them in the resulting data.frame.
    red_df$pos_next <- c(red_df$pos[seq(2, nrow(red_df))], red_df$pos[nrow(red_df)] + window_sz + 1)
    red_df$pos_prev <- c(red_df$pos[1] - window_sz - 1, red_df$pos[seq(1, nrow(red_df) - 1)])
    
    # We add two columns indicating if a probe is within the window size
    # range with its previous and with its next probe
    red_df$in_prev <- red_df$pos - red_df$pos_prev <= window_sz
    red_df$in_next <- red_df$pos_next - red_df$pos <= window_sz
    
    # We drop all the probes that do not have a previous nor a next
    # probe within the range of interest
    red_df <- red_df[red_df$in_prev | red_df$in_next, ]
    if(nrow(red_df) < N) {
      return(data.frame(chr=NA, pos=NA, region=NA))
    }
    
    # Using the cumsum function and by negating the content of the "in_next"
    # column we can define the regions of CpGs within the range since they
    # will be tagged with the same number
    red_df$cum <- cumsum(!red_df$in_next)
    
    # Correct the base position of the change in the region
    red_df$cum2 <- red_df$cum
    for(ii in seq(2, nrow(red_df))) {
      if(red_df$cum[ii] != red_df$cum[ii - 1] & red_df$in_prev[ii] & !red_df$in_next[ii]) {
        red_df$cum2[ii] <- red_df$cum[ii] - 1
      }
    }
    
    # Computing the frequency of each "number" assign to the region we can 
    # know how may probes are in it. We can use this frequency to filter out
    # those regions with less probes than given by N.
    # We also give to the regions a proper name.
    fr <- data.frame(table(red_df$cum2), stringsAsFactors = FALSE)
    fr <- as.numeric(as.character(fr$Var1[fr$Freq >= N]))
    if(length(fr) > 0) {
      fr <- data.frame(current = fr, new = paste0(pref, "_", chr, "_", seq_len(length(fr))))
      red_df <- red_df[red_df$cum2 %in% fr$current, ]
      rownames(fr) <- paste0("O", fr$current)
      red_df$region <- fr[paste0("O", red_df$cum2), "new"]
      
      
      # Since the first and last probe in a chromosome will have TRUE in
      # prev or next distance we need to be sure to drop them if they
      # are not in the window
      red_df$dist_next <- red_df$pos_next - red_df$pos
      red_df$dist_next[length(red_df$dist_next)] <- 0
      red_df <- red_df[red_df$dist_next <= window_sz, ]
      
      # We drop the columns with the flags used for the outlier and region
      # detection
      red_df <- red_df[ , c("chr", "pos", "region", "cpg",
                            "Sign_individuals_t0.15_n3_w1k_BOTH", "Sign_direction_t0.15_n3_w1k_BOTH")]
      return(red_df)
    } else {
      return(data.frame(chr=NA, pos=NA, region=NA))
    }
  }))
  
  return(x[!is.na(x$chr), ])
}
barbosa_ori$chr <- barbosa_ori$seqnames
barbosa_ori$pos <- barbosa_ori$start
barbosa_ori$cpg <- unlist(barbosa_ori[, 1])

barbosa_reg <- get_regions(barbosa_ori[barbosa_ori$Sign_window_t0.15_n3_w1k_BOTH == 1, ])
barbosa_comb <- do.call(rbind, lapply(unique(barbosa_reg$region), function(reg) {
  x <- barbosa_reg[barbosa_reg$region == reg, ]
  if(nrow(x) > 0) {
    data.frame(
      chromosome = x$chr[1],
      start = min(x$pos),
      end = max(x$pos),
      length = max(x$pos) - min(x$pos),
      sz = nrow(x),
      samp = x$Sign_individuals_t0.15_n3_w1k_BOTH[1],
      cpg_ids = paste(x$`rownames.gset.`, collapse = ",", sep = ""),
      outlier_score = NA,
      outlier_significance = NA,
      adj_pvalue = NA,
      outlier_direction = x$Sign_direction_t0.15_n3_w1k_BOTH[1]
    )
  } else {
    empty
  }
}))
barbosa_comb <- barbosa_comb %>%
  mutate(samp_l = strsplit(samp, ","),
         dir_l = strsplit(outlier_direction, ",")) %>%
  unnest(c(samp_l, dir_l))



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
                        ifelse(n == 1, "1",
                          ifelse(n < 6, "2-5",
                               ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+"))) %>%
  count(method, n_cat) %>% 
  # complete(method, n_cat, fill = list(n = 0)) %>%
  distinct() %>%
  ggplot(aes(fill = n_cat, color = n_cat, y = n, x=method)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Epimutations per sample") +
  scale_color_discrete(name = "Epimutations per sample") +
  theme_bw() +
  ylab("N samples") +
  ggtitle("epimutacions") +
  theme(plot.title = element_text(hjust = 0.5))

png("figures/GSE84727_epimut_samp_dist_stacked.png", width = 800, height = 400)
p_samp_dist
dev.off()

## Overlap between methods ####
candRegsGR <- epimutacions:::get_candRegsGR()

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
upset1 <- upset(fromList(upset.list), sets = methods, order.by = "freq",
      nintersects = 16,
      mainbar.y.label = "", 
      sets.x.label = "N Epimutations") 
dev.off()

# Add ramr results ####
mapToEpimutacionsFormat <- function(ramr_gr){
  
  df <- as.data.frame(ramr_gr)
  df <- df[, c("sample", "seqnames", "start", "end", "revmap")]
  df$seqnames <- as.character(df$seqnames)
  colnames(df) <- c("sample", "chromosome", "start", "end", "cpg_ids")
  out_samps <- setdiff(unique(GSE84727.comb$sample), unique(df$sample))
  rbind(df, data.frame(chromosome = 0,
                       start = 0, end = 0, cpg_ids = NA, 
                       sample = no_samps))
  
}

ramr_df_list <- lapply(res.GSE84727.ramr.list, mapToEpimutacionsFormat)
ramr_df <- Reduce(rbind, ramr_df_list) %>%
  mutate(method = rep(c("ramr-IQR", "ramr-beta", "ramr-wbeta"), sapply(ramr_df_list, nrow)))

GSE84727.comb2 <- rbind(GSE84727.comb[, sel_cols], ramr_df[, sel_cols])
GSE84727.comb2$method <- factor(GSE84727.comb2$method, 
                               levels = c("quantile-perl", "quantile-R", "beta", "manova", "mlm", "iForest", "mah-dist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))

## Plots ####
p_samp_dist_ramr <- GSE84727.comb2 %>%
  filter(method %in% c("quantile-R", "beta", "mlm", "ramr-beta", "ramr-IQR", "ramr-wbeta")) %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n == 1, "1",
                               ifelse(n < 6, "2-5",
                                      ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+"))) %>%
  count(method, n_cat) %>% 
  # complete(n_cat, fill = list(n = 0)) %>%
  distinct() %>%
  ggplot(aes(fill = n_cat, color = n_cat, y = n, x=method)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "Epimutations per sample") +
  scale_color_discrete(name = "Epimutations per sample") +
  theme_bw() +
  ylab("N samples") +
  ggtitle("ramr") +
  theme(plot.title = element_text(hjust = 0.5))

png("figures/GSE84727_epimut_samp_dist_stacked.png", width = 800, height = 400)
p_samp_dist
dev.off()


## Overlap ####
GSE84727.epi2 <- subset(GSE84727.comb2, chromosome != 0)
rstGR2 <- GenomicRanges::makeGRangesFromDataFrame(GSE84727.epi2)
ensembldb::seqlevelsStyle(rstGR2) <- "UCSC" ## Ensure chromosomes have the same format
over <- GenomicRanges::findOverlaps(rstGR2, candRegsGR)
GSE84727.epi2$epi_region_id <- NA
GSE84727.epi2$epi_region_id[S4Vectors::from(over)] <- names(candRegsGR[S4Vectors::to(over)])

GSE84727.epi2$reg_id <- paste(GSE84727.epi2$sample , GSE84727.epi2$epi_region_id )
methods <- levels(GSE84727.epi2$method)
names(methods) <- methods
upset.list2 <- lapply(methods[c("quantile-R", "beta", "mlm", "ramr-beta", "ramr-IQR")], 
                     function(x) subset(GSE84727.epi2, method == x)$reg_id)

png("figures/GSE84727_epimut_method_overlaps.png", width = 1000, height = 700)
upset2 <- upset(fromList(upset.list2), sets = c("quantile-R", "beta", "mlm", "ramr-beta", "ramr-IQR"),
      order.by = "freq", nintersects  = 11,
      mainbar.y.label = "", 
      sets.x.label = "N Epimutations") 
dev.off()

upset1_cw <- cowplot::plot_grid(NULL, upset1$Main_bar, upset1$Sizes, upset1$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,2),
                           rel_widths = c(2,3))

upset2_cw <- cowplot::plot_grid(NULL, upset2$Main_bar, upset2$Sizes, upset2$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,2),
                           rel_widths = c(2,3))

grid <- plot_grid(p_samp_dist, p_samp_dist_ramr, upset1_cw, upset2_cw, 
          ncol = 2, rel_heights = c(2, 3), rel_widths = c(1.15, 1), 
          labels = c("A", "C", "B", "D"))
ggsave("figures/GSE84727_epimut_method_panel.eps", width = 8000, height = 4000, dpi = 600, units = "px")

dev.off()


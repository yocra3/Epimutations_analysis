#'#################################################################################
#'#################################################################################
#' Produce plots and tables of epimutations in GSE84727
#'#################################################################################
#'#################################################################################

# Load libraries and data ####
library(epimutacions)
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(UpSetR)
library(cowplot)

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
# Preprocess Garg regions ####
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

barbosa_r <- res.GSE84727.list$quantile
barbosa_r_filt <- subset(barbosa_r, chromosome != 0)


barbosa_ori_GR <- makeGRangesFromDataFrame(barbosa_comb, keep.extra.columns = TRUE)
barbosa_r_GR <- makeGRangesFromDataFrame(barbosa_r_filt, keep.extra.columns = TRUE)

overlap_ori <- sapply(seq_len(length(barbosa_ori_GR)), function(i){
  gr <- barbosa_ori_GR[i]
  tab_GRm <- barbosa_r_GR[barbosa_r_GR$sample == gr$samp_l]
  length(findOverlaps(gr, tab_GRm)) > 0
})

overlap_R <- sapply(seq_len(length(barbosa_r_GR)), function(i){
  gr <- barbosa_r_GR[i]
  tab_GRm <- barbosa_ori_GR[barbosa_ori_GR$samp_l == gr$sample]
  length(findOverlaps(gr, tab_GRm)) > 0
})

# Merge results ####

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

## Merge
sel_cols <- c("sample", "chromosome", "start", "end", "cpg_ids", "method")
GSE84727.comb <- rbind(GSE84727.comb[, sel_cols], barbosa_comb[, sel_cols])
GSE84727.comb$method <- factor(GSE84727.comb$method, 
                               levels = c("quantile-perl", "quantile", "beta", "manova", "mlm", "isoforest", "mahdistmcd"))

levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="quantile"] <- "quantile-R"
levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="isoforest"] <- "iForest"
levels(GSE84727.comb$method)[levels(GSE84727.comb$method)=="mahdistmcd"] <- "mah-dist"

# Plots ####
## Total epimutations per method ####
GSE84727.comb %>%
  group_by(method) %>%
  summarize(n = sum(chromosome != 0))

## Proportion of samples with epimutations ####
p_samp_prop <- GSE84727.comb %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  group_by(method) %>%
  summarize(prop = mean(n > 0)) %>%
  ggplot(aes(x = method, y = prop*100)) +
  geom_bar(stat = "identity") + theme_bw() +
  scale_y_continuous(name = "Samples with epimutations (%)", limits = c(0, 100))
png("figures/GSE84727_epimut_samp_prop.png")
p_samp_prop
dev.off()



## Distribution epimutations per sample ####
GSE84727.comb %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  glm(n ~ method, ., family = "poisson") %>%
  summary()

GSE84727.comb %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  filter(n > 0) %>%
  group_by(method) %>%
  summarize(median = median(n),
            low = quantile(n, 0.25),
            up = quantile(n, 0.75))
  

GSE84727.comb %>%
  group_by(method, sample) %>%
  summarize(n = sum(chromosome != 0)) %>%
  filter(n > 0) %>%
  glm(n ~ method, ., family = "poisson") %>%
  summary()


# png("figures/GSE84727_epimut_samp_dist.png")
# GSE84727.comb %>%
#   group_by(method, sample) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ggplot(aes(x = method, y = n, fill = method)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()
# 
# 
# p_samp_dist <- GSE84727.comb %>%
#   group_by(method, sample) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ungroup() %>%
#   mutate(n_wind = ifelse(n > quantile(n, 0.98), quantile(n, 0.98), n)) %>%
#   ggplot(aes(x = method, y = n_wind, fill = method)) +
#   geom_boxplot() +
#   theme_bw()
# png("figures/GSE84727_epimut_samp_dist_winsor.png")
# p_samp_dist
# dev.off()
#   
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

# png("figures/GSE84727_epimut_samp_dist_winsor_episamps.png")
# GSE84727.comb %>%
#   group_by(method, sample) %>%
#   summarize(n = sum(chromosome != 0)) %>%
#   ungroup() %>%
#   mutate(n_wind = ifelse(n > quantile(n, 0.98), quantile(n, 0.98), n)) %>%
#   filter(n > 0) %>%
#   ggplot(aes(x = method, y = n_wind, fill = method)) +
#   geom_boxplot() +
#   theme_bw()
# dev.off()

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


# upset <- ggdraw() + draw_image("./figures/GSE84727_epimut_method_overlaps.png")
# 
# png("figures/GSE84727_epimut_panel.png", width = 1000, height = 700)
# plot_grid(p_n, p_samp_prop, p_samp_dist, upset, labels = "AUTO", nrow = 2)
# dev.off()
# 
# sapply(names(upset.list), function(x) 
#   length(intersect(upset.list[[x]], unlist(upset.list[names(upset.list) != x])))/length(unique(upset.list[[x]])))
# 
# sapply(names(upset.list)[4:7], function(x) 
#   length(intersect(upset.list[[x]], unlist(upset.list[names(upset.list) != x])))/length(unique(upset.list[[x]])))
# 
# png("figures/GSE84727_epimut_method_overlaps_top.png", width = 1000, height = 700)
# upset(fromList(upset.list), sets = methods[4:7], order.by = "freq",
#       mainbar.y.label = "Common epimutations", 
#       sets.x.label = "Epimutations per method", 
#       text.scale = c(1.7, 1.7, 1.2, 1.2, 2, 1.4))
# dev.off()

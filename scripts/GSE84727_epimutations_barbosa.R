#'#################################################################################
#'#################################################################################
#' Prepare GSE84727 to run barbosa
#'#################################################################################
#'#################################################################################

library(minfi)
library(BiocParallel)
library(epimutacions)
library(readxl)
load("data/GSE84727/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Create file
annot <- data.frame(rowRanges(gset))
out <- cbind(rownames(gset), annot, getBeta(gset))
write.table(out, file = "results/barbosa_original/GSE84727.input.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## Perl code
perl scripts/barbosa_paper2.pl -f results/barbosa_original/GSE84727.input.txt \
-c 1 -s 2 -e 3 -a 6 -1 0.005 -9 0.995 -t 0.15 -p 3 \
-o results/barbosa_original/GSE84727.barbosa_result.txt 

## Load results
library(tidyverse)
library(GenomicRanges)
library(epimutacions)
barbosa_ori <- read.delim("results/barbosa_original/GSE84727.barbosa_result.txt", 
                          header = TRUE)
barbosa_filt <- barbosa_ori[, -c(7:855)]
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
      red_df <- red_df[ , c("chr", "pos", "region", "rownames.gset.",
                            "Sign_individuals_t0.15_n3_w1k_BOTH", "Sign_direction_t0.15_n3_w1k_BOTH")]
      return(red_df)
    } else {
      return(data.frame(chr=NA, pos=NA, region=NA))
    }
  }))
  
  return(x[!is.na(x$chr), ])
}
barbosa_filt$chr <- barbosa_filt$seqnames
barbosa_filt$pos <- barbosa_filt$start

barbosa_reg <- get_regions(barbosa_filt[barbosa_filt$Sign_window_t0.15_n3_w1k_BOTH == 1, ])
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
  
load("results/epimutations/GSE84727.epimutations.allSamples.barbosaori.Rdata")
barbosa_r_ori <- res.GSE84727.list$barbosa
barbosa_r_ori_filt <- subset(barbosa_r_ori, chromosome != 0)

load("results/epimutations/GSE84727.epimutations.allSamples.Rdata")
barbosa_r <- res.GSE84727.list$barbosa
barbosa_r_filt <- subset(barbosa_r, chromosome != 0)


barbosa_ori_GR <- makeGRangesFromDataFrame(barbosa_comb, keep.extra.columns = TRUE)
barbosa_r_GR <- makeGRangesFromDataFrame(barbosa_r_filt, keep.extra.columns = TRUE)
barbosa_r_1.0_GR <- makeGRangesFromDataFrame(barbosa_r_ori_filt, keep.extra.columns = TRUE)

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


overlap_ori1 <- sapply(seq_len(length(barbosa_ori_GR)), function(i){
  gr <- barbosa_ori_GR[i]
  tab_GRm <- barbosa_r_1.0_GR[barbosa_r_1.0_GR$sample == gr$samp_l]
  length(findOverlaps(gr, tab_GRm)) > 0
})

overlap_R1 <- sapply(seq_len(length(barbosa_r_1.0_GR)), function(i){
  gr <- barbosa_r_1.0_GR[i]
  tab_GRm <- barbosa_ori_GR[barbosa_ori_GR$samp_l == gr$sample]
  length(findOverlaps(gr, tab_GRm)) > 0
})

r_excl <- barbosa_r_filt[!overlap_R, ]
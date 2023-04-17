#'#################################################################################
#'#################################################################################
# Simulate epimutations based on GSE84727 data
#'#################################################################################
#'#################################################################################

# Prepare data ####
## Load libraries ####
library(ramr)
library(epimutacions)
library(minfi)
library(tidyverse)

## Load data ####
load("data/GSE84727.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Convert GRSet to ramr format
data.ranges <- granges(gset)
data.betas  <- getBeta(gset)
sample.ids  <- colnames(gset)
mcols(data.ranges) <- data.betas

## Simulate data ####
set.seed(27)

# unique random AMRs
amrs.unique_020 <-
  simulateAMR(data.ranges, nsamples = 100, regions.per.sample = 10,
              min.cpgs = 3, merge.window = 1000, dbeta = 0.20)

amrs.unique_040 <-
  simulateAMR(data.ranges, nsamples = 100, regions.per.sample = 10,
              min.cpgs = 3, merge.window = 1000, dbeta = 0.40)

# random noise outside of AMR regions
noise_020 <-
  simulateAMR(data.ranges, nsamples = 100, regions.per.sample = 10,
              exclude.ranges = amrs.unique_020,
              min.cpgs = 1, max.cpgs = 1, merge.window = 1, dbeta = 0.5)

noise_040 <-
  simulateAMR(data.ranges, nsamples = 100, regions.per.sample = 10,
              exclude.ranges = amrs.unique_040,
              min.cpgs = 1, max.cpgs = 1, merge.window = 1, dbeta = 0.5)


sim.data020 <-
  simulateData(data.ranges, nsamples = 100,
               amr.ranges = c(amrs.unique_020, noise_020), cores = 2)

save(sim.data020, amrs.unique_020, noise_020,
     file = "results/simulations/sim020_ramr_format.Rdata")
simGset020 <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.data020)))
save(simGset020, file = "results/simulations/sim020_GRSet.Rdata")

sim.data040 <-
  simulateData(data.ranges, nsamples = 100,
               amr.ranges = c(amrs.unique_040, noise_040), cores = 2)

save(sim.data040, amrs.unique_040,  noise_040,
     file = "results/simulations/sim040_ramr_format.Rdata")
simGset040 <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.data040)))
save(simGset040, file = "results/simulations/sim040_GRSet.Rdata")

# Effect size ####
## Run epimutacions methods using LOO ####
epimut_res020 <- lapply(names(epi_parameters()), function(met){
  epimutations_one_leave_out(simGset020, method = met)
})
names(epimut_res020) <- names(epi_parameters())


epimut_res040 <- lapply(names(epi_parameters()), function(met){
  epimutations_one_leave_out(simGset040, method = met)
})
names(epimut_res040) <- names(epi_parameters())
save(epimut_res020, epimut_res040, file = "results/simulations/sim_epimut_res.Rdata")

## Run RAMR ####
ramr_res020 <- lapply(c("IQR", "beta", "wbeta"), function(met){
  getAMR(sim.data020, ramr.method = met, min.cpgs = 3,
         qval.cutoff = 1e-3,
         merge.window = 1000, cores = 2)
})

ramr_res040 <- lapply(c("IQR", "beta", "wbeta"), function(met){
  getAMR(sim.data040, ramr.method = met, min.cpgs = 3,
         qval.cutoff = 1e-3,
         merge.window = 1000, cores = 2)
})
names(ramr_res040) <- names(ramr_res020) <- c("ramr-IQR", "ramr-beta", "ramr-wbeta")
save(ramr_res020, ramr_res040, file = "results/simulations/sim_ramr_res.Rdata")

## Evaluate ####
all.ranges <- getUniverse(sim.data020, min.cpgs = 3, merge.window = 1000)

epimut_res020GR <- lapply(epimut_res020, function(x) {
  GR <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  GR[seqnames(GR) != 0]
})
epimut_res040GR <- lapply(epimut_res040, function(x) {
  GR <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  GR[seqnames(GR) != 0]
})

# true positives
getValues <- function(res, epi_regs, all_ranges){
  
  
  all_neg <- all_ranges[all_ranges %outside% epi_regs]
  
  if(length(res) == 0){
    return(  c(TP = 0,
               FP = 0, 
               TN = length(all_neg),
               FN = length(epi_regs))
            
    ) 
  } else {
  c(TP = sum(sapply(unique(res$sample), function(x) sum(subset(res, sample == x) %over% subset(epi_regs, sample == x)))), 
    FP = sum(res %outside% epi_regs),
    TN = sum(all_neg %outside%  res),
    FN = sum(sapply(unique(res$sample), function(x) sum(subset(epi_regs, sample == x) %outside% subset(res, sample == x))))
  ) 
  }
}

summary_tab20 <- sapply(c(ramr_res020, epimut_res020GR), 
                        getValues, epi_regs = amrs.unique_020,
                        all_ranges = all.ranges)

summary_tab40 <- sapply(c(ramr_res040, epimut_res040GR), 
                        getValues, epi_regs = amrs.unique_040,
                        all_ranges = all.ranges)

rmar_df <- rbind(summary_tab20 %>% 
  t() %>%
  data.frame() %>%
  mutate(Method = rownames(.), 
         Dataset = "0.20"),
  summary_tab40 %>% 
    t() %>%
    data.frame() %>%
    mutate(Method = rownames(.), 
           Dataset = "0.40")
)

plot_main <- rmar_df %>%
  mutate(TPR = TP / (TP + FN), 
         FDR = FP / (TP + FP)
         ) %>%
  gather(Measure, Value, 7:8) %>%
  mutate(Measure = factor(Measure, levels = c("TPR", "FDR")),
         Method = factor(Method, levels = c("quantile", "beta", "manova", "mlm", "iForest", "mahdist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))) %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(Dataset ~ Measure) +
  scale_fill_manual(values = c("#E69F00", "#F0E442","cyan", "deepskyblue2", "deepskyblue4", "blue", "grey10", "grey50", "grey90")) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))

png("figures/simulations_effectsize.png", res = 300, width = 2500, height = 1000)
plot_main 
dev.off()


# Sample size ####
## Simulate data a general dataset with only unique epimutations and subset ####
set.seed(27)
Nsamp <- seq(20, 90, 10)
sim_sample_res <- parallel::mclapply(Nsamp, function(N){
  
    samps <- sample(colnames(mcols(sim.data040)), N)
    epimut_samp <- lapply(names(epi_parameters()), function(met){
      epimutations_one_leave_out(simGset040[, samps], method = met)
    })
    names(epimut_samp) <- names(epi_parameters())
    
    epimut_sampGR <- lapply(epimut_samp, function(x) {
      GR <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
      GR[seqnames(GR) != 0]
    })
    

    ramr_samp <- lapply(c("IQR", "beta", "wbeta"), function(met){
      subsamp <- sim.data040[, samps]
      getAMR(subsamp, ramr.method = met, min.cpgs = 3,
             merge.window = 1000, cores = 2, qval.cutoff = 1e-3)
    })
    names(ramr_samp) <- c("IQR", "rmar-beta", "wbeta")
    
    list(samples = samps, results = c(epimut_sampGR, ramr_samp))
  
}, mc.cores = 8)
save(sim_sample_res, file = "results/simulations/sim_res_samplesize.Rdata")


## Create tables
sim_sample_tabs <- lapply(sim_sample_res, function(iterL){
    regs <- subset(amrs.unique_040, sample %in% iterL$samples)
    names(iterL$res)[7:9] <- c("ramr-IQR", "ramr-beta", "ramr-wbeta")
    vals <- sapply(iterL$res, getValues, epi_regs = regs,
                   all_ranges = all.ranges)
    
    t(vals) %>% 
      data.frame() %>%
      mutate(Method = rownames(.)) %>%
      mutate(TPR = TP / (TP + FN), 
             FDR = FP / (TP + FP)) %>%
      as_tibble()
})
sim_sample_tabs_N <- Reduce(rbind, sim_sample_tabs) %>%
  mutate(N = rep(Nsamp, sapply(sim_sample_tabs, nrow)))

## epimutations

sample_size_plot <- sim_sample_tabs_N %>%
  select(Method, TPR, FDR, N) %>%
  gather(Measure, Value, 2:3) %>%
  mutate(Measure = factor(Measure, levels = c("TPR", "FDR")),
         Method = factor(Method, levels = c("quantile", "beta", "manova", "mlm", "iForest", "mahdist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))) %>%
  ggplot(aes(x = N, y = Value, color = Method)) +
  scale_color_manual(values = c("#E69F00", "#F0E442","cyan", "deepskyblue2", "deepskyblue4", "blue", "grey10", "grey50", "grey80")) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_grid(~ Measure)

png("figures/simulations_samplesize.png", res = 300, width = 2500, height = 1000)
sample_size_plot 
dev.off()




# Epimut. frequency ####
### Subset input data to reigons that can have epimutations
# data.ranges.epi <- subsetByOverlaps(data.ranges, all.ranges)

### Only add non-unique epimutations
simulate_non_unique <- function(N, diff, Nregs){
  
  epi_regs <- sample(all.ranges, 1000)
  data.ranges.epi <- subsetByOverlaps(data.ranges, epi_regs)
  nonunique <-
    simulateAMR(data.ranges.epi, nsamples = N, 
                regions.per.sample = Nregs, samples.per.region = N, min.cpgs = 3,
                merge.window = 1000, dbeta = diff)
   sim.data <-
    simulateData(data.ranges.epi, nsamples = 100,
                 amr.ranges = nonunique, cores = 2)
   gset <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.data)))
   
  list(regs = nonunique, ramr_data = sim.data, epimut_data = gset)
}

set.seed(27)
nonunique_sim <- lapply(2:10, simulate_non_unique, diff = 0.4, Nregs = 100)
save(nonunique_sim, file = "results/simulations/nonunique_data.Rdata")

## Run methods ####
nonunique_sim_ramr <- lapply(nonunique_sim, function(iter){
  lapply(c("IQR", "beta", "wbeta"), function(met){
    getAMR(iter$ramr_data, ramr.method = met, min.cpgs = 3,
           qval.cutoff = 1e-3,
           merge.window = 1000, cores = 2)
  })
})
save(nonunique_sim_ramr, file = "results/simulations/nonunique_ramr.Rdata")

nonunique_sim_cc <- lapply(nonunique_sim, function(iter){
  lapply(names(epi_parameters()), function(met){
    gset <- iter$epimut_data
    samps <-  unique(iter$regs$sample)
    epimutations(gset[, samps], gset[, !colnames(gset) %in% samps], method = met)
  })
})
save(nonunique_sim_cc, file = "results/simulations/nonunique_casecontrol.Rdata")

nonunique_sim_loo <- lapply(nonunique_sim, function(iter){
  lapply(names(epi_parameters()), function(met){
    epimutations_one_leave_out(iter$epimut_data, method = met)
  })
})
save(nonunique_sim_loo, file = "results/simulations/nonunique_loo.Rdata")


processIter <- function(iter, sim){
  
  res <- sapply(iter, function(x) 
    getValues(x, epi_regs =  sim$regs, 
              all_ranges = getUniverse(sim$ramr_data, min.cpgs = 3, merge.window = 1000)
))
  res %>%
    t() %>%
    data.frame() %>%
    mutate(Method = c("ramr-IQR", "ramr-beta", "ramr-wbeta")) %>%
    select(Method, TP, FN, TN, FP)
}

nonunique_sim_ramr_l <- Map(processIter, nonunique_sim_ramr, nonunique_sim)
nonunique_sim_ramr_df <- Reduce(rbind, nonunique_sim_ramr_l) %>%
  mutate(N = rep(2:10, sapply(nonunique_sim_ramr_l, nrow)),
         Mode = "LOO",
         Package = "ramr")


processIterEpimut <- function(iter, sim){
  
  res <- sapply(iter, function(x) {
    regs <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
    regs <- regs[seqnames(regs) != "0"]
    getValues(regs, epi_regs =  sim$regs,
              all_ranges = getUniverse(sim$ramr_data, min.cpgs = 3, merge.window = 1000)
    )
    })
  res %>%
    t() %>%
    data.frame() %>%
    mutate(Method = names(epi_parameters())) %>%
    select(Method, TP, FN, TN, FP)
}

nonunique_sim_cc_l <- Map(processIterEpimut, nonunique_sim_cc, nonunique_sim)
nonunique_sim_cc_df <- Reduce(rbind, nonunique_sim_cc_l) %>%
  mutate(N = rep(2:10, sapply(nonunique_sim_cc_l, nrow)),
         Mode = "Case-Control", 
         Package = "epimutacions")

nonunique_sim_loo_l <- Map(processIterEpimut, nonunique_sim_loo, nonunique_sim)
nonunique_sim_loo_df <- Reduce(rbind, nonunique_sim_loo_l) %>%
  mutate(N = rep(2:10, sapply(nonunique_sim_loo_l, nrow)),
         Mode = "LOO",
         Package = "epimutacions")

nonunique_df <- rbind(nonunique_sim_ramr_df, nonunique_sim_cc_df, nonunique_sim_loo_df) %>%
  mutate(TPR = TP / (TP + FN))

plot_nonunique <- nonunique_df %>%
  mutate(Method = factor(Method, levels = c("quantile", "beta", "manova", "mlm", "iForest", "mahdist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))) %>%
  ggplot(aes(x = N, y = TPR, color = Method)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#E69F00", "#F0E442","cyan", "deepskyblue2", "deepskyblue4", "blue", "grey10", "grey50", "grey90")) +
  theme_bw() +
  xlab("Epimut. Freq (%)") +
  facet_grid(. ~ Package +  Mode  )


png("figures/simulations_freq.png", res = 300, width = 2500, height = 1000)
plot_nonunique 
dev.off()

# Time ####
### Create reference without epimutations
set.seed(27)
sim.data_control <-
  simulateData(sim.data040, nsamples = 100,
               amr.ranges = noise_040, cores = 2)
save(sim.data_control, file = "results/simulations/sim.reference_ramr_format.Rdata")

sim.data_ctrlGRS <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.data_control)))

## Compute time ####
install.packages("microbenchmark")
library(microbenchmark)


ramr_tests <- lapply(c("IQR", "beta", "wbeta"), function(x) {
  bquote( getAMR(sim.data040, ramr.method = .(x), min.cpgs = 3,
                 qval.cutoff = 1e-3,
                 merge.window = 1000, cores = 1))
})
names(ramr_tests) <- paste0("ramr-", c("IQR", "beta", "wbeta"))

cc_tests <- lapply(names(epi_parameters()), function(met){
  bquote(epimutations(simGset040, sim.data_ctrlGRS, method = .(met)))
})
names(cc_tests) <- paste("cc", names(epi_parameters()))

loo_tests <- lapply(names(epi_parameters()), function(met){
  bquote(epimutations_one_leave_out(simGset040, method = .(met)))
})
names(loo_tests) <- paste("loo", names(epi_parameters()))


time_report <- microbenchmark(
  list = c(cc_tests, loo_tests, ramr_tests), times = 5)
save(time_report, file = "results/simulations/sim_res_time.Rdata")

plot_time <- data.frame(time_report) %>%
  mutate(Method = recode(expr, IQR = "ramr-IQR" , beta = "ramr-beta", wbeta = "ramr-wbeta"), 
         Package = ifelse(grepl("ramr", Method), "ramr", "epimutacions"),
         Mode = ifelse(grepl("cc", Method), "Case-Control", "LOO"),
         Method = gsub("cc ", "", Method),
         Method = gsub("loo ", "", Method),
         Method = factor(Method, levels = c("quantile", "beta", "manova", "mlm", "iForest", "mahdist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))) %>%
  ggplot(aes(x = Method, y = time*10^-9, fill = Method)) +
  geom_boxplot() +
  ylab("Time (s)") +
  scale_fill_manual(values = c("#E69F00", "#F0E442","cyan", "deepskyblue2", "deepskyblue4", "blue", "grey10", "grey50", "grey90")) +
  facet_grid(~ Package + Mode, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))


png("figures/simulations_time.png", res = 300, width = 2500, height = 1000)
plot_time 
dev.off()

# Memory profiling ####
install.packages("bench")
library(bench)

mem_ramr <- mark(exprs = ramr_tests, 
                 iterations = 1, check = FALSE, filter_gc = FALSE)
save(mem_ramr, file = "results/simulations/sim_res_mem_ramr.Rdata")

mem_cc <- mark(exprs = cc_tests,
               iterations = 1, check = FALSE, filter_gc = FALSE)
save(mem_cc, file = "results/simulations/sim_res_mem_cc.Rdata")

lapply(names(epi_parameters()), function(met){
  mem_res <- mark(exprs = loo_tests[paste("loo", met)], 
       iterations = 1, check = FALSE, filter_gc = FALSE)
  save(mem_res, file = paste0("results/simulations/sim_res_mem_loo", met, ".Rdata"))
})

## Manually do quantile 
f <- tempfile()
utils::Rprofmem(f, threshold = 1)
eval(loo_tests[["loo quantile"]])
utils::Rprofmem(NULL)
cot <- read_table(f, col_names = FALSE, col_types = cols_only(X1 = col_double()))
save(cot, file= "results/simulations/sim_res_mem_loo_quantile_df.Rdata")


## Make plot
load("results/simulations/sim_res_mem_cc.Rdata")
load("results/simulations/sim_res_mem_ramr.Rdata")

mem_res_loo <- lapply(c("manova", "mlm", "iForest", "mahdist", "beta"), function(met){
  load(file = paste0("results/simulations/sim_res_mem_loo", met, ".Rdata"))
  mem_res
})
load("results/simulations/sim_res_mem_loo_quantile_df.Rdata")

mem_df <- data.frame(names = as.character(names(
  c(mem_ramr$expression, mem_cc$expression, sapply(mem_res_loo, function(x) x$expression), "loo quantile" = "a"))),
                     memory = c(mem_ramr$mem_alloc, mem_cc$mem_alloc, sapply(mem_res_loo, function(x) x$mem_alloc), sum(cot$X1, na.rm = TRUE))/(1024^3))


plot_mem <- mem_df %>%
  mutate(Method = recode(names, IQR = "ramr-IQR" , beta = "ramr-beta", wbeta = "ramr-wbeta"), 
         Package = ifelse(grepl("ramr", Method), "ramr", "epimutacions"),
         Mode = ifelse(grepl("cc", Method), "Case-Control", "LOO"),
         Method = gsub("cc ", "", Method),
         Method = gsub("loo ", "", Method),
         Method = factor(Method, levels = c("quantile", "beta", "manova", "mlm", "iForest", "mahdist", "ramr-IQR", "ramr-beta", "ramr-wbeta"))) %>%
  ggplot(aes(x = Method, y = memory, fill = Method)) +
  geom_bar(stat = "identity") +
  ylab("Max memory (Gb)") +
  scale_fill_manual(values = c("#E69F00", "#F0E442","cyan", "deepskyblue2", "deepskyblue4", "blue", "grey10", "grey50", "grey90")) +
  facet_grid(~ Package + Mode, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))

png("figures/simulations_memory.png", res = 300, width = 2500, height = 1000)
plot_mem 
dev.off()

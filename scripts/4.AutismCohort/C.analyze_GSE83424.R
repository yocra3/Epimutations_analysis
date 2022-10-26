#'#################################################################################
#'#################################################################################
#' Analyze GSE83424 dataset
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)
library(epimutacions)
library(BiocParallel)
library(meffil)
library(tidyverse)
library(robustbase)
library(pheatmap)
library(cowplot)
library(ExperimentHub)

## Load data ####
load("results/preprocess/GSE83424/GSE83424.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.PCAresiduals.GenomicRatioSet.Rdata")
load("results/epimutations/GSE83424.epimutations.allSamples.residuals.Rdata")
load("results/epimutations/GSE83424.epimutations.casecontrol.residuals.Rdata")


## Add region annotation ####
addAnnotation <- function(rst){
  rstGR <- GenomicRanges::makeGRangesFromDataFrame(rst)
 ## Ensure chromosomes have the same format
 seqlevelsStyle(rstGR) <- "UCSC"
 #Get candidate regions
 candRegsGR <- epimutacions:::get_candRegsGR()
 over <- GenomicRanges::findOverlaps(rstGR, candRegsGR)
 #variables (avoid long code)
 ids <- names(candRegsGR[S4Vectors::to(over)])
 cre <- candRegsGR[S4Vectors::to(over)]$CRE
 cre_type <- candRegsGR[S4Vectors::to(over)]$CRE_type

 rst$epi_region_id[S4Vectors::from(over)] <- ids
 rst$CRE[S4Vectors::from(over)] <- cre
 rst$CRE_type[S4Vectors::from(over)] <- cre_type
 rst
}
gse83424.residuals.loo <- addAnnotation(gse83424.residuals.loo)

# Process epimutations ####
plotDisease <- function(set, range){

  miniset <- subsetByOverlaps(set, range)

  df <- getBeta(miniset)

  df <- t(df) %>% data.frame()
  df$id <- colnames(miniset)

  df.gath <- gather(df, cpg, methylation, seq_len(nrow(miniset)))
  df.gath$disease <- factor(ifelse(miniset$status == "Case", "ASD", "Control"))

  df.gath$Coordinates <- start(rowRanges(miniset)[df.gath$cpg])

  ggplot(df.gath, aes(x = Coordinates, y = methylation, group = id, col = disease)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(name = "DNA methylation level", limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(name = "Status",
                       values = c("red", "grey"))
}

## Case control ####
### Preprocess data ####
gse83424.cc.df <- gse83424.residuals.cc %>%
    left_join(colData(gset) %>%
                data.frame() %>%
                dplyr::select(status, title, geo_accession) %>%
                mutate(sample = geo_accession))

gse83424.cc.sum <- gse83424.cc.df %>%
    group_by(sample, status) %>%
    summarize(n = sum(start != 0 & cpg_n > 2)) %>%
    mutate(n_cat = ifelse(n == 0, "0",
                          ifelse(n == 1, "1",
                                 ifelse(n < 6, "2-5",
                                        ifelse(n < 20, "6-20", "20+")))),
           n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))

### Get recurrent epimutations ####
recur.disease.epi <- gse83424.cc.df  %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(epi_region_id, disease_n) %>%
  summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
arrange(recur.disease.epi, desc(freq))

gse83424.recur.epi <- filter(gse83424.cc.df, epi_region_id %in% subset(recur.disease.epi, n >= 3)$epi_region_id) %>%
  group_by(epi_region_id) %>%
  summarize(samples = paste(title, collapse = ";"),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            cpgs = paste(unique(unlist(strsplit(cpg_ids, ","))), collapse = ";"),
            cpg_n = length(unique(unlist(strsplit(cpg_ids, ","))))) %>%
  mutate(size = end - start) %>%
  dplyr::select(epi_region_id, samples, chromosome, start, end, size, cpgs, cpg_n)

write.table(gse83424.recur.epi, file = "figures/GSE83424.recurrEpi.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

## Recurrent epimutations ####
## NUP210L - no autism
plotDisease(gset, GRanges("chr1:154127138-154127537"))

## Intrónico y mucho ruido
plotDisease(gset, GRanges("chr10:134045514-134045609"))

## LINC01167
plotDisease(gset, GRanges("chr10:134775542-134775949"))

## Intrónico BRSK2 y mucho ruido - autism
plotDisease(gset, GRanges("chr11:1463541-1463662"))

## Intrónico NINJ2 y mucho ruido - no autism
plotDisease(gset, GRanges("chr12:739280-740338"))

## Intergénico
plotDisease(gset, GRanges("chr12:7781004-7781431"))

## OBI1 pero dudoso
plotDisease(gset, GRanges("chr13:79233506-79234435"))

## Intrónico TIMP2 - no autism
plotDisease(gset, GRanges("chr17:76875678-76876239"))

## RP1-97D16.9
plotDisease(gset, GRanges("chr6:27730016-27730563"))

## HLA
plotDisease(gset, GRanges("chr6:31275718-31275881"))


#### Explore individuals with recurrent epimutations ####
recu.inds <- lapply(gse83424.recur.epi$epi_region_id, function(reg){
  subset(gse83424.cc.df, epi_region_id == reg)$sample
}) %>%
  unlist( ) %>%
  table()

### Explore individual epimutations ####
gse83424.cc.reg.df <- subset(gse83424.cc.df, chromosome != 0 & cpg_n >= 3)

annot <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

isTSS <- lapply(gse83424.cc.reg.df$cpg_ids, function(x){
  cpgs <- strsplit(x, ",")[[1]]

  any(sapply(annot[cpgs, ]$UCSC_RefGene_Group, function(i) grepl("TSS", i)))
})
gse83424.cc.reg.df$TSS <- unlist(isTSS)
subset(gse83424.cc.reg.df, abs(delta_beta) > 0.3 & TSS) %>%
  dplyr::select(-starts_with("CRE"), -ends_with("pvalue"), -cpg_ids) %>%
  arrange(sample, epi_region_id) %>% data.frame()

subset(gse83424.cc.reg.df, abs(magnitude) > 0.20 & TSS) %>%
  makeGRangesFromDataFrame() %>% as.character()

plotDisease(gset, GRanges("chr5:157079440-157079520")) ## SOX30 - no autism
plotDisease(gset, GRanges("chr11:18433500-18433745")) ## LDHC - no autism
plotDisease(gset, GRanges("chr10:27702774-27703377")) ## PTCHD3 - no autism
plotDisease(gset, GRanges("chr13:79233506-79234435" )) ## OBI1 - recurrente - no autism
plotDisease(gset, GRanges("chr6:151646540-151646601" )) ## AKAP12 - no autism
plotDisease(gset, GRanges("chr5:80597196-80597828"  )) ## ZCCHC9  - reported in Aida paper
plotDisease(gset, GRanges("chr5:178986291-178986728"  )) ## RUFY1 - no autism
plotDisease(gset, GRanges("chr1:202172778-202172912")) ## LGR6 - dudoso - no autism
plotDisease(gset, GRanges("chr17:6899310-6899577")) ## ALOX12 - dudoso - no autism
plotDisease(gset, GRanges("chr16:2770814-2770901" )) ## PRSS27 - no autism
plotDisease(gset, GRanges("chr20:30134929-30135362"  )) ## MCTS2P - no autism

plotDisease(gset, GRanges("chr2:30669597-30669863")) ## LCLAT1 - dudoso - no autism
plotDisease(gset, GRanges("chr20:5485144-5485294")) ## LINC00654 - dudoso - no autism
plotDisease(gset, GRanges("chr5:78985432-78985592"   )) ## CMYA5 - dudoso - no autism


## Epimutation Aida paper
plotDisease(gset, GRanges("chr1:183439365-183439493")) ## Muy corto
plotDisease(gset, GRanges("chr6:13274151-13274354")) ## Replicada




plot_epimutations_group <- function(dmr, epi_samples, methy, genome = "hg19", genes_annot = FALSE,
                              regulation = FALSE, from = NULL, to = NULL){
  ## NULL arguments
  if(is.null(dmr)){
    stop("The argument 'dmr' must be introduced")
  }
  if(is.null(methy)){
    stop("The argument 'beta' must be introduced")
  }
  if(is.null(genome)){
    stop("The argument 'genome' must be introduced")
  }
  ##Unique DMR
  if(nrow(dmr) > 1){
    warning("more than one DMR introduced (nrow > 1)
            only the first element will be used")
    dmr <- dmr[1,]
  }
  ##Genome assembly
  if(genome != "hg38" & genome != "hg19" & genome != "hg18"){
    stop("Argument 'genome' must be 'hg38', 'hg19' or 'hg18'")
  }
  ##Epimutation start('from') and end('to') possitions
   ## * 'from' and 'to' introduced together
  if(is.null(from) & !is.null(to) | !is.null(from) & is.null(to)){
    stop("Arguments 'from' and 'to' must be provided together")
  }
   ## * 'from' is smaller than 'to'
  if(!is.null(from) & !is.null(to)){
    if(from > to){
      stop("The value of argument 'from' must be smaller than 'to'")
    }
  }

  if (!requireNamespace("grDevices"))
    stop("'grDevices' package not avaibale")

  # DMR column names must be always
  # the same (set the common column names)
  dmr  <- epimutacions:::cols_names(dmr, cpg_ids_col = TRUE)  #epi_plot

  # Set 'from' and 'to' arguments value
  if(is.null(from) & is.null(to)){
    from <- dmr$start - 1000
    to <- dmr$end + 1000
  }

  #Generate GenomicRanges object to contain in the same object:
  ## * Genomic ranges of each CpG in the DMR
  ## * Beta values
  gr <- epimutacions:::create_GRanges_class(methy, dmr[,"cpg_ids"]) #epi_plot

  betas_sd_mean2 <- function(gr, controls){

    if(!is(gr, "GRanges")){
      stop("The argument 'gr' must be ")
    }
    df <- as.data.frame(gr)
    colnames(df) <- c("seqnames", "start",
                      "end", "width", "strand",
                      colnames(GenomicRanges::elementMetadata(gr)))
    #Compute:
    # * beta values
    # * Population mean
    # * 1, 1.5, and 2 standard deviations from the mean of the distribution
    betas <- as.data.frame(S4Vectors::values(gr))
    betas <- betas[, controls]
    mean <- rowMeans(betas)
    sd <- apply(betas,1,sd)
    sd_1_lower <- abs(mean - sd)
    sd_1_upper <- mean + sd
    sd_1.5_lower <- mean - 1.5*sd
    sd_1.5_lower <- ifelse(sd_1.5_lower > 0, sd_1.5_lower, 0)
    sd_1.5_upper <- mean +  1.5*sd
    sd_2_lower <- mean - 2*sd
    sd_2_lower <- ifelse(sd_2_lower > 0, sd_2_lower, 0)
    sd_2_upper <- mean +  2*sd

    sd <- cbind(df[,c("seqnames", "start","end","width","strand")],
                sd_1_lower, sd_1_upper,
                sd_1.5_lower,sd_1.5_upper,sd_2_lower,sd_2_upper)
    mean <- cbind(df[,c("seqnames", "start","end","width","strand")], mean)


    #Melt beta values, mean and sd object (necessary for the ggplot)


    if (requireNamespace("reshape2", quietly = TRUE)){
    beta_values <- reshape2::melt(df, id = c("seqnames",
                                             "start",
                                             "end",
                                             "width",
                                             "strand"))
    mean <- reshape2::melt(mean, id = c("seqnames",
                                        "start",
                                        "end",
                                        "width",
                                        "strand",
                                        "mean"))
    sd <- reshape2::melt(sd, id = c("seqnames",
                                    "start",
                                    "end",
                                    "width",
                                    "strand",
                                    "sd_1_lower",
                                    "sd_1_upper",
                                    "sd_1.5_lower",
                                    "sd_1.5_upper",
                                    "sd_2_lower",
                                    "sd_2_upper"))
    }else{
      stop("'reshape2' package not available")
    }

    #Create the output list
    output <- list("beta_values" = beta_values, "mean" = mean, "sd" = sd)
    return(output)
  }



  betas_sd_mean <- betas_sd_mean2(gr, colnames(methy)[methy$status == "Control"]) #epi_plot

  #Generate variables in 'beta_values' data frame containing:
  # * status: case sample name/'control'
  # * color: 'red' for case sample and 'black' for control sample
  # * lines: 'longdash' for controls and
  #          'solid' for case and population mean
  status <- ifelse(betas_sd_mean$beta_values$variable %in% epi_samples,
                   as.character(betas_sd_mean$beta_values$variable),
                   colData(gset)[betas_sd_mean$beta_values$variable, ]$status)
  betas_sd_mean$beta_values$status <- status
  rm(status)
  lines <- ifelse(!betas_sd_mean$beta_values$status %in% epi_samples,
                  "longdash","solid")
  betas_sd_mean$beta_values$lines <- lines
  rm(lines)
  # colors <- c("Control" = "black",
  #             "mean" = "darkblue",
  #             "Case" = "firebrick")
  #
  # extra_cols <- c("hotpink", "blueviolet", "lightpink", "purple", "lavender")
  # colors <- c(colors, extra_cols[seq_len(length(epi_samples))])
  # names(colors) <- c("Control", "mean", "Case", epi_samples)

  #Generate a variable with the CpGs names
  variable <- betas_sd_mean$beta_values$variable
  names <- betas_sd_mean$beta_values[variable == dmr$sample,]
  rm(variable)
  names$id <- names(gr)

  #Plot epimutations

  plot_betas <- ggplot2::ggplot() +
    ggplot2::geom_line(data = betas_sd_mean$beta_values,
                       ggplot2::aes(x = start,
                                    y = value,
                                    group = variable,
                                    color = status),
                       linetype = betas_sd_mean$beta_values$lines) +
    ggplot2::geom_point(data = betas_sd_mean$beta_values,
                        ggplot2::aes(x = start,
                                     y = value,
                                     group = variable,
                                     color = status))
  plot_sd <- plot_betas +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd,
                         ggplot2::aes(x = start,
                                      ymin = sd_2_lower,
                                      ymax = sd_2_upper),
                         fill = "gray39", alpha = 0.4) +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd,
                         ggplot2::aes(x = start,
                                      ymin = sd_1.5_lower,
                                      ymax = sd_1.5_upper),
                         fill = "gray40", alpha = 0.4) +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd,
                         ggplot2::aes(x = start,
                                      ymin = sd_1_lower,
                                      ymax = sd_1_upper),
                         fill = "gray98", alpha = 0.4)

  plot_mean <-  plot_sd +
    ggplot2::geom_line(data = betas_sd_mean$mean,
                       ggplot2::aes(x = start,
                                    y = mean,
                                    color = "mean")) +
    ggplot2::geom_point(data = betas_sd_mean$mean,
                        ggplot2::aes(x = start, y = mean),
                        show.legend = TRUE)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot_cpg_names <- plot_mean +
      ggrepel::geom_text_repel() +
      ggplot2::annotate(geom = "text",
                        x = names$start,
                        y = names$value + 0.05,
                        label = names$id,
                        color = "black")
  } else {
    stop("'ggrepel' package not avaibale")
  }


  plot <- plot_cpg_names +
    ggplot2::lims(y = c(0,1)) +
    # ggplot2::scale_colour_manual(name = "Status", values = colors) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0(dmr$sample,": ",
                            dmr$seqnames, ":",
                            dmr$start,
                            " - ", dmr$end)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = "Coordinates") +
    ggplot2::labs(y = "DNA methylation level")

}
## Selected epimutations
# BRSK2 y mucho ruido - autism db
reg1 <- filter(gse83424.cc.df, epi_region_id == "chr11_1463541")
prior.epi1 <- plot_epimutations_group(reg1, reg1$sample, gset) + ggtitle("BRSK2")  +
                  scale_color_manual(name = "Status", values =  c("black",  "darkblue", "firebrick",
                                                "hotpink", "blueviolet", "lightpink",
                                                "purple", "lavender"),
                                     breaks = c("Control", "mean", "Case", filter(gse83424.cc.df, epi_region_id == "chr11_1463541")$sample),
                                     labels = c("Control", "Control mean", "Case", filter(gse83424.cc.df, epi_region_id == "chr11_1463541")$title)) +
                  xlim(min(reg1$start) - 10, max(reg1$end) + 10)

## Gene NUP20L PMID: 30970224
reg2 <- filter(gse83424.cc.df, epi_region_id == "chr1_154127138")
prior.epi2 <- plot_epimutations_group(reg2[2, ], reg2$sample, gset) + ggtitle("NUP20L TSS")  +
                  scale_color_manual(name = "Status", values =  c("black",  "darkblue", "firebrick",
                                                "hotpink", "blueviolet", "lightpink",
                                                "purple", "lavender"),
                                     breaks = c("Control", "mean", "Case", reg2$sample),
                                     labels = c("Control", "Control mean", "Case", reg2$title)) +
                  xlim(c(min(reg2$start) - 20, max(reg2$end) + 20))

## ZCCHC9  - reported in Aida paper
reg3 <- filter(gse83424.cc.df, epi_region_id == "chr5_80594700")
prior.epi3 <- plot_epimutations_group(reg3, reg3$sample, gset) + ggtitle("ZCCHC9 TSS")  +
                  scale_color_manual(name = "Status", values =  c("black",  "darkblue", "firebrick",
                                                "hotpink", "blueviolet", "lightpink",
                                                "purple", "lavender"),
                                     breaks = c("Control", "mean", "Case", reg3$sample),
                                     labels = c("Control", "Control mean", "Case", reg3$title)) +
                  xlim(c(min(reg3$start) - 40, max(reg3$end) + 40))

## Replicada Aida
reg4 <- filter(gse83424.cc.df, epi_region_id == "chr6_13274151")
prior.epi4 <- plot_epimutations_group(reg4, reg4$sample, gset) + ggtitle("PHACTR1")  +
                  scale_color_manual(name = "Status", values =  c("black",  "darkblue", "firebrick",
                                                "hotpink", "blueviolet", "lightpink",
                                                "purple", "lavender"),
                                     breaks = c("Control", "mean", "Case", reg4$sample),
                                     labels = c("Control", "Control mean", "Case", reg4$title)) +
                  xlim(c(min(reg4$start) - 20, max(reg4$end) + 20))

## Sup Figure 29
png("figures/GSE83424.priorEpimutations.png", width = 800, height = 600)
plot_grid(prior.epi1, prior.epi2, prior.epi3, prior.epi4, ncol = 2, labels = LETTERS[1:4])
dev.off()

gse83424.cc.out.df <- select(gse83424.cc.reg.df, epi_region_id, title, chromosome, start, end, sz, cpg_n,
  cpg_ids, delta_beta)
colnames(gse83424.cc.out.df)[2] <- "SampleID"
write.table(gse83424.cc.out.df, file = "figures/GSE83424.casecontrol_epis.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

## Leave-one-out ####
### Preprocess data ####
gse83424.loo.df <-  gse83424.residuals.loo  %>%
  left_join(colData(gset) %>%
            data.frame() %>%
            dplyr::select(status, title, geo_accession) %>%
            mutate(sample = geo_accession))

gse83424.loo.sum <- gse83424.loo.df %>%
  group_by(sample, status) %>%
  summarize(n = sum(!is.na(cpg_n) & cpg_n >= 3)) %>%
  mutate(n_cat = ifelse(n == 0, "0",
                        ifelse(n == 1, "1",
                               ifelse(n < 6, "2-5",
                                      ifelse(n < 20, "6-20", "20+")))),
         n_cat = factor(n_cat, levels = c("0", "1", "2-5", "6-20", "20+")))


disease.burden <- gse83424.loo.sum %>%
  group_by(status, n_cat) %>%
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  complete(status, n_cat, fill = list(n = 0, p = 0)) %>%
  ggplot(aes(x = status, y = p*100, color = n_cat, fill = n_cat)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(name = "Proportion of individuals") +
  scale_x_discrete(name = "Disease") +
  scale_color_discrete(name = "Epimutations per sample") +
  scale_fill_discrete(name = "Epimutations per sample")

## Sup Figure 30
png("figures/GSE83424.diseaseburden.png", width = 500, height = 300)
disease.burden
dev.off()

### Recurrent epimutations - for comparison with case-control ####
recur.disease.epi.loo <- gse83424.loo.df %>%
  group_by(status) %>%
  mutate(disease_n = length(unique(sample))) %>%
  filter(chromosome != 0 & cpg_n > 2) %>%
  group_by(status, epi_region_id, disease_n) %>%
  summarize(n = length(unique(sample)),
            freq = length(unique(sample))/disease_n) %>%
  ungroup() %>%
  distinct()
## None of the epimutations from case-control identified with loo


### Factors influencing having an epimutation ####
### Having epimutation vs not having an epimutation
#### Crude
epi_model_crude_risk <- gse83424.loo.sum %>%
    mutate(out = ifelse(n > 1, 1, n)) %>%
  glm(out ~ status, ., family = "binomial")
summary(epi_model_crude_risk)

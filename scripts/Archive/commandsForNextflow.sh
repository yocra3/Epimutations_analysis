#'#################################################################################
#'#################################################################################
#' Run methylation QC data
Rscript 'writeLines(meffil::meffil.snp.names(), con="results/preprocess/methylation/snp-names.txt")'

plink --bfile data/INMA_genos/INMA_HRC.merged.sex.rs --extract results/preprocess/methylation/snp-names.txt --recode A --out results/preprocess/methylation/INMA.methSNPs
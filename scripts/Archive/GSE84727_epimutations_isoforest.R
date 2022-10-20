#'#################################################################################
#'#################################################################################
#' Use results from GSE84727 to refine isoforest results
#'#################################################################################
#'#################################################################################

load("results/epimutations/GSE84727.epimutations.allSamples.Rdata")
thres <- seq(0.5, 1, 0.05)
n_epi <- sapply(thres, function(x) sum(res.GSE84727.list$isoforest$outlier_score > x))

# > sapply(res.GSE84727.list, function(x) sum(x$chr != 0))
# manova        mlm  isoforest mahdistmcd    barbosa       beta 
# 14899      14137      84663       2643        488       2185 
# > data.frame(n_epi, thres)
# n_epi thres
# 1  84663  0.50
# 2  57554  0.55
# 3  37481  0.60
# 4  22423  0.65
# 5  11253  0.70
# 6   4188  0.75
# 7   1022  0.80
# 8    152  0.85
# 9      6  0.90
# 10     0  0.95
# 11     0  1.00

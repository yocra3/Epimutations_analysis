## Normalization methods
methods=(Illumina Noob Quantile Raw SWAN Functional)
for i in "${methods[@]}"
do
    Rscript workflows/bin/createGenomicRatioSets_Raw.R Esteller.minfi${i}Normalization.normalizedRaw.GenomicRatioSet.Rdata Esteller.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds Esteller.minfi${i}Normalization
done

Rscript workflows/bin/createGenomicRatioSets_Raw.R Esteller.BMIQNormalization.normalizedRaw.GenomicRatioSet.Rdata Esteller.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds Esteller.BMIQNormalization

## No ComBat
Rscript bin/createGenomicRatioSets_Raw.R INMA_comb.normalizedRaw.GenomicRatioSet.corrected.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA_comb.normalizedRaw
Rscript bin/createGenomicRatioSets_Raw.R MeDALL_all.normalizedRaw.GenomicRatioSet.Rdata MeDALL_all.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds MeDALL_all.normalizedRaw
Rscript bin/createGenomicRatioSets_Raw.R INMA0combined.normalizedRaw.GenomicRatioSet.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA0combined.normalizedRaw
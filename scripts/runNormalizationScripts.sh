## Normalization methods
methods=(Illumina Noob Quantile Raw SWAN Functional)
for i in "${methods[@]}"
do
    Rscript workflows/bin/createGenomicRatioSets.R Esteller.minfi${i}Normalization.normalizedRaw.GenomicRatioSet.Rdata Esteller.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds Esteller.minfi${i}Normalization
done

for i in "${methods[@]}"
do
    Rscript workflows/bin/INMA0_epimutations_different_normalization.R Esteller.minfi${i}Normalization.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata Esteller.minfi${i}Normalization INMA.commonControlSamples.Rdata
done


## Normalization approaches
approach=(INMA0combined INMA_comb)
for i in "${approach[@]}"
do
    Rscript workflows/bin/INMA0_epimutations_different_batch.R ${i}.normalizedComBat.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata ${i} INMA.commonControlSamples.Rdata
done

Rscript bin/createGenomicRatioSets.R INMA_comb.normalizedComBat.GenomicRatioSet.corrected.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA_comb
Rscript bin/createGenomicRatioSets.R INMA0combined.normalizedComBat.GenomicRatioSet.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA0combined

## Medall
Rscript workflows/bin/createGenomicRatioSets.R MeDALL_all.normalizedComBat.GenomicRatioSet.Rdata MeDALL_all.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds MeDALL_all

## No ComBat
Rscript bin/createGenomicRatioSets_Raw.R INMA_comb.normalizedRaw.GenomicRatioSet.corrected.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA_comb.normalizedRaw
Rscript bin/createGenomicRatioSets_Raw.R MeDALL_all.normalizedRaw.GenomicRatioSet.Rdata MeDALL_all.qc.objects.clean.Rdata data/HM450.hg19.manifest.rds MeDALL_all.normalizedRaw
Rscript bin/createGenomicRatioSets_Raw.R INMA0combined.normalizedRaw.GenomicRatioSet.Rdata INMA_comb.qc.objects.round1.Rdata data/HM450.hg19.manifest.rds INMA0combined.normalizedRaw

### GSE168739
Rscript workflows/bin/createGenomicRatioSets_Raw.R results/preprocess/GSE168739/GSE168739.normalizedRaw.GenomicRatioSet.Rdata results/preprocess/GSE168739/GSE168739.qc.objects.round1.Rdata data/EPIC.hg19.manifest.rds GSE168739.normalizedRaw
Rscript workflows/bin/createGenomicRatioSets_Raw.R results/preprocess/GSE168739/GSE168739.normalizedComBat.GenomicRatioSet.Rdata results/preprocess/GSE168739/GSE168739.qc.objects.round1.Rdata data/EPIC.hg19.manifest.rds GSE168739.normalizedComBat

### GSE87650
Rscript workflows/bin/createGenomicRatioSets_Raw.R results/preprocess/GSE87650/GSE87650.wholeblood.GenomicRatioSet.Rdata results/preprocess/GSE87650/GSE87650.wholeblood.qc.objects.round4.Rdata data/HM450.hg19.manifest.rds GSE87650.wholeblood

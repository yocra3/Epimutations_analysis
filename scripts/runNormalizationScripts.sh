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

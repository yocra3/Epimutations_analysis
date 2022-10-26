#'#################################################################################
#'#################################################################################
# Set project folders' structure and link files
## Base folder -> $HOME/data/WS_HELIX/HELIX_analyses/epimutations_CR
#'#################################################################################
#'#################################################################################

## Make folders
mkdir data
mkdir results
mkdir results/preprocess
mkdir results/preprocess/phenotypes
mkdir results/preprocess/methylation
mkdir results/candidateRegions/

mkdir reports

## Link raw data
### INMA 
#### IDATs
ln -s $HOME/data/WS_INMA/Methylation_INMA/450K_blood/DATA/IDAT data/IDATs_MeDALL
ln -s $HOME/data/WS_INMA/Methylation_INMA/450K_blood/DATA2/IDAT/ data/IDATs_Esteller

#### Phenotypes
ln -s $HOME/data/WS_INMA/Methylation_INMA/450K_blood/DATA2/samples/BD_INMASab_metilacio.dta data/INMA_meth_pheno.dta

#### Genotypes
ln -s $HOME/data/WS_INMA/GWAS_INMA/QC/IMPUTATION_HRC_chrX_Michigan/ImpPLINK_QC data/INMA_genos

### HELIX
ln -s $HOME/data/WS_HELIX/HELIX_preproc/methylation/Final_data/methylome_subcohort_ComBatSlide_6cells_v4.Rdata data/HELIX.GenomicRatioSet.Rdata
ln -s $HOME/data/WS_HELIX/HELIX_preproc/methylation/QC_CRA/detectionvals.RData data/HELIX.detectionPvals.Rdata
ln -s /PROJECTES/HELIX_OMICS/data_final/trans/child/8y/blood_HTAv2_QChelix_20180701/transcriptome_subcohort_f1_v3.RData data/HELIX.genexp.Rdata
ln -s /PROJECTES/INMA/inma_expset_all_20180615.rda data/INMA4.genexp.Rdata
ln -s $HOME/data/WS_HELIX/HELIX_preproc/methylation/Final_data/methylome_subcohort_v4.RData data/HELIX.noCombat.GenomicRatioSet.Rdata

ln -s $HOME/data/WS_HELIX/HELIX_preproc/exposome/FinalDataset/imppostnatal_v3.Rdata data/postExposome.Rdata
ln -s $HOME/data/WS_HELIX/HELIX_preproc/exposome/FinalDataset/imppregnancy_v3.Rdata data/pregExposome.Rdata

## Add eQTM catalogue
ln -s ../expr_met_SM/HELIX_blood_eQTM_WebCat/eQTM_autosome_adj.cells_SIG.txt.gz data/eqtm.txt.gz

## Download CREs from ENCODE
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed

## Download epimutations list
wget https://ars.els-cdn.com/content/image/1-s2.0-S0002929720302883-mmc4.xlsx
mv 1-s2.0-S0002929720302883-mmc4.xlsx Epimutations.PMID32937144.xlsx
chmod 444 Epimutations.PMID32937144.xlsx

## Extract GSE168739
mkdir data/GSE168739
tar -xvf data/GSE168739_RAW.tar -C data/GSE168739

## Extract GSE51032
mkdir data/GSE51032
tar -xvf data/GSE51032_RAW.tar -C data/GSE51032

## Extract GSE112611
mkdir data/GSE112611
tar -xvf data/GSE112611_RAW.tar -C data/GSE112611

rsync -azvh --progress . ~/data/WS_HELIX/HELIX_analyses/epimutations_CR/
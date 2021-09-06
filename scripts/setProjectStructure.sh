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



## Download CREs from ENCODE
wget https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed

## Download epimutations list
wget https://ars.els-cdn.com/content/image/1-s2.0-S0002929720302883-mmc4.xlsx
mv 1-s2.0-S0002929720302883-mmc4.xlsx Epimutations.PMID32937144.xlsx
chmod 444 Epimutations.PMID32937144.xlsx

rsync -azvh --progress . ~/data/WS_HELIX/HELIX_analyses/epimutations_CR/
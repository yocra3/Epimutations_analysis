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
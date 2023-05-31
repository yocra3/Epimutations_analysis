# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("tidyverse", "cowplot", "ramr"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("dgarrimar/mlm")' ## Epimutations dependency
RUN R -e 'devtools::install_github("isglobal-brge/epimutacions@e176598")' ## Epimutations

# Install required Bioconductor packages for revision 2
RUN R -e 'BiocManager::install(c("ChAMP", "DMRcate"))'
RUN git clone https://github.com/brentp/combined-pvalues && \
  cd combined-pvalues && \
  python3 setup.py install
  
RUN apt-get update && \
    apt-get install -y bedtools


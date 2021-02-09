# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install required Bioconductor packages
RUN R -e 'BiocManager::install("minfi")'

# Install tidyverse
RUN R -e 'install.packages("tidyverse")'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
RUN R -e 'devtools::install_github("isglobal-brge/EpiMutations:develop")' ## Epimutations

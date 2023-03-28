# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("tidyverse", "cowplot", "ramr"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("dgarrimar/mlm")' ## Epimutations dependency
RUN R -e 'devtools::install_github("isglobal-brge/epimutacions@e176598")' ## Epimutations

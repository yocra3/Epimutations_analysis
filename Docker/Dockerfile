# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("minfi", "Gviz", "a4Base", "tidyverse"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
RUN R -e 'devtools::install_github("dgarrimar/mlm")' ## Epimutations dependency
RUN R -e 'devtools::install_github("isglobal-brge/EpiMutations", ref = "develop")' ## Epimutations

# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("minfi", "Gviz", "a4Base", "tidyverse"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
RUN R -e 'devtools::install_github("dgarrimar/mlm")' ## Epimutations dependency
RUN R -e 'devtools::install_github("isglobal-brge/EpiMutations", ref = "develop")' ## Epimutations

# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
ENV PATH /opt/conda/bin:$PATH

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/epimutations/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name epimutations > epimutations.yml


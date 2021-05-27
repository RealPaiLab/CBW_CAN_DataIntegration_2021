# create Rstudio docker with supporting packages to suit CBW 2021 AWS environment
FROM rocker/rstudio

RUN apt-get update && apt-get install -y git vim
RUN apt-get install -y libxml2 zlib1g-dev libglpk40
RUN apt-get install openjdk-8-jre-headless

# install BioConductor
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")';
RUN R -e 'BiocManager::install(version="3.12")';
#install additional packages
RUN install2.r \
   -d TRUE \
   -r "https://cran.rstudio.com" \
   ROCR \
   pracma \
   ggplot2 \
   glmnet \
   igraph \
   reshape2 \
   parallel \ 
   bigmemory \
   doParallel \
   foreach \
   combinat \
   rappdirs \
   RColorBrewer \
   Rtsne \
   rmarkdown \
   knitr

# Install required Bioconductor package
RUN R -e 'BiocManager::install(c("RCy3","biomaRt","MultiAssayExperiment","curatedTCGAData","BiocFileCache","GenomicRanges","GenomeInfoDb","IRanges","scater","clusterExperiment","netSmooth","BiocStyle","TCGAutils"))'

RUN git clone https://github.com/BaderLab/netDx.git

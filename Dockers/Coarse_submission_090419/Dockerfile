## Start from this Docker image
FROM rocker/tidyverse:3.6.0

## Install R packages in Docker image

## RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
## RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

## RUN Rscript -e "install.packages('BiocManager')"
## RUN Rscript -e "BiocManager::install('WGCNA')"

RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('tidyr')"
RUN Rscript -e "install.packages('readr')"
RUN Rscript -e "install.packages('purrr')"
RUN Rscript -e "install.packages('tibble')"
RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('glmnet')"

## install.packages("BiocManager") 
## BiocManager::install("WGCNA") 

## Copy your files into Docker image
COPY bin/* ./

## Make script executable
RUN chmod a+x decon_prediction.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "decon_prediction.R"]
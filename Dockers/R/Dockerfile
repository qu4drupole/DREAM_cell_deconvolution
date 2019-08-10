## Start from this Docker image
FROM rocker/tidyverse:3.6.0

## Install R packages in Docker image

## RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
## RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('WGCNA')"

## install.packages("BiocManager") 
## BiocManager::install("WGCNA") 

## Copy your files into Docker image
COPY bin/docker_test.R ./

## Make script executable
RUN chmod a+x docker_test.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "docker_test.R"]
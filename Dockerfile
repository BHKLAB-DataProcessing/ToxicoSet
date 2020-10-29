# hash:sha256:5211f36fab9b514f0d1ad9005f3a612dea08555a685ef0a15c24272656504e8a
FROM registry.codeocean.com/codeocean/r-base:4.0.0-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_github( \
        "bhklab/CoreGx", \
        ref = "febf598cb579333e8bc1f1c5fdad6ac30b709416")' \
    && Rscript -e 'remotes::install_github( \
        "bhklab/ToxicoGx", \
        ref = "5cc8641e454ef33c93aabf5ab4ad226f9a9280f3")'

RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools); install_version('BiocManager', version = '1.30.10', repos = 'http://cran.us.r-project.org')"

RUN Rscript -e "BiocManager::install(version='3.11', ask = FALSE)"
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "Biobase", \
        "BiocManager", \
        "PharmacoGx", \
        "SummarizedExperiment", \
        "abind", \
        "affy", \
        "affyio", \
        "biomaRt", \
        "car", \
        "data.table", \
        "dplyr", \
        "gdata", \
        "ggplot2", \
        "readxl", \
        "xml2" \
    ))' # Original versions: 2.48.0 1.30.10 2.0.1 1.18.1 1.4-5 1.66.0 1.58.0 2.44.0 3.0-8 1.12.8 0.8.5 2.18.0 3.3.0 1.3.1 1.3.2

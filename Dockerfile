FROM rocker/r-base

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libcurl4-openssl-dev \
    pandoc \
    python-pip \
    scons \
    ssh
RUN install2.r --error --deps TRUE \
    ape \
    expm \
    ggplot2 \
    jsonlite \
    mcmc \
    optparse \
    pander \
    rmarkdown \
    stringdist \
    && rm -rf /tmp/downloaded_packages/
RUN install2.r -r http://bioconductor.org/packages/3.0/bioc --deps TRUE \
    BiocInstaller \
    Biostrings \
    && rm -rf /tmp/downloaded_packages/
RUN pip install nestly

COPY . /context
WORKDIR /context
CMD ./build-and-test.sh

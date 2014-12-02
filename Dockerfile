FROM rocker/r-base

RUN apt-get update && apt-get install -y --no-install-recommends \
    python-pip \
    libcurl4-openssl-dev

RUN install2.r --error --deps TRUE \
    expm \
    mcmc \
    stringdist \
    optparse \
    jsonlite \
    ape \
    rmarkdown \
    ggplot2 && \
    rm -rf /tmp/downloaded_packages/

RUN pip install nestly

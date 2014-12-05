FROM rocker/r-base

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libcurl4-openssl-dev \
    python-pip \
    ssh

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

RUN git clone git@github.com:petrelharp/context.git
WORKDIR /data/context/nestly
RUN scons simple
RUN scons seed

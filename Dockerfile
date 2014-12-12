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
    IRanges \
    && rm -rf /tmp/downloaded_packages/
RUN pip install nestly


# set up auth
RUN mkdir -p /root/.ssh \
    && chmod 700 /root/.ssh
ADD bunnyhutch_id_rsa /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa \
    && ssh-keyscan github.com >> /root/.ssh/known_hosts

# run!
CMD git clone git@github.com:petrelharp/context.git \
    && cd /context/nestly \
    && scons simple \
    && scons -j 6 seed \
    && cd ../json-cpg && ./minimal-workflow.sh \
    && cd ../json-tasep && ./workflow.sh \
    && cd ../json-ising && ./minimal-workflow.sh

FROM rocker/r-base

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libcurl4-openssl-dev \
    python-pip \
    scons \
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

# set up auth
RUN mkdir -p /root/.ssh && \
    chmod 700 /root/.ssh
ADD bunnyhutch_id_rsa /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa && \
    ssh-keyscan github.com >> /root/.ssh/known_hosts

# run!
RUN git clone git@github.com:petrelharp/context.git
WORKDIR /data/context/nestly
RUN scons simple
RUN scons seed

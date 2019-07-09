FROM ubuntu:18.04

LABEL maintainer="well309@gmail.com"

ENV ADMIXTURE_VERSION 1.3.0

RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
    ca-certificates \
    curl
RUN curl -O http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-${ADMIXTURE_VERSION}.tar.gz
RUN tar xzvf admixture_linux-${ADMIXTURE_VERSION}.tar.gz
RUN cp admixture_linux-${ADMIXTURE_VERSION}/admixture /usr/local/bin/
RUN rm -r admixture_linux*
RUN apt-get remove --yes \
    ca-certificates \
    curl
RUN apt-get autoremove --yes
RUN rm -rf /var/lib/apt/lists/*

CMD admixture
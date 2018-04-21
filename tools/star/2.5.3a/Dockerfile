FROM debian:stretch

LABEL maintainer="well309@gmail.com"

ENV STAR_VERSION 2.5.3a

RUN apt-get update && \
  apt-get install --yes --no-install-recommends \
    ca-certificates \
    g++ \
    make \
    wget \
    zlib1g-dev && \
  wget -q https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
  tar xf ${STAR_VERSION}.tar.gz && \
  cd STAR-${STAR_VERSION}/source && \
  make STAR && \
  make STARlong && \
  mv STAR STARlong /usr/local/bin/ && \
  cd ../.. && \
  rm -r ${STAR_VERSION}.tar.gz STAR-${STAR_VERSION}

CMD STAR

FROM broadinstitute/picard

LABEL maintainer="well309@gmail.com"

ENV PICARD_VERSION 2.6.0

RUN apt-get update && \
  apt-get install --yes --no-install-recommends \
    ca-certificates \
    curl && \
  curl -fsSL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -o /usr/picard/picard.jar && \
  apt-get remove --yes \
    curl \
    ca-certificates && \
  apt-get autoremove --yes && \
  rm -rf /var/lib/apt/lists/*

CMD [ "picard" ]
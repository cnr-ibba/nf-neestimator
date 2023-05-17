#
# VERSION 0.2
# DOCKER-VERSION  23.0.4
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    Nestimator v2.x image
# TO_BUILD:       docker build --rm -t bunop/neestimator2x .
# TO_TAG:         docker tag bunop/neestimator2x:latest bunop/neestimator2x:0.2
#

# Get the last rocker/shiny image
FROM debian:11-slim

# The mainteiner of this image
LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it"

# install package dependencies
RUN apt-get update && apt-get install -y \
        procps \
    && apt-get clean && rm -rf /var/lib/apt/lists/

ADD Ne2-1L /usr/local/bin/Ne2-1L

RUN chmod 755 /usr/local/bin/Ne2-1L

CMD ["bash"]

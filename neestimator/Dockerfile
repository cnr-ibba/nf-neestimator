#
# VERSION 0.1
# DOCKER-VERSION  23.0.4
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    Nestimator v2.x image
# TO_BUILD:       docker build --rm -t bunop/neestimator2x .
# TO_TAG:         docker tag bunop/neestimator2x:latest bunop/neestimator2x:0.1
#

# Get the last rocker/shiny image
FROM debian:11-slim

# The mainteiner of this image
LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it"

# install package dependencies
RUN apt-get update && apt-get install -y \
        procps \
        build-essential \
        git \
    && apt-get clean && rm -rf /var/lib/apt/lists/

WORKDIR /opt

# clone NeEstimator software
RUN git clone https://github.com/bunop/NeEstimator2.X.git \
    && cd NeEstimator2.X \
    && git checkout issue-1

WORKDIR /opt/NeEstimator2.X

# compile source code and install
RUN make
RUN ln -s /opt/NeEstimator2.X/Ne2x /usr/local/bin/Ne2-1L

CMD ["bash"]

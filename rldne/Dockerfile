#
# VERSION 0.2
# DOCKER-VERSION  23.0.4
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    A docker r-base container with RLDNe installed
# TO_BUILD:       docker build --rm -t bunop/rldne .
# TO_TAG:         docker tag bunop/rldne:latest bunop/rldne:0.2
#

# Get the last rocker/shiny image
FROM r-base:4.3.0

# The mainteiner of this image
LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it"

# install package dependencies
RUN apt-get update && apt-get install -y \
        procps \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libcurl4-openssl-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/

# Install R packages. Packages are listed in packages.txt
COPY install.r /root/
RUN Rscript --slave --vanilla /root/install.r

CMD ["R"]

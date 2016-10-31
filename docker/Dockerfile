FROM continuumio/anaconda3:4.1.1
MAINTAINER "Greg Way" <gregway@mail.med.upenn.edu>

# Install base packages and R
RUN apt-get update --fix-missing && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libssh2-1-dev \
    vim \
    r-base \
    r-base-dev \
    r-cran-rcpp

# Update to R 3.3.1
RUN echo "deb http://cran.irsn.fr/bin/linux/debian jessie-cran3/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480
RUN apt-get update && apt-get install -y r-base

# Clone TAD pathways GitHub repository
RUN git clone https://github.com/greenelab/tad_pathways.git  
RUN mkdir ~/.checkpoint


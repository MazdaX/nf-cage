FROM alpine:3.11 
ARG  ARIA2_VER=1.35.0
RUN  apk add --no-cache ca-certificates make g++ gcc file cppunit-dev zlib-dev openssl-dev expat-dev sqlite-dev c-ares-dev libssh2-dev \
	&&   mkdir /aria2c  \
	&&   wget -P /aria2c   https://github.com/aria2/aria2/releases/download/release-${ARIA2_VER}/aria2-${ARIA2_VER}.tar.gz  \
	&&   tar  -zxvf  /aria2c/aria2-${ARIA2_VER}.tar.gz  -C    /aria2c  \
	&&   cd  /aria2c/aria2-${ARIA2_VER}  \
	&&   sed  -i  's/"1", 1, 16/"1", 1, 128/g'          src/OptionHandlerFactory.cc    \
	&&   sed  -i  's/"20M", 1_m, 1_g/"4k", 1_k, 1_g/g'  src/OptionHandlerFactory.cc    \
	&&   ./configure --without-libxml2  --without-gnutls --with-openssl  --host=x86_64-alpine-linux-musl   \
	&&   make install-strip   \
	&&   ldd /usr/local/bin/aria2c   |cut -d ">" -f 2|grep lib|cut -d "(" -f 1|xargs tar -chvf /aria2c/aria2c.tar  \
	&&   mkdir /aria2  \
	&&   tar  -xvf /aria2c/aria2c.tar  -C /aria2   \
	&&   cp --parents /usr/local/bin/aria2c /aria2


FROM continuumio/miniconda3
#FROM ubuntu:xenial  
#FROM r-base
#EXPOSE 5000
#EXPOSE 443

LABEL \
  author="Mazdak Salavati" \
  description="nf-cage base image for use in nf-core pipelines" \
  maintainer="Mazdak.Salavati@roslin.ed.ac.uk"

# Essentials
RUN apt-get update && \
	apt install -y procps wget gzip pigz && \
	apt-get autoclean && \
	apt-get autoremove && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda && \
	conda config --add channels default

# Setting up the nf_DSL2
COPY environment.yml /
RUN conda install --yes nomkl mamba -c conda-forge && conda clean -afy
RUN mamba env create -n nf_DSL2 -f /environment.yml && conda clean -afy
ENV PATH /opt/conda/envs/nf_DSL2/bin:$PATH


#adding modules to the docker image
COPY modules /modules

# For posterity debugging
RUN conda env export -n nf_DSL2 > /nf_DLS2_docker.yml

#Aria2c from the first build
COPY --from=0 /aria2 /
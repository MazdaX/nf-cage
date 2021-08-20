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
RUN conda env create -n nf_DSL2 && conda install -n nf_DSL2 mamba -c conda-forge && conda clean -a
SHELL ["conda", "activate", "-n", "nf_DSL2", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "nf_DSL2"]
RUN mamba install -n nf_DSL2 --file /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf_DSL2/bin:$PATH

#adding modules to the docker image
COPY modules /modules

# For posterity debugging
RUN conda env export -n nf_DSL2 > /nf_DLS2_docker.yml


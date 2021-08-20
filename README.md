# BovReg nf-cage pipeline 
This is a CAGE analysis pipeline used in BovReg consoritum (https://www.bovreg.eu/) studies. This pipeline was wrapped using NextFlow DSL2 syntax from demultiplexing to base pair resolution and strand specific read counting <br> (i.e. compatible for import to CAGEfightR https://bioconductor.org/packages/release/bioc/html/CAGEfightR.html)

Pipeline has the following steps: 
1. Demultiplexing raw CAGE sequence data using _FASTXtoolkit_ (allowed mismatch 1)
2. Trimming CAGE tags using _tagDust2_ (HMM read architechture provided)
3. QC before and after (_FastQC_ and _MultiQC_)
4. Mapping against reference genome (i.e. ARS-UCD1.2_btau5_0.1_Y 1000bull project) using bowtie2
5. BAM >>> bp resolution bedGraph for +ve and -ve strands >>>> bigWig

# Installation
You would need to use the docker container for reproducibility of the analysis. Please follow the installation guidelines of Docker on your system : https://docs.docker.com/get-docker/<br>
or the desktop client (installs the engine and GUI as the same time) https://docs.docker.com/desktop/<br>
Download the code repo and the respective docker image to your system (Tested on WSL2 Windows 10 and Linux Ubuntu 18.04 LTS) 

```
git clone https://github.com/MazdaX/nf-cage.git
docker pull mazdax/nf-cage:minimal

```
_NB. Running this pipeline without docker would require modification of module nf scripts and its not recommended_

# Quick Start

Sample FASTQ files are in the root folder (fastq_files) along with the barcodes per samples (barcode_files). Replace the files inside these 2 folders with your own experimental data in order to run the pipeline on your dataset. 

```
nextflow run . --with-docker mazdax/nf-cage:minimal 

```

# Container and toolset
This pipeline uses a docker container for all the tools required and the mamba environment. Please find the details at Docker hub public repository mazdax/nf-cage:
```
docker pull mazdax/nf-cage
```
The list of the tools and versions are available in the __environment.yml__ and the __Dockerfile__. 





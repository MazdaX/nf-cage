<p float="right">
  <img align="left" width="300" ![BovReg Logo] src="/images/BV_logo.png">
  <img align="right" width="200" ![NextFlow Logo] src="/images/NF_logo.png">  <br><br><br><br><br>
  <img align="right" width="180" ![Docker Logo] src="/images/docker_logo.png">
</p>

# BovReg nf-cage pipeline 
This is a CAGE analysis pipeline used in BovReg consortium (https://www.bovreg.eu/) studies. This pipeline was wrapped using NextFlow DSL2 syntax from demultiplexing to base pair resolution and strand specific read counting <br> (i.e. compatible for import to CAGEfightR https://bioconductor.org/packages/release/bioc/html/CAGEfightR.html)

Pipeline has the following steps: 
1. Demultiplexing raw CAGE sequence data using _FASTXtoolkit_ (allowed mismatch 1)
2. Trimming CAGE tags using _tagDust2_ (HMM read architecture provided) 
3. QC before and after (_FastQC_ and _MultiQC_)
4. Mapping against reference genome [ARS-UCD1.2_Btau5.0.1Y 1000 bull project](https://sites.ualberta.ca/~stothard/1000_bull_genomes/) using bowtie2
5. BAM >>> bp resolution bedGraph for +ve and -ve strands >>>> bigWig


# Installation

Install the latest NextFlow engine on your system following the instructions available at https://www.nextflow.io/<br>
After a successful install you should be able to query the following without any errors: 

```
~/nf-cage$ nextflow info
  Version: 21.04.0 build 5552
  Created: 02-05-2021 16:22 UTC (17:22 BST)
  System: Linux 5.10.43.3-microsoft-standard-WSL2
  Runtime: Groovy 3.0.7 on OpenJDK 64-Bit Server VM 11.0.9.1-internal+0-adhoc..src
  Encoding: UTF-8 (UTF-8)
```

You would need to use the docker container for reproducibility of the analysis. Please follow the installation guidelines of Docker on your system : https://docs.docker.com/get-docker/ or the desktop client (installs the engine and GUI as the same time) https://docs.docker.com/desktop/<br>

Download this code repository and the respective docker image (Tested on WSL2 Windows 10 and Linux Ubuntu 18.04 LTS) 

```

git clone https://github.com/MazdaX/nf-cage.git
docker pull mazdax/nf-cage:latest

```
_NB. Running this pipeline without docker would require modification of module nf scripts and its not recommended_

# Quick Start

Sample FASTQ files are in the root folder (fastq_files) along with the barcodes per samples (barcode_files). Replace the files inside these 2 folders with your own experimental data in order to run the pipeline on your dataset. 

```
nextflow run . --with-docker mazdax/nf-cage:latest

```

# Container and toolset
This pipeline uses a docker container for all the tools required and the mamba environment. Please find the details at Docker hub public repository https://hub.docker.com/r/mazdax/nf-cage :
```
docker pull mazdax/nf-cage:latest
```
The list of the tools and versions are available in the __environment.yml__ and the __Dockerfile__. 
The aria2c compilation was modified using https://registry.hub.docker.com/r/johngong/aria2/dockerfile .

```
docker build -t mazdax/nf-cage:aria2c .
```





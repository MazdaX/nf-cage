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
nextflow run . -with-docker mazdax/nf-cage:latest 

```

or 

```
nextflow run -profile docker .
N E X T F L O W  ~  version 21.04.3
Launching `./main.nf` [cranky_cantor] - revision: 9e652e725b

===============================================
        nf-cage BovReg's pipeline
        github.com/mazdax/nf-cage
        docker pull mazdax/nf-cage
===============================================
Input : Raw Illumina CAGE sequences
Input : Barcode list TSV (i.e. sample   barcode)
Input : Reference Genome FASTA (bowtie2 index)
Output : Strand specific bp resolution bigWig 
Running task: AIO CAGEfightR import ready
-----------------------------------------------

executor >  local (103)
[9b/f36040] process > qc_pre (FastQC...)                      [100%] 2 of 2 ✔
[f2/a8d4ee] process > demux (Demultiplexing...)               [100%] 2 of 2 ✔
[53/7323fa] process > merger (Gathering demux...)             [100%] 20 of 20 ✔
[37/e53393] process > qc_post (FastQC...)                     [100%] 20 of 20 ✔
[35/4590e7] process > multiqc (MultiQC...)                    [100%] 1 of 1 ✔
[87/7af88a] process > trimmer (Trimming by TagDust2...)       [100%] 19 of 19 ✔
[26/c6d571] process > mapKeeper (Sourcing the reference...)   [100%] 1 of 1 ✔
[0d/b0ce36] process > mapper (Mapping using bowtie2...)       [100%] 19 of 19 ✔
[c6/e7ddfb] process > bG2bW (BAM >>> bedGraph >>> BigWig ...) [100%] 19 of 19 ✔
Completed at: 04-Sep-2021 13:42:08
Duration    : 1d 2h 35m 41s
CPU hours   : 212.7
Succeeded   : 103

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

NB. nf-core and pyinquirer packages were manually removed from the docker builds. 
NB. The main branch contains tested pipeline and for nf-core compatibility a "nfcore_opt" branch was created




<p float="right">
  <img align="left" width="300" ![BovReg Logo] src="/images/BV_logo.png">
  <img align="right" width="200" ![NextFlow Logo] src="/images/NF_logo.png">  <br><br><br><br><br>
  <img align="right" width="180" ![Docker Logo] src="/images/docker_logo.png">
</p>

[![DOI](https://zenodo.org/badge/397618218.svg)](https://zenodo.org/badge/latestdoi/397618218)

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
nextflow info
  Version: 22.04.3 build 5703
  Created: 18-05-2022 19:22 UTC (20:22 BST)
  System: Linux 5.10.102.1-microsoft-standard-WSL2
  Runtime: Groovy 3.0.10 on OpenJDK 64-Bit Server VM 11.0.9.1-internal+0-adhoc..src
  Encoding: UTF-8 (UTF-8)
```

You would need to use the docker container for reproducibility of the analysis. Please follow the installation guidelines of Docker on your system : https://docs.docker.com/get-docker/ or the desktop client (installs the engine and GUI as the same time) https://docs.docker.com/desktop/<br>

Download this code repository and the respective docker image (Tested on WSL2 Windows 10 and Linux Ubuntu 18.04 LTS) 

```

git clone https://github.com/MazdaX/nf-cage.git
docker pull mazdax/nf-cage:latest

```
_NB. Running this pipeline without docker would require modification of module nf scripts and its not recommended_
_NB. The docker image will be converted to a singularity sif using the -profile singularity _

# Quick Start

The help manu can be access as following: 

```
nextflow run . --help

N E X T F L O W  ~  version 22.04.3
Launching `./main.nf` [nasty_bose] DSL2 - revision: 17b2193c3b

======================================================================
         

        ███╗   ██╗███████╗       ██████╗ █████╗  ██████╗ ███████╗
        ████╗  ██║██╔════╝      ██╔════╝██╔══██╗██╔════╝ ██╔════╝
        ██╔██╗ ██║█████╗  █████╗██║     ███████║██║  ███╗█████╗  
        ██║╚██╗██║██╔══╝  ╚════╝██║     ██╔══██║██║   ██║██╔══╝  
        ██║ ╚████║██║           ╚██████╗██║  ██║╚██████╔╝███████╗
        ╚═╝  ╚═══╝╚═╝            ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝                                                         

======================================================================
    nf-cage pipeline
    github.com/mazdax/nf-cage
    docker pull mazdax/nf-cage
    Written by Mazdak Salavati (Twitter @MazdakS)

This is a CAGE analysis pipeline used in the BovReg consortium project (https://www.bovreg.eu/).
This pipeline was wrapped using NextFlow DSL2 syntax from demultiplexing to base pair resolution 
and strand specific read counting. The output for each single end FASTQ file are 2 bedGraph
(+,-strands) bp resolution CAGE tag counts. The bedGraph outputs can be directly used in the CAGEfightR
package (https://bioconductor.org/packages/release/bioc/html/CAGEfightR.html)

Pipeline has the following steps:

- Demultiplexing raw single end CAGE sequence data using FASTXtoolkit (allowed mismatch 1)
- Trimming CAGE tags using tagDust2 
(HMM read architecture provided: -1 B:${barcode} -2 F:CAGNNNG -3 R:N -4 P:ATCTCGTATGCCGTCTTCTGCTT -dust 100)
- QC before and after (FastQC and MultiQC)
- Mapping against reference genome ARS-UCD1.2_Btau5.0.1Y 1000 bull project using bowtie2
- BAM >>> bp resolution bedGraph for +ve and -ve strands >>>> bigWig
======================================================================
      

    Pipeline version: 0.0.1
    
    Typical usage:
            nextflow run . -profile singularity 
    
    Parameters:
    
    --help              Prints this message

    --ref_fasta         URL address to the reference FASTA File 
                        (e.g. "http://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz")

    --raw_fastq         Address of the folder containing gz-compressed CAGE (SE libraries before demux) FASTQ files
                        (e.g. --raw_fastq fastq_files/*.gz)

    --barcodes          Tab separated file with no header for the input barcodes to be used in demux 
                        (i.e. sampleName\tbarcode)
                        
    --allowed_mismatch  The number of barcode mismatches allowed by the fastx_barcode_splitter.pl script
                         (default: 1)
```


Sample FASTQ files are in the root folder (fastq_files) along with the barcodes per samples (barcode_files). Replace the files inside these 2 folders with your own experimental data in order to run the pipeline on your dataset. 

```
nextflow run . -profile singularity

```

or 

```
nextflow run . -c configs/conf/eddie.config

N E X T F L O W  ~  version 22.04.3
Launching `./main.nf` [silly_davinci] DSL2 - revision: 1908217c19

===============================================
        nf-cage BovReg's pipeline
        github.com/mazdax/nf-cage
        docker pull mazdax/nf-cage
===============================================
Input : Raw Illumina CAGE sequences
Input : Barcode list TSV (i.e. sample   barcode)
Input : Reference Genome FASTA (URL or local)
Output : Strand specific bp resolution bigWig 
Running task: AIO CAGEfightR import ready
Workflow version: 0.0.1
-----------------------------------------------

executor >  sge (19)
[d6/d0ed33] process > qc_pre (Running FastQC...)              [100%] 2 of 2, cached: 2 ✔
[c9/2f0d69] process > DEMUX (Demultiplexing...)               [100%] 2 of 2, cached: 2 ✔
[46/301625] process > MERGER (Gathering demux...)             [100%] 20 of 20, cached: 20 ✔
[a4/d52782] process > qc_post (Running FastQC...)             [100%] 20 of 20, cached: 20 ✔
[28/5435f1] process > MULTIQC (Running MultiQC...)            [100%] 1 of 1, cached: 1 ✔
[d9/a7e4d3] process > TRIMMER (Trimming by TagDust2...)       [100%] 19 of 19, cached: 19 ✔
[dc/9dce8d] process > DOWNLOADREF (Sourcing the reference...) [100%] 1 of 1, cached: 1 ✔
[5a/414702] process > BT2BUILD (Bowtie2 build...)             [100%] 1 of 1, cached: 1 ✔
[f6/839e10] process > BT2MAPPER (Mapping using bowtie2...)    [100%] 19 of 19, cached: 19 ✔
[10/9c29e1] process > BG2BW (BAM >>> bedGraph >>> BigWig ...) [100%] 19 of 19 ✔
Completed at: 08-Jun-2022 12:01:09
Duration    : 3h 7m 58s
CPU hours   : 216.0 (2.7% cached)
Succeeded   : 19
Cached      : 85

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
NB. The main branch contains tested pipeline and for nf-core compatibility a "nfcore_opt" branch was created.
NB. The dev and dev_hpc branches were merged in to main for the final version. 

# Citation

Salavati M, Caulton A, Clark R, Gazova I, Smith TPL, Worley KC, Cockett NE, Archibald AL, Clarke SM, Murdoch BM, Clark EL. Global Analysis of Transcription Start Sites in the New Ovine Reference Genome (Oar rambouillet v1.0). Front Genet. 2020 Oct 23;11:580580. doi: 10.3389/fgene.2020.580580. PMID: 33193703; PMCID: PMC7645153.


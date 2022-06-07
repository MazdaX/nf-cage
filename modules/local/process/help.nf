#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl = 2


def helpMessage() {
  log.info '''
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

This is a CAGE analysis pipeline used in the BovReg consortium project (https://www.bovreg.eu/). This pipeline was wrapped using NextFlow DSL2 syntax from
demultiplexing to base pair resolution and strand specific read counting. The output for each single end FASTQ file are 2 bedGraph/bigWig (+,-strands) bp resolution 
CAGE tag counts. The bedGraph/bigWig outputs can be directly used in the CAGEfightR package (https://bioconductor.org/packages/release/bioc/html/CAGEfightR.html)

Pipeline has the following steps:

- Demultiplexing raw single end CAGE sequence data using FASTXtoolkit (allowed mismatch 1)
- Trimming CAGE tags using tagDust2 (HMM read architecture provided: -1 B:${barcode} -2 F:CAGNNNG -3 R:N -4 P:ATCTCGTATGCCGTCTTCTGCTT -dust 100)
- QC before and after (FastQC and MultiQC)
- Mapping against reference genome ARS-UCD1.2_Btau5.0.1Y 1000 bull project using bowtie2
- BAM >>> bp resolution bedGraph for +ve and -ve strands >>>> bigWig
======================================================================
      '''
    log.info"""
    Pipeline version: ${workflow.manifest.version}
    
    Typical usage:
            nextflow run . -profile singularity 
    
    Parameters:
    
    --help              Prints this message
    --ref_fasta         URL address to the reference FASTA File (e.g. "http://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz")
    --raw_fastq         Address of the folder containing gz-compressed CAGE (SE libraries before demux) FASTQ files (e.g. --raw_fastq fastq_files/*.gz)
    --barcodes          Tab separated file with no header for the input barcodes to be used in demux (i.e. sampleName\\tbarcode)
    --allowed_mismatch  The number of barcode mismatches allowed by the fastx_barcode_splitter.pl script (default: 1)
    """
}
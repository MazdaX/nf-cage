#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

//Should be coded for ARG in the future
params.allowed_mismatch = 1

process DEMUX {
    tag "Demultiplexing..."
    label "small"
    //Redundancy to save space
    //publishDir "$projectDir/demux/${sample_id}" , pattern: "*.fastq", mode: 'copy', overWrite: true
    cache true    

    input:
        tuple val( sample_id ), path( fastq_ch )
        path( barcode )
        
    output:
        path "CAGE_*.fastq" , emit: OUT_demux
                
    script:
    """
        mkdir ${sample_id}
        zcat $fastq_ch |\
        /modules/local/scripts/fastx_0.0.13/bin/fastx_barcode_splitter.pl \
        --bcfile $barcode \
        --bol --mismatch $params.allowed_mismatch \
        --prefix "${sample_id}_" \
        --suffix ".fastq"
        
    """
       
}


process MERGER {
    tag "Gathering demux..."
    label 'small'
    //label "proccess_wsl"
    publishDir "$projectDir/demux/",pattern: "*.fastq.gz", mode: 'copy' , overWrite: true
    cache true
    
    input:
        tuple val(sample_id), path(directory)
    output:
        path( '**.f*q.gz' ), emit: OUT_merger
    script:
    
    //Needs fixing based on the input tuple
    """
    cat ${directory} > ${sample_id}.fastq
    pigz --force ${sample_id}.fastq
            
    """
}
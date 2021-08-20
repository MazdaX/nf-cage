#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

params.out="$projectDir/bams"
params.out2="$projectDir/ref"

process mapKeeper {
    tag "Sourcing the reference..."
    publishDir params.out2 , mode: 'copy', overWrite: true
    cpus 6
    maxForks 100
    cache true

    output:
        path '*ARS-UCD1.2*'
    script:
    """
    mkdir -p /ref

    if [ -f /ref/ARS-UCD1.2.fa ];then
        echo "Reference exists ..."
    else
        echo "Downloading Bos_taurus.ARS-UCD1.2 from Ensembl v103..."
        wget http://ftp.ensembl.org/pub/release-103/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
    fi;

    if [ -f /ref/ARS-UCD1.2.1.bt2 ];then
            echo "Reference exists and indices are in the right folder."
        else
            pigz -d -p 6 Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
            mv Bos_taurus.ARS-UCD1.2.dna.toplevel.fa /ref/ARS-UCD1.2.fa
            samtools faidx /ref/ARS-UCD1.2.fa
            echo "Indexing Bos_taurus.ARS-UCD1.2 for the bowtie2..."
            bowtie2-build --threads 6 /ref/ARS-UCD1.2.fa /ref/ARS-UCD1.2
    fi
    #the 1000bull genome MT is longer than ENSEMBL and this file should be reproduced for the 1KB runs
    awk '{print \$1,\$2+2}' /ref/ARS-UCD1.2.fa.fai > /ref/ref_cov
    #The relative output issue with line 19
    cp -r /ref/ /
    """
}



process mapper {
    tag "Mapping using bowtie2..."
    publishDir params.out , mode: 'copy', overWrite: true
    cpus 6
    maxForks 100
    cache true
    
    input:
        tuple val(name) , path(trimmed_fastq)
    output:
        path '*.bam' , emit: OUT_mapped
        path '*.metrics', emit: OUT_mapped_metrics
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 

    //Docker addresses fir ref has the same problem as the trimmer tagdust from the modules
    //Docker optimisation is needed
    
    script:
    """
        mkdir -p /bams && \
        bowtie2 -p 6 --met-file ${name}.metrics --very-sensitive \
        --rg-id ${name} --rg LB:${name} --rg PL:ILLUMINA --rg SM:${name} \
        -x /ref/ARS-UCD1.2 \
        -U ${trimmed_fastq} | \
        samtools view -@ 6 -bS -F 4 | \
        samtools sort -@ 6 -o ${name}.bam
        samtools index ${name}.bam

    """
       
}
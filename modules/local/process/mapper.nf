#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl = 2

params.out  = "$projectDir/bams"
params.out2 = "$projectDir/"

process downloadRef {
    tag "Sourcing the reference..."
    publishDir "${params.out2}/ref" , mode: 'copy', overWrite: true
    cpus params.all_threads
    maxForks 100
    cache true

    input:
        val fasta
        
    output:
        path '*.fa', emit: fasta

    script:
    """
    #echo "Downloading Bos_taurus.ARS-UCD1.2 from Ensembl v103..."
    #wget $fasta
    aria2c -x 16 $fasta

    pigz -d -p $params.all_threads *.fa.gz
    """
}

process bowtie2Build {
    tag "Bowtie2 build..."
    publishDir "${params.out2}/ref" , mode: 'copy', overWrite: true
    cpus params.all_threads
    maxForks 100
    cache true

    input:
        path fasta
    output:
        path 'ref', emit: bowtie_index
        path 'ref_cov', emit: ref_cov

    script:
    """
    mkdir -p ref

    samtools faidx $fasta

    echo "Indexing ${fasta} for bowtie2..."
    bowtie2-build --threads $params.all_threads ${fasta} ref/${fasta.baseName}
    #the 1000bull genome MT is longer than ENSEMBL and this file should be reproduced for the 1KB runs
    awk '{print \$1,\$2+2}' ${fasta}.fai > ref_cov
    """
}

process mapper {
    tag "Mapping using bowtie2..."
    publishDir params.out , mode: 'copy', overWrite: true
    cpus params.all_threads
    maxForks 100
    cache true
    
    input:
        tuple val(name) , path(trimmed_fastq)
        path bowtie_index
    output:
        path '*.bam' , emit: OUT_mapped
        path '*.metrics', emit: OUT_mapped_metrics
        path '*.bai'
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 

    //Docker addresses for ref has the same problem as the trimmer tagdust from the modules
    //Docker optimisation is needed
    
    script:
    """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`

        mkdir -p bams && \\
        bowtie2 -p $params.all_threads --met-file ${name}.metrics --very-sensitive \\
        --rg-id ${name} --rg LB:${name} --rg PL:ILLUMINA --rg SM:${name} \\
        -x \$INDEX \\
        -U ${trimmed_fastq} | \\
        samtools view -@ $params.all_threads -bS -F 4 | \\
        samtools sort -@ $params.all_threads -o ${name}.bam
        samtools index ${name}.bam
    """
}
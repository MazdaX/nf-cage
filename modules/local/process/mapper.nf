#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl = 2

process DOWNLOADREF {
    tag "Sourcing the reference..."
    publishDir "$projectDir/ref" , mode: 'copy', overWrite: true
    maxForks 100
    cache true

    input:
        val fasta_url
        
    output:
        path '*.fa', emit: fasta
    
    script:
    """
        aria2c -x 16 $fasta_url
        pigz -d -p $task.cpus \$(basename '$fasta_url')
    """
    
}

process BT2BUILD {
    tag "Bowtie2 build..."
    publishDir "$projectDir/ref" , mode: 'copy', overWrite: true
    maxForks 100
    cache true

    input:
        path 'fasta'
    output:
        path 'ref', emit: bowtie_index
        path 'ref_cov', emit: ref_cov

    script:
    """
    mkdir -p ref

    samtools faidx $fasta

    if [ -s "$projectDir/ref/*.bt2" ]; then
	    echo "Bowtie2 indexes already exist"
    else
	    echo "Indexing ${fasta} for bowtie2..."
	    bowtie2-build --threads $task.cpus ${fasta} ref/${fasta.baseName}
    fi
    #the 1000bull genome MT is longer than ENSEMBL and this file should be reproduced for the 1KB runs
    awk '{print \$1,\$2+2}' ${fasta}.fai > ref_cov
    """
}

process BT2MAPPER {
    tag "Mapping using bowtie2..."
    publishDir "$projectDir/bams" , mode: 'copy', overWrite: true
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
    bowtie2 -p $task.cpus --met-file ${name}.metrics --very-sensitive \\
    --rg-id ${name} --rg LB:${name} --rg PL:ILLUMINA --rg SM:${name} \\
    -x \$INDEX \\
    -U ${trimmed_fastq} | \\
    samtools view -@ $task.cpus -bS -F 4 | \\
    samtools sort -@ $task.cpus -o ${name}.bam
    samtools index ${name}.bam
    """
}

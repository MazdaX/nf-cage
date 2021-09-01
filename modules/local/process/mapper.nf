#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl = 2

params.out  = "$projectDir/bams"
params.out2 = "$projectDir/ref"

process mapKeeper {
    tag "Sourcing the reference..."
    publishDir params.out2 , mode: 'copy', overWrite: true
    cpus params.all_threads
    maxForks 100
    cache true
    /*
    docker containers run by NF can NOT operate anything in the root of the instance. The bin and tmp folders that are mounted automatically by the NF are also only good for usage in PATH. The only current solution to mounting external (projectDir) folder inside the containers is to use readonly mounting points within the home folder of the running docker. 
    */

    containerOptions "-v $projectDir/ref:/home/ref:ro"

    //The optional output prevents any unwanted error by the NF when the ref folder already exists
    output:
        path 'ref/ARS-UCD1.2*' optional true
        path 'ref/ref_cov'  optional true
    
    script:
    """
    mkdir -p ref

    if [ -s /home/ref/ARS-UCD1.2.fa ];then
        echo "Reference exists ..."
    else
        #echo "Downloading Bos_taurus.ARS-UCD1.2 from Ensembl v103..."
        #wget http://ftp.ensembl.org/pub/release-103/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
        aria2c -x 16 https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa.gz
    fi;

    if [ -s /home/ref/ARS-UCD1.2.1.bt2 ];then
            echo "Reference exists and indices are in the right folder."
    else
        #pigz -d -p $params.all_threads Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
        pigz -d -p $params.all_threads ARS-UCD1.2_Btau5.0.1Y.fa.gz
        mv ARS-UCD1.2_Btau5.0.1Y.fa ref/ARS-UCD1.2.fa
        samtools faidx ref/ARS-UCD1.2.fa
        echo "Indexing Bos_taurus.ARS-UCD1.2 for the bowtie2..."
        bowtie2-build --threads $params.all_threads ref/ARS-UCD1.2.fa ref/ARS-UCD1.2
        #the 1000bull genome MT is longer than ENSEMBL and this file should be reproduced for the 1KB runs
        awk '{print \$1,\$2+2}' ref/ARS-UCD1.2.fa.fai > ref/ref_cov
    fi;
    """
}



process mapper {
    tag "Mapping using bowtie2..."
    publishDir params.out , mode: 'copy', overWrite: true
    cpus params.all_threads
    maxForks 100
    cache true
    containerOptions "-v $projectDir/ref:/home/ref:ro"    
    
    input:
        tuple val(name) , path(trimmed_fastq)
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
        mkdir -p bams && \
        bowtie2 -p $params.all_threads --met-file ${name}.metrics --very-sensitive \
        --rg-id ${name} --rg LB:${name} --rg PL:ILLUMINA --rg SM:${name} \
        -x /home/ref/ARS-UCD1.2 \
        -U ${trimmed_fastq} | \
        samtools view -@ $params.all_threads -bS -F 4 | \
        samtools sort -@ $params.all_threads -o ${name}.bam
        samtools index ${name}.bam

    """
       
}
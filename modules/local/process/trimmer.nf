#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

process TRIMMER {
    tag "Trimming by TagDust2..."
    label 'medium'
    
    publishDir "$projectDir/trimmed" , mode: 'copy', overWrite: true
    maxForks 100
    cache true
    debug false
    
    input:
        tuple val(sample), val(barcode), path(directory)
    output:
        path( '*_BC_*.fq.gz' ), emit: OUT_trimmed
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 
    script:

    // For the docker runs the tools inside the scripts local are not in the right Path due to $projectDir which is not defined in the docker PATH. Needs costumising the docker PATH to explicitly declare the modules path
    """
        /modules/local/scripts/tagdust_2.33/src/tagdust ${directory} \
        -t $task.cpus \
        -1 B:${barcode} \
        -2 F:CAGNNNG \
        -3 R:N \
        -4 P:ATCTCGTATGCCGTCTTCTGCTT \
        -dust 100  \
        -o ${sample}
        
        pigz --force -p $task.cpus ${sample}_BC_${barcode}.fq
        rm -f ${sample}_un.fq
    """
       
}



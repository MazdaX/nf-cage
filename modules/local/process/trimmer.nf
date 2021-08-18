#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

def trimKeeper() {
  cmd="""mkdir -p $projectDir/trimmed"""
  result=cmd.execute().text
}

params.out="$projectDir/trimmed"

process trimmer {
    tag "Trimming by TagDust2..."
    label "proccess_wsl"
    publishDir params.out , mode: 'copy', overWrite: true
    cpus 6
    maxForks 100
    cache true
    
    input:
        tuple val(sample), val(barcode), path(directory)
    output:
        path 'trimmed/*_BC_*.fq.gz' , emit: OUT_trimmed
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 
    script:
    """
        $projectDir/modules/local/scripts/tagdust_2.33/src/tagdust ${directory} \
        -t 6 \
        -1 B:${barcode} \
        -2 F:CAGNNNG \
        -3 R:N \
        -4 P:ATCTCGTATGCCGTCTTCTGCTT \
        -dust 100  \
        -o $projectDir/trimmed/${sample}
        
        pigz --force -p 6 $projectDir/trimmed/${sample}_BC_${barcode}.fq
        rm -f $projectDir/trimmed/${sample}_un.fq
    """
       
}



#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

def demuxMaker() {
  cmd="""mkdir -p $projectDir/demux"""
  result=cmd.execute().text
}

//Should be coded for ARG in the future
params.allowed_mismatch = 1
//params.in="$projectDir/fastq_files/*.fastq.gz"
params.out="$projectDir/demux"
//params.out2="$projectDir/demux/CAGE*/*.fastq"

process demux {
    tag "Demultiplexing..."
    label "proccess_wsl"
    publishDir params.out , mode: 'copy', overWrite: true
    cpus 6
    maxForks 100
    cache true
    //afterScript 'echo "Done!!" > reporter.txt'
    

    input:
        path IN_fastq               //the name of path variable should differ from workflow variable declarations
        path IN_barcodes
        params.allowed_mismatch

    output:
        path '**/*.fastq' , emit: OUT_demux
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 
    script:
    """
        NAME=\$(basename -s .fastq.gz $IN_fastq)
        mkdir \${NAME}
        zcat $IN_fastq |\
        $projectDir/modules/local/scripts/fastx_0.0.13/bin/fastx_barcode_splitter.pl \
        --bcfile $IN_barcodes \
        --bol --mismatch $params.allowed_mismatch \
        --prefix \${NAME}/ \
        --suffix ".fastq" 
    """
       
}


process merger {
    tag "Gathering demux..."
    label "proccess_wsl"
    publishDir params.out , mode: 'copy' , overWrite: true
    cpus 6
    maxForks 100
    cache true
    
    input:
        tuple val(sample), path(directory)
    output:
        path '**.f*q.gz' , emit: OUT_merger
    script:
    """
    for d in ${directory};do 
        cat \${d}/${sample}.fastq >> ${sample}.fastq
    done;
    pigz --force -p 6 ${sample}.fastq
            
    """
}
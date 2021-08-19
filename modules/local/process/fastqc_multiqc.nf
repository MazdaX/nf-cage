#!/usr/bin/env nextflow

//demux, merge, fastqc and multiqc of Cattle CAGE datasets
// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

def reportsMaker() {
  cmd="""mkdir -p $projectDir/reports"""
  result=cmd.execute().text
}


//params.in still doesn't work as intended
//params.in="$projectDir/fastq_files/*.fastq.gz"
params.out="$projectDir/reports"


process fastqc {
    tag "FastQC..."
    label "process_wsl"
    publishDir params.out , mode:'copy', overwrite: true
    cpus 6
    maxForks 100
    cache true

    input:
      path IN_fastq
    output:
      path "*.zip" , emit: fastqc_zip_file
      path "*.html", emit: fastqc_html_file
   
    script:
      """
      fastqc --no-extract $IN_fastq 
      """    
}


process multiqc {
    tag "MultiQC..."
    label "process_wsl"
    //Publish is the final output directory
    publishDir "$projectDir/reports", mode:'copy'
    cpus 1
    maxForks 1
    cache true

    input:
      path OUT_fastqc    
    output:
      path "multiqc_*.html"
      path "multiqc_data"
    script:
      """
      multiqc  $OUT_fastqc
      """
}



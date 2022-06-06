#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

def tmpMaker() {
  cmd="""mkdir -p $projectDir/tmp"""
  result=cmd.execute().text
}


process FASTQC {
    tag "Running FastQC..."
    label 'small'

    publishDir "$projectDir/reports/${sample_id}" , mode:'copy', overwrite: true
    cache true
    debug false
    
    input:
      tuple val(sample_id), path(fastq_ch)
    output:
      path( "*.zip" ), emit: fastqc_zips
      path( "*.html" ), emit: fastqc_htmls
   
    script:
      """
      fastqc --no-extract $fastq_ch 
      """    
}


process MULTIQC {
    tag "Running MultiQC..."
    label 'small'

    //Publish is the final output directory
    publishDir "$projectDir/reports", mode:'copy'
    cache true
    debug false

    input:
      path fastqc_zips

    output:
      path( "multiqc_*.html" )
      path( "multiqc_data" )

    script:
      """
      multiqc  $fastqc_zips
      """
}



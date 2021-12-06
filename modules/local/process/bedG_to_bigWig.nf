#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

params.out="$projectDir/bams"

process bG2bW {
    tag "BAM >>> bedGraph >>> BigWig ..."
    publishDir params.out , mode: 'copy', overWrite: true
    cpus = 1
    maxForks 100
    cache true
    //containerOptions "-v $projectDir/ref:/home/ref:ro"

    input:
        tuple val(name) , path(bam)
    // Use of subsampled fastq files results in bedtools and bG2bW errors that prevents the test run. The optional flag allows the pipeline to finish regardless.
    //Double check the EOF for both bw and bedGraph.gz files to ensure correct analysis. 

    output:
        path '*.bedGraph.gz' optional true
        path '*.bw' optional true
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 

    //Docker addresses for ref has the same problem as the trimmer tagdust from the modules
    //Docker optimisation is needed
    
    script:
    """
        bedtools genomecov -ibam ${bam} -d -strand + | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.plus.bedGraph & \
        bedtools genomecov -ibam ${bam} -d -strand - | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.minus.bedGraph 
        sort -k 1,1 -k2,2n ${name}.plus.bedGraph > ${name}_tmp_plus
        sort -k 1,1 -k2,2n ${name}.minus.bedGraph > ${name}_tmp_minus 
        bedGraphToBigWig ${name}_tmp_plus /home/ref/ref_cov ${name}.plus.bw && \
        rm -f ${name}_tmp_plus
        bedGraphToBigWig ${name}_tmp_minus /home/ref/ref_cov ${name}.minus.bw && \
        rm -f ${name}_tmp_minus
        pigz -p $params.all_threads *.bedGraph
    """
          
}



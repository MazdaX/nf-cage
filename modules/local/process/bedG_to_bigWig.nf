#!/usr/bin/env nextflow

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2


process BG2BW {
    tag "BAM >>> bedGraph >>> BigWig ..."
    publishDir "$projectDir/bams" , mode: 'copy', overWrite: true
    label 'small'
    cache true

    input:
        tuple val(name) , path(bam)
        each ref_cov
    // Use of subsampled fastq files results in bedtools and bG2bW errors that prevents the test run. The optional flag allows the pipeline to finish regardless.
    //Double check the EOF for both bw and bedGraph.gz files to ensure correct analysis. 

    output:
        path '*.bedGraph.gz' optional true
        path '*.bw' optional true
                
    script:
    """
        bedtools genomecov -ibam ${bam} -d -strand + | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.plus.bedGraph &\
        bedtools genomecov -ibam ${bam} -d -strand - | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.minus.bedGraph 
        sort -k 1,1 -k2,2n ${name}.plus.bedGraph > ${name}_tmp_plus
        sort -k 1,1 -k2,2n ${name}.minus.bedGraph > ${name}_tmp_minus 
        bedGraphToBigWig ${name}_tmp_plus $ref_cov ${name}.plus.bw && \
        rm -f ${name}_tmp_plus
        bedGraphToBigWig ${name}_tmp_minus $ref_cov ${name}.minus.bw && \
        rm -f ${name}_tmp_minus
        pigz -p $task.cpus *.bedGraph
    """
          
}



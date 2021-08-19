#!/usr/bin/env nextflow

// The $baseDir is an env variable in DSL2

//Enabling the DSL2 syntax 
nextflow.enable.dsl=2

params.out="$projectDir/bams"

process bG2bW {
    tag "BAM >>> bedGraph >>> BigWig ..."
    label "proccess_wsl"
    publishDir params.out , mode: 'copy', overWrite: true
    cpus 1
    maxForks 1
    cache true
    
    input:
        tuple val(name) , path(bam)
    output:
        path '*.bedGraph.gz' , emit: OUT_bedGraph
        path '*.bw', emit: OUT_bigWig
                
    //In order to use system \$vars as well as DSL $vars
    // The issue is the for loop behaviour which cannot be invoked (process)
    // more than once in a workflow. Need to solve this issue. 

    //Docker addresses fir ref has the same problem as the trimmer tagdust from the modules
    //Docker optimisation is needed
    
    script:
    """
        bedtools genomecov -ibam ${bam} -d -strand + | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.plus.bedGraph & \
        bedtools genomecov -ibam ${bam} -d -strand - | awk -v width=1 '!(\$1~/^NW/)&&(\$3!=0) {print \$1,\$2,\$2+width,\$3}' > ${name}.minus.bedGraph 
        sort -k 1,1 -k2,2n ${name}.plus.bedGraph > ${name}_tmp_plus
        sort -k 1,1 -k2,2n ${name}.minus.bedGraph > ${name}_tmp_mins 
        /modules/local/scripts/UCSC/bedGraphToBigWig ${name}_tmp_plus /ref/ref_cov ${name}.plus.bw && \
        rm -f ${name}_tmp_plus
        /modules/local/scripts/UCSC/bedGraphToBigWig ${name}_tmp_minus /ref/ref_cov ${name}.minus.bw && \
        rm -f ${name}_tmp_minus
        pigz *.bedGraph
    """
          
}


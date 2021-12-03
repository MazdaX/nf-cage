nextflow.enable.dsl = 2
include {bG2bW} from './modules/local/process/bedG_to_bigWig.nf'

workflow {
    input_bams = channel.fromPath("$projectDir/bams/*.bam")
    convert_single_channel = input_bams
                                .flatten()
                                .map { it->tuple( it
                                                    .simpleName
                                                    .split('_')[0]                                                    
                                                    ,it) }
                                
    convert_single_channel.view()
    //bG2bW(convert_single_channel)

}
#!/usr/bin/env nextflow

/*
========================================================================================
                                    nf-CAGE
========================================================================================
nf-CAGE analysis pipeline 
Developed and maintained by Dr Mazdak Salavati
#### Homepage / Documentation
https://github.com/MazdaX/nf-cage.git
docker pull mazdax/nf-cage:latest
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
params.help = ''

if (params.help) {
    include {helpMessage} from './modules/local/process/help.nf'
    helpMessage()
    exit 0
}

////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

log.info """\

===============================================
        nf-cage BovReg's pipeline
        github.com/mazdax/nf-cage
        docker pull mazdax/nf-cage
===============================================
Input : Raw Illumina CAGE sequences
Input : Barcode list TSV (i.e. sample\tbarcode)
Input : Reference Genome FASTA (URL or local)
Output : Strand specific bp resolution bigWig 
Running task: $params.all_in_one
Workflow version: ${workflow.manifest.version}
-----------------------------------------------
"""

//include {module as module_alias; module as module_alias2}
include {DEMUX;MERGER} from './modules/local/process/demux.nf'
include {TRIMMER} from './modules/local/process/trimmer.nf'
include {tmpMaker;FASTQC as qc_pre;FASTQC as qc_post;MULTIQC} from './modules/local/process/fastqc_multiqc.nf'
include {DOWNLOADREF; BT2BUILD; BT2MAPPER} from './modules/local/process/mapper.nf'
include {BG2BW} from './modules/local/process/bedG_to_bigWig.nf'

params.ref_fasta = 'https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa.gz'
params.raw_fastq = "$projectDir/fastq_files/*.gz"
params.barcodes = Channel.value("$projectDir/barcode_files/barcode.list")


workflow {
        //Queue channel allows for parallel execution of all files in a FIFO queue
        raw_fastq_ch = Channel
                .fromPath(params.raw_fastq)
                .map({it -> [it.simpleName,it]})    
        
        //barcodes should be read unlimited times >>> value.channel 
        //barcodes = Channel.value("$projectDir/barcode_files/barcode.list")

        
        //AIO workflow
	tmpMaker()
        qc_pre(raw_fastq_ch) 

        
        
        DEMUX(raw_fastq_ch,params.barcodes)
        demux_single_ch=DEMUX.out.OUT_demux
                                .flatten()
                                .map { it-> [it.simpleName.replaceAll('CAGE_[0-9][0-9]_',''),it] }
                                .groupTuple(by: 0 , size: -1)
                                                                    
        MERGER(demux_single_ch)
        
        
        //Need to emit the merger output so it can be pegged to the qc_post
        fastqc_single_ch=MERGER.out.OUT_merger
                                .flatten()
                                .map({it -> [it.simpleName,it]})   

        // This post QC should get the merged files not individuals
        qc_post(fastqc_single_ch)
        
        demux_single_channel=qc_pre.out.fastqc_zips
                                .mix(qc_post.out.fastqc_zips)
                                .flatten()
                        
        MULTIQC(demux_single_channel.collect())
        
        trimmer_single_ch=MERGER.out.OUT_merger
                                 .flatten()
                                 .map { it->tuple( it.simpleName,it ) }
                                 .join(channel
                                        .fromPath("$projectDir/barcode_files/barcode.list")
                                        .splitCsv(header:false, sep: ' '))
                                 .map { it->tuple(it[0],it[2],it[1])}
                                 
        
        TRIMMER(trimmer_single_ch)

        
        //Mapping and BedGraph to BigWig conversion
        mapping_single_ch=TRIMMER.out.OUT_trimmed
                                .flatten()
                                .map{ it-> tuple(it.simpleName,it)}
        
        
        // Downloading reference fasta file
        DOWNLOADREF(Channel.of(params.ref_fasta))
        
        // Build bowtie2 index
        BT2BUILD(DOWNLOADREF.out)

        mapping_index_ch=BT2BUILD.out.bowtie_index
                                .flatten()
                                .map({ it -> "$it".replaceAll(/.[1-4].bt2|.rev.[1-2].bt2/,"") })
                                .unique()
                                .map({ it -> new File(it) })
                                //.view()

        // Mappping reads
        BT2MAPPER(mapping_single_ch,mapping_index_ch)
        
        convert_single_channel=BT2MAPPER.out.OUT_mapped
                                .flatten()
                                .map{ it-> tuple(it.simpleName,it)}
        
        BG2BW(convert_single_channel,BT2BUILD.out.ref_cov)
        
}


////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

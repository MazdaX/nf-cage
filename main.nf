#!/usr/bin/env nextflow
/*
========================================================================================
                                    nf-CAGE
========================================================================================
 nf-CAGE analysis pipline 
 Developed and maintained by Dr Mazdak Salavati
 #### Homepage / Documentation
 https://git.ecdf.ed.ac.uk/ruminant_functional_genomics/nextflow.git
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2
// for wsl runs

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
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



include {demuxMaker;demux;merger} from './modules/local/process/demux.nf'
include {trimKeeper;trimmer} from './modules/local/process/trimmer.nf'
include {reportsMaker;fastqc as qc_pre;fastqc as qc_post;multiqc} from './modules/local/process/fastqc_multiqc.nf'
//include {module as module_alias; module as module_alias2}

workflow {
        //path doesn't take relative paths but file does
        //Queue channel allows for parallel executation of all files in a FIFO queue
        raw_fastq=channel.fromPath("$projectDir/fastq_files/*.fastq.gz")     
        
        //barcodes should be read unlimited times >>> value.channel  >>> channel.value(1)
        //a proccess is needed to read in the list of barcodes to a value list where they could be iterated
        barcodes=channel.value("$projectDir/barcode_files/barcode.list")

        demux_fastq=channel.fromPath("$projectDir/demux/*.fastq.gz")
        
        //The output from demultiplexing to be merged in to single sample fastq
        trimmed_files=channel.fromPath("$projectDir/trimmed/CAGE*_BC_*.fq.gz")
        
        
        //AIO workflow
        reportsMaker()
              
        qc_pre(raw_fastq) 
        // emit index 1 is the html files
        // emit index 0 is the zip files
        
        
        demux(raw_fastq,barcodes)
        //demux.out.OUT_demux.view()
        
        demux_single_channel=demux.out.OUT_demux
                                .flatten()
                                .map { it->tuple(it.simpleName,it.parent) }
                                .groupTuple(by:0)
        
        //demux_single_channel.view()
        
        merger(demux_single_channel)
        
        //Need to emit the merger output so it can be pegged to the qc_post
        fastqc_single_channel=merger.out.OUT_merger
                                .flatten()
        // This post QC should get the merged files not individuals
        qc_post(fastqc_single_channel)
        
        demux_single_channel=qc_pre.out.fastqc_zip_file
                                .mix(qc_post.out.fastqc_zip_file)
                                .flatten()
                        
        multiqc(demux_single_channel.collect())
        

        trimmer_single_channel=demux_fastq
                                 .flatten()
                                 .map { it->tuple( it.simpleName,it ) }
                                 .join(channel
                                        .fromPath("$projectDir/barcode_files/barcode.list")
                                        .splitCsv(header:false, sep: ' '))
                                 .map { it->tuple(it[0],it[2],it[1])}
        
        trimmer_single_channel.view()
        trimKeeper()
        //trimmer(trimmer_single_channel)

        // second call to the same proccess with different input files
        //fastqc(trimmed_files)
        //multiqc(fastqc.out[0].collect())

}


////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
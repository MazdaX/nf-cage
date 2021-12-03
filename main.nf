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
Input : Reference Genome FASTA (bowtie2 index)
Output : Strand specific bp resolution bigWig 
Running task: $params.all_in_one
-----------------------------------------------
"""

//include {module as module_alias; module as module_alias2}
include {demuxMaker;demux;merger} from './modules/local/process/demux.nf'
include {trimKeeper;trimmer} from './modules/local/process/trimmer.nf'
include {tmpMaker;reportsMaker;fastqc as qc_pre;fastqc as qc_post;multiqc} from './modules/local/process/fastqc_multiqc.nf'
include {downloadRef; bowtie2Build; mapper} from './modules/local/process/mapper.nf'
include {bG2bW} from './modules/local/process/bedG_to_bigWig.nf'

workflow {
        //path doesn't take relative paths but file does
        //Queue channel allows for parallel execution of all files in a FIFO queue
        raw_fastq=channel.fromPath("$projectDir/fastq_files/*.gz")     
        
        //barcodes should be read unlimited times >>> value.channel  >>> channel.value(1)
        //a process is needed to read in the list of barcodes to a value list where they could be iterated
        barcodes=channel.value("$projectDir/barcode_files/barcode.list")

        // reference fasta file
        fasta=channel.value('https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa.gz')
	//fasta=channel.value('http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
 
        //AIO workflow
        reportsMaker()
	tmpMaker()
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
        
        trimmer_single_channel=merger.out.OUT_merger
                                 .flatten()
                                 .map { it->tuple( it.simpleName,it ) }
                                 .join(channel
                                        .fromPath("$projectDir/barcode_files/barcode.list")
                                        .splitCsv(header:false, sep: ' '))
                                 .map { it->tuple(it[0],it[2],it[1])}
                                 
        
        //trimmer_single_channel.view()
        trimKeeper()
        trimmer(trimmer_single_channel)

        
        //Mapping and BedGraph to BigWig conversion
        mapping_single_channel=trimmer.out.OUT_trimmed
                                .flatten()
                                .map{ it-> tuple(it.simpleName,it)}
        
        // Downloading reference fasta file
        downloadRef(fasta)
        // Build bowtie2 index
        bowtie2Build(downloadRef.out.fasta)
        // Mappping reads
        mapper(mapping_single_channel, bowtie2Build.out.bowtie_index)

        convert_single_channel=mapper.out.OUT_mapped
                                .flatten()
                                .map{ it-> tuple(it.simpleName,it)}
        convert_single_channel.view()

        bG2bW(convert_single_channel)
        
}


////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

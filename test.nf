nextflow.enable.dsl = 2

//At the end overwrite the main.nf
include {DEMUX;MERGER} from './modules/local/process/demux.nf'
include {TRIMMER} from './modules/local/process/trimmer.nf'
include {FASTQC as qc_pre;FASTQC as qc_post;MULTIQC} from './modules/local/process/fastqc_multiqc.nf'
include {DOWNLOADREF; bowtie2Build; mapper} from './modules/local/process/mapper.nf'
include {bG2bW} from './modules/local/process/bedG_to_bigWig.nf'

raw_fastq_ch = Channel
                .fromPath("$projectDir/fastq_files/*.gz")
                .map({it -> [it.simpleName,it]})
//                .view()

barcodes = Channel.value("$projectDir/barcode_files/barcode.list")
/*
workflow{
    //qc_pre(raw_fastq_ch)
    DEMUX(raw_fastq_ch,barcodes)
    demux_single_ch=DEMUX.out.OUT_demux
                                .flatten()
                                .map { it-> [it.simpleName.replaceAll('CAGE_[0-9][0-9]_',''),it] }
                                .groupTuple(by: 0 , size: -1)
                                                                    
    MERGER(demux_single_ch)
    //multiqc(preqc.out.fastqc_zips.collect())
}
*/


params.ref_fasta = 'https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa.gz'
Channel
    .of( ["$projectDir/ref/",params.ref_fasta] )
    .map( {folder, url -> [folder+url.tokenize('/')[-1]] })
    .view()

include { INDEX } from './index'
include { QUANT } from './quant'
include { FASTQC } from './fastqc'
include { FASTP } from './nf-core/fastp'

workflow RNASEQ {
    take:
    read_pairs_ch
    transcriptome

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)

    // Adapt the flat tuple [id, fastq_1, fastq_2] to the nf-core meta-map
    // format expected by FASTP: [meta, [reads], adapter_fasta]
    fastp_input_ch = read_pairs_ch.map { id, fastq_1, fastq_2 ->
        [[id: id, single_end: false], [fastq_1, fastq_2], []]
    }

    // Run fastp for adapter/quality trimming
    FASTP(fastp_input_ch, false, false, false)

    // Convert FASTP trimmed reads back to flat tuple [id, R1, R2] for QUANT
    trimmed_ch = FASTP.out.reads.map { meta, reads ->
        [meta.id, reads[0], reads[1]]
    }

    // Quantify using the trimmed reads
    quant_ch = QUANT(trimmed_ch, index)
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    samples = samples_ch
}

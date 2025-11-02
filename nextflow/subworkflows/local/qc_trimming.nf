/*
========================================================================================
    QC and Trimming Subworkflow
========================================================================================
    Performs quality control and adapter trimming on raw FASTQ files
----------------------------------------------------------------------------------------
*/

include { FASTQC as FASTQC_RAW      } from '../../modules/nf-core/fastqc/main'
include { FASTP                     } from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_TRIMMED  } from '../../modules/nf-core/fastqc/main'

workflow QC_TRIMMING {
    take:
    ch_reads        // channel: [ val(meta), [ reads ] ]
    adapter_fasta   // path: adapter sequences (optional)
    skip_fastqc     // val: skip FastQC
    skip_trimming   // val: skip trimming

    main:
    ch_versions = Channel.empty()
    ch_fastqc_raw_html = Channel.empty()
    ch_fastqc_raw_zip = Channel.empty()
    ch_fastqc_trimmed_html = Channel.empty()
    ch_fastqc_trimmed_zip = Channel.empty()
    ch_fastp_json = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_trim_reads = ch_reads

    //
    // MODULE: Run FastQC on raw reads
    //
    if (!skip_fastqc) {
        FASTQC_RAW (
            ch_reads
        )
        ch_fastqc_raw_html = FASTQC_RAW.out.html
        ch_fastqc_raw_zip  = FASTQC_RAW.out.zip
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    //
    // MODULE: Run fastp for adapter trimming and quality filtering
    //
    if (!skip_trimming) {
        FASTP (
            ch_reads,
            adapter_fasta,
            false,  // save_trimmed_fail
            false   // save_merged
        )
        ch_trim_reads = FASTP.out.reads
        ch_fastp_json = FASTP.out.json
        ch_fastp_html = FASTP.out.html
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        //
        // MODULE: Run FastQC on trimmed reads
        //
        if (!skip_fastqc) {
            FASTQC_TRIMMED (
                ch_trim_reads
            )
            ch_fastqc_trimmed_html = FASTQC_TRIMMED.out.html
            ch_fastqc_trimmed_zip  = FASTQC_TRIMMED.out.zip
            ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())
        }
    }

    emit:
    reads = ch_trim_reads                           // channel: [ val(meta), [ reads ] ]

    fastqc_raw_html = ch_fastqc_raw_html           // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip  = ch_fastqc_raw_zip            // channel: [ val(meta), [ zip ] ]

    fastp_html = ch_fastp_html                      // channel: [ val(meta), [ html ] ]
    fastp_json = ch_fastp_json                      // channel: [ val(meta), [ json ] ]

    fastqc_trimmed_html = ch_fastqc_trimmed_html   // channel: [ val(meta), [ html ] ]
    fastqc_trimmed_zip  = ch_fastqc_trimmed_zip    // channel: [ val(meta), [ zip ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}

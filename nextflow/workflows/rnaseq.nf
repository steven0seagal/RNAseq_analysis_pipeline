/*
========================================================================================
    Standard RNA-seq Differential Expression Workflow
========================================================================================
    Complete workflow for RNA-seq analysis from FASTQ to differential expression

    Steps:
    1. Quality control (FastQC)
    2. Adapter trimming (fastp)
    3. Genome alignment (STAR)
    4. Gene quantification (featureCounts)
    5. QC aggregation (MultiQC)
    6. Differential expression (DESeq2)
----------------------------------------------------------------------------------------
*/

include { QC_TRIMMING               } from '../subworkflows/local/qc_trimming'
include { ALIGNMENT_QUANTIFICATION  } from '../subworkflows/local/alignment_quantification'
include { DESEQ2_ANALYSIS           } from '../modules/local/deseq2/analysis'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'

workflow RNASEQ {

    take:
    ch_input       // channel: [ val(meta), [ reads ] ]
    ch_fasta       // channel: [ val(meta), path(fasta) ]
    ch_gtf         // channel: [ val(meta), path(gtf) ]
    ch_star_index  // channel: path(star_index) or empty
    ch_metadata    // path: metadata file

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Quality control and trimming
    //
    QC_TRIMMING (
        ch_input,
        [],  // adapter_fasta (empty for auto-detection)
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(QC_TRIMMING.out.versions)

    // Collect QC files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(
        QC_TRIMMING.out.fastqc_raw_zip.map { meta, zip -> zip }.ifEmpty([])
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        QC_TRIMMING.out.fastp_json.map { meta, json -> json }.ifEmpty([])
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        QC_TRIMMING.out.fastqc_trimmed_zip.map { meta, zip -> zip }.ifEmpty([])
    )

    //
    // SUBWORKFLOW: Alignment and quantification
    //
    ALIGNMENT_QUANTIFICATION (
        QC_TRIMMING.out.reads,
        ch_fasta,
        ch_gtf,
        ch_star_index,
        false,  // star_ignore_sjdbgtf
        params.seq_platform ?: '',
        params.seq_center ?: ''
    )
    ch_versions = ch_versions.mix(ALIGNMENT_QUANTIFICATION.out.versions)

    // Collect alignment files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(
        ALIGNMENT_QUANTIFICATION.out.star_log.map { meta, log -> log }
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        ALIGNMENT_QUANTIFICATION.out.counts_summary.map { meta, summary -> summary }
    )

    //
    // MODULE: MultiQC - Aggregate QC reports
    //
    if (!params.skip_multiqc) {
        MULTIQC (
            ch_multiqc_files.collect(),
            [],  // multiqc_config
            [],  // extra_multiqc_config
            []   // multiqc_logo
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    //
    // MODULE: DESeq2 differential expression analysis
    //
    ch_deseq2_results = Channel.empty()
    ch_deseq2_plots = Channel.empty()

    if (!params.skip_deseq2 && ch_metadata) {
        DESEQ2_ANALYSIS (
            ALIGNMENT_QUANTIFICATION.out.counts.map { meta, counts -> counts },
            ch_metadata,
            params.control_condition,
            params.treatment_condition,
            params.padj_threshold,
            params.log2fc_threshold
        )
        ch_deseq2_results = DESEQ2_ANALYSIS.out.results
        ch_deseq2_plots = Channel.empty().mix(
            DESEQ2_ANALYSIS.out.pca_plot,
            DESEQ2_ANALYSIS.out.volcano_plot,
            DESEQ2_ANALYSIS.out.heatmap
        )
        ch_versions = ch_versions.mix(DESEQ2_ANALYSIS.out.versions)
    }

    emit:
    // QC outputs
    fastqc_raw_html    = QC_TRIMMING.out.fastqc_raw_html
    fastp_html         = QC_TRIMMING.out.fastp_html
    fastqc_trimmed_html = QC_TRIMMING.out.fastqc_trimmed_html

    // Alignment outputs
    bam                = ALIGNMENT_QUANTIFICATION.out.bam
    bai                = ALIGNMENT_QUANTIFICATION.out.bai
    star_log           = ALIGNMENT_QUANTIFICATION.out.star_log

    // Quantification outputs
    counts             = ALIGNMENT_QUANTIFICATION.out.counts
    counts_summary     = ALIGNMENT_QUANTIFICATION.out.counts_summary

    // MultiQC output
    multiqc_report     = params.skip_multiqc ? Channel.empty() : MULTIQC.out.report

    // DESeq2 outputs
    deseq2_results     = ch_deseq2_results
    deseq2_plots       = ch_deseq2_plots

    // Version information
    versions           = ch_versions
}

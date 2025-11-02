/*
========================================================================================
    Alignment and Quantification Subworkflow
========================================================================================
    Performs genome alignment with STAR and gene quantification with featureCounts
----------------------------------------------------------------------------------------
*/

include { STAR_GENOMEGENERATE       } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                } from '../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { SUBREAD_FEATURECOUNTS     } from '../../modules/nf-core/subread/featurecounts/main'

workflow ALIGNMENT_QUANTIFICATION {
    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_gtf              // channel: [ val(meta), path(gtf) ]
    ch_star_index       // channel: path(star_index) or empty
    star_ignore_sjdbgtf // val: ignore GTF for STAR
    seq_platform        // val: sequencing platform
    seq_center          // val: sequencing center

    main:
    ch_versions = Channel.empty()
    ch_genome_bam = Channel.empty()
    ch_genome_bai = Channel.empty()
    ch_star_index_out = ch_star_index

    //
    // MODULE: Generate STAR index if not provided
    //
    if (!ch_star_index) {
        STAR_GENOMEGENERATE (
            ch_fasta,
            ch_gtf
        )
        ch_star_index_out = STAR_GENOMEGENERATE.out.index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // MODULE: Align reads to reference genome with STAR
    //
    STAR_ALIGN (
        ch_reads,
        ch_star_index_out,
        ch_gtf,
        star_ignore_sjdbgtf,
        seq_platform,
        seq_center
    )
    ch_genome_bam = STAR_ALIGN.out.bam_sorted
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // MODULE: Index BAM files
    //
    SAMTOOLS_INDEX (
        ch_genome_bam
    )
    ch_genome_bai = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Gene-level quantification with featureCounts
    //
    // Prepare input for featureCounts - collect all BAM files
    ch_featurecounts_input = ch_genome_bam
        .map { meta, bam -> bam }
        .collect()
        .map { bams ->
            def meta = [id: 'all_samples']
            [ meta, bams ]
        }

    SUBREAD_FEATURECOUNTS (
        ch_featurecounts_input,
        ch_gtf.map { meta, gtf -> gtf }
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    emit:
    bam           = ch_genome_bam                      // channel: [ val(meta), path(bam) ]
    bai           = ch_genome_bai                      // channel: [ val(meta), path(bai) ]

    star_log      = STAR_ALIGN.out.log_final           // channel: [ val(meta), path(log) ]

    counts        = SUBREAD_FEATURECOUNTS.out.counts   // channel: [ val(meta), path(counts) ]
    counts_summary = SUBREAD_FEATURECOUNTS.out.summary // channel: [ val(meta), path(summary) ]

    star_index    = ch_star_index_out                  // channel: path(star_index)

    versions      = ch_versions                        // channel: [ versions.yml ]
}

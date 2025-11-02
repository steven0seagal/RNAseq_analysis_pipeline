#!/usr/bin/env nextflow
/*
========================================================================================
    RNAseq Analysis Pipeline - Nextflow Implementation
========================================================================================
    Main workflow entry point

    Supports:
    - Standard RNA-seq differential expression analysis
    - Circular RNA (circRNA) discovery and quantification
    - Multimodal cancer genomics analysis

    Usage:
        nextflow run main.nf --input samplesheet.csv --mode rnaseq
        nextflow run main.nf --input samplesheet.csv --mode circrna
        nextflow run main.nf --input samplesheet.csv --mode multimodal
        nextflow run main.nf --input samplesheet.csv --mode integrated
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
========================================================================================
*/

include { RNASEQ       } from './workflows/rnaseq'
// include { CIRCRNA      } from './workflows/circrna'
// include { MULTIMODAL   } from './workflows/multimodal'
// include { INTEGRATED   } from './workflows/integrated'

include { WorkflowRnaseq } from './lib/WorkflowRnaseq'

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// Print help message
def helpMessage() {
    log.info"""
    ============================================================================
     RNAseq Analysis Pipeline v${workflow.manifest.version}
    ============================================================================

    Usage:
        nextflow run main.nf --input samplesheet.csv [options]

    Mandatory arguments:
        --input [file]              Path to comma-separated sample sheet
        --mode [str]                Analysis mode: 'rnaseq', 'circrna', 'multimodal', or 'integrated'
        --outdir [path]             The output directory where results will be saved
        --genome_fasta [file]       Path to genome FASTA file
        --annotation_gtf [file]     Path to gene annotation GTF file

    Optional arguments:
        --metadata [file]           Path to metadata file for DE analysis (default: config/metadata.tsv)
        --strand_specificity [int]  Strandness: 0=unstranded, 1=forward, 2=reverse (default: 2)
        --skip_fastqc               Skip FastQC
        --skip_trimming             Skip adapter trimming
        --skip_multiqc              Skip MultiQC report generation
        --skip_deseq2               Skip DESeq2 differential expression analysis

    Profiles:
        -profile conda              Use conda for dependency management
        -profile docker             Use docker containers
        -profile singularity        Use singularity containers
        -profile test               Run with test data

    For more information, see the documentation at:
        https://github.com/username/RNAseq_analysis_pipeline
    ============================================================================
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate input parameters
WorkflowRnaseq.validateInputParams(params)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Print parameter summary
    log.info WorkflowRnaseq.paramsSummaryLog(workflow, params)

    //
    // Create input channels
    //
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = row.fastq_2 ? false : true
            meta.condition = row.condition ?: 'unknown'

            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null

            if (fastq_2) {
                return [meta, [fastq_1, fastq_2]]
            } else {
                return [meta, [fastq_1]]
            }
        }

    // Create reference channels
    ch_fasta = params.genome_fasta
        ? Channel.fromPath(params.genome_fasta, checkIfExists: true)
            .map { fasta -> [[id: 'genome'], fasta] }
        : Channel.empty()

    ch_gtf = params.annotation_gtf
        ? Channel.fromPath(params.annotation_gtf, checkIfExists: true)
            .map { gtf -> [[id: 'annotation'], gtf] }
        : Channel.empty()

    ch_star_index = params.star_index
        ? Channel.fromPath(params.star_index, checkIfExists: true)
        : Channel.empty()

    ch_metadata = params.metadata
        ? file(params.metadata, checkIfExists: true)
        : null

    //
    // Route to appropriate workflow based on mode
    //
    if (params.mode == 'rnaseq') {
        log.info "Running standard RNA-seq workflow..."

        RNASEQ (
            ch_input,
            ch_fasta,
            ch_gtf,
            ch_star_index,
            ch_metadata
        )

    } else if (params.mode == 'circrna') {
        log.info "Running circRNA analysis workflow..."
        // CIRCRNA()
        log.warn "circRNA workflow not yet implemented. Please see workflows/circrna.nf"

    } else if (params.mode == 'multimodal') {
        log.info "Running multimodal analysis workflow..."
        // MULTIMODAL()
        log.warn "Multimodal workflow not yet implemented. Please see workflows/multimodal.nf"

    } else if (params.mode == 'integrated') {
        log.info "Running integrated analysis workflow (all modes)..."
        // INTEGRATED()
        log.warn "Integrated workflow not yet implemented. Please see workflows/integrated.nf"

    } else {
        error "Invalid mode: ${params.mode}. Must be 'rnaseq', 'circrna', 'multimodal', or 'integrated'"
    }
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    if (workflow.success) {
        log.info """
        ============================================================================
         Pipeline completed successfully!
        ============================================================================
         Completed at : ${workflow.complete}
         Duration     : ${workflow.duration}
         Success      : ${workflow.success}
         Work directory : ${workflow.workDir}
         Exit status  : ${workflow.exitStatus}
        ============================================================================
         Results saved to: ${params.outdir}
        ============================================================================
        """.stripIndent()
    } else {
        log.error """
        ============================================================================
         Pipeline completed with errors!
        ============================================================================
         Completed at : ${workflow.complete}
         Duration     : ${workflow.duration}
         Success      : ${workflow.success}
         Work directory : ${workflow.workDir}
         Exit status  : ${workflow.exitStatus}
         Error report : ${workflow.errorReport}
        ============================================================================
        """.stripIndent()
    }
}

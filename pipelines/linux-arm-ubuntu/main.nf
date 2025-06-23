#!/usr/bin/env nextflow

// ------------------
// CHANNEL DEFINITIONS
// ------------------

reads_ch = Channel
    .fromFilePairs(params.input, size: 1)
    .map { name, files -> [ name.replaceAll(/_R?[12].*/, ''), files[0] ] }

genome_dir_ch = Channel.fromPath(params.genome_dir, checkIfExists: true)

// ------------------
// PROCESS: Get versions
// ------------------

process GET_VERSIONS {
    publishDir "${params.outdir}/versions", mode: 'copy'

    output:
    path "versions.txt"

    script:
    """
    echo "FastQC: \$(fastqc --version)" > versions.txt
    echo "STAR: \$(STAR --version)" >> versions.txt
    """
}

// ------------------
// PROCESS: FastQC
// ------------------

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    path "*_fastqc.{html,zip}"

    script:
    """
    fastqc ${read} --threads ${task.cpus}
    """
}

// ------------------
// PROCESS: STAR align
// ------------------

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(read)
    path genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")
    path "${sample_id}_Log.final.out"

    script:
    """
    STAR \\
        --genomeDir ${genome_dir} \\
        --readFilesIn ${read} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --runThreadN ${task.cpus}
    """
}

// ------------------
// WORKFLOW
// ------------------

workflow {
    GET_VERSIONS()
    FASTQC(reads_ch)
    STAR_ALIGN(reads_ch, genome_dir_ch)
}

#!/usr/bin/env nextflow

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{html,zip}"

    script:
    """
    fastqc ${reads} --threads ${task.cpus}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")
    path "${sample_id}_Log.final.out"

    script:
    """
    STAR --genomeDir ${genome_dir} \\
         --readFilesIn ${reads} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${sample_id}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN ${task.cpus}
    """
}

workflow {
    reads_ch = Channel
        .fromFilePairs(params.input, size: 1)
        .map { id, file -> [id.replaceAll(/_R?[12].*/, ""), file[0]] }

    genome_dir_ch = Channel.fromPath(params.genome_dir, checkIfExists: true)

    FASTQC(reads_ch)
    STAR_ALIGN(reads_ch, genome_dir_ch)
}

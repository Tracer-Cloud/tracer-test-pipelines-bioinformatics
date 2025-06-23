#!/usr/bin/env nextflow

process FASTQC {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc $reads
    """
}

process STAR_ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path genome_dir

    output:
    path "${sample_id}_Aligned.sortedByCoord.out.bam"
    path "${sample_id}_Log.final.out"

    script:
    """
    STAR --genomeDir $genome_dir \
         --readFilesIn $reads \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_ \
         --outSAMtype BAM SortedByCoordinate
    """
}

workflow {
    reads_ch = Channel.fromFilePairs("${params.input}", size: 1)
        .map { id, file -> tuple(id.replaceAll(/_R[12].*/, ""), file) }

    genome_dir_ch = Channel.fromPath(params.genome_dir)

    FASTQC(reads_ch)
    STAR_ALIGN(reads_ch, genome_dir_ch)
}

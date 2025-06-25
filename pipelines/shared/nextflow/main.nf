#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import shared workflows
include { version_check } from './workflows/version-check.nf'
include { fasta_analysis } from './workflows/fasta-analysis.nf'

// Default workflow - version check
workflow {
    version_check()
}

// Alternative workflows can be called with parameters
// workflow {
//     fasta_analysis()
// }

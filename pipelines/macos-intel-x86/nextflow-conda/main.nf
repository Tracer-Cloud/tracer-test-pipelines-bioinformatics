nextflow.enable.dsl = 2
include { version_check } from '../../shared/nextflow/workflows/macos-version-check.nf'
workflow { version_check() } 

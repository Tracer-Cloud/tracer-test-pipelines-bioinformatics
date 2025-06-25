nextflow.enable.dsl = 2
include { workflow as VERSION_CHECK_WORKFLOW } from '../../shared/nextflow/workflows/version-check.nf'
workflow { VERSION_CHECK_WORKFLOW() } 
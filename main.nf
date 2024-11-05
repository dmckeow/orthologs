#!/usr/bin/env nextflow


// Include modules
include { BAM_INDEX } from './subworkflows/bam_index'
include { ORTHOFINDER } from './modules/nf-core/orthofinder/main'

workflow {
    if (params.run_bam_index) {

        BAM_INDEX(params.sample_sheet)

    }
}
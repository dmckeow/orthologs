#!/usr/bin/env nextflow


// Include modules
include { BAM_INDEX } from './subworkflows/bam_index'
include { RUN_ORTHOFINDER } from './subworkflows/initial_orthology'


workflow {
    if (params.run_bam_index) {
        BAM_INDEX(params.sample_sheet)
    }

    RUN_ORTHOFINDER(params.fasta_dir, params.prior_run)

}
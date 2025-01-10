#!/usr/bin/env nextflow

include { PREFILTER } from './subworkflows/prefilter'
include { INIT_ORTHO_ORTHOFINDER } from './subworkflows/init_ortho_orthofinder'

workflow {

    params.samplesheet = null
    params.search_params = null

    if (params.samplesheet == null) {
        error "Please provide a samplesheet with the '--samplesheet' option."
    }

    if (params.run.prefilter_hmmsearch && params.search_params == null) {
        log.warn "search_params needs '--search_params' option - running on whole proteome"
    }

    prefilter_metamap = Channel.empty()

    // Run the orthology workflow
    // Always run prefilter step because it also preprocesses fastas even if no seach step is run
    PREFILTER (
        params.samplesheet,
        params.search_params,
        params.outdir,
        params.run.prefilter_hmmsearch
    )
    prefilter_metamap = PREFILTER.out.metamap


    if (params.run.init_ortho_orthofinder) {
        INIT_ORTHO_ORTHOFINDER (
            params.samplesheet,
            params.search_params,
            params.outdir,
            params.mcl_inflation,
            prefilter_metamap
        )
    }
}
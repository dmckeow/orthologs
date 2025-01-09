#!/usr/bin/env nextflow

include { PREFILTER_PROTEOMES } from './subworkflows/prefilter_proteomes'
include { INIT_ORTHO_ORTHOFINDER } from './subworkflows/init_ortho_orthofinder'

workflow {

    params.samplesheet = null
    params.search_params = null

    if (params.samplesheet == null) {
        error "Please provide a samplesheet with the '--samplesheet' option."
    }

    if (params.run.prefilter_proteomes && params.search_params == null) {
        error "prefilter_proteomes requires a samplesheet with the '--search_params' option."
    }

    prefilter_proteomes_fastas = Channel.empty()
    prefilter_proteomes_metamap = Channel.empty()

    // Run the orthology workflow
    if (params.run.prefilter_proteomes) {
        PREFILTER_PROTEOMES (
            params.samplesheet,
            params.search_params,
            params.outdir
        )
        prefilter_proteomes_fastas = PREFILTER_PROTEOMES.out.fastas
        prefilter_proteomes_metamap = PREFILTER_PROTEOMES.out.metamap
    }

    if (params.run.init_ortho_orthofinder) {
        INIT_ORTHO_ORTHOFINDER (
            params.samplesheet,
            params.search_params,
            params.outdir,
            params.mcl_inflation,
            prefilter_proteomes_fastas,
            prefilter_proteomes_metamap
        )
    }
}
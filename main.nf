#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.preview.output = true

include { CHECK_SPECIES_TREE } from './modules/local/check_species_tree'
include { PREFILTER } from './subworkflows/prefilter'
include { INIT_ORTHO as INIT_ORTHO_ORTHOFINDER } from './subworkflows/init_ortho'
include { INIT_ORTHO as INIT_ORTHO_BROCCOLI } from './subworkflows/init_ortho'

workflow {

    params.samplesheet = null
    params.search_params = null

    if (params.samplesheet == null) {
        error "Please provide a samplesheet with the '--samplesheet' option."
    }

    if (params.run.prefilter_hmmsearch && params.search_params == null) {
        log.warn "search_params needs '--search_params' option - running on whole proteome"
    }
    
    if (params.species_tree) {
        CHECK_SPECIES_TREE(
            file(params.samplesheet),
            params.species_tree
        )

        // Handle the log output
        CHECK_SPECIES_TREE.out
            .map { log ->
                def content = log.text
                println "Species Tree Check Log:"
                println content
                if (content.contains("ERROR:")) {
                    error "Species tree validation failed. Pipeline terminated due to missing species in the tree."
                }
            }
    }

    // Run the orthology workflow
    // Always run prefilter step because it also preprocesses fastas even if no seach step is run
    PREFILTER (
        params.samplesheet,
        params.search_params,
        params.outdir,
        params.run.prefilter_hmmsearch
    )

    if (params.run.orthofinder) {
        INIT_ORTHO_ORTHOFINDER (
            params.samplesheet,
            "orthofinder_results",
            params.mcl_inflation,
            PREFILTER.out.fasta_info_metamap,
            PREFILTER.out.cleanfastas_collected,
            "orthofinder"
        )
    }

    if (params.run.broccoli) {
        INIT_ORTHO_BROCCOLI (
            params.samplesheet,
            "broccoli_results",
            params.mcl_inflation,
            PREFILTER.out.fasta_info_metamap,
            PREFILTER.out.cleanfastas_collected,
            "broccoli"
        )
    }

}


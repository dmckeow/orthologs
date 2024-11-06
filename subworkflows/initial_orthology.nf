#!/usr/bin/env nextflow

// Include modules
include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main'

workflow RUN_ORTHOFINDER {
    take:
    fasta_dir
    prior_run

    main:
    // Create a channel for the FASTA files
    fasta_ch = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
        .collect()
        .map { files -> [[id: "orthofinder"], files] }
        .view { "FASTA_CH: $it" }

    // Create a channel for the optional prior run
    prior_run_ch = prior_run
        ? Channel.value([[id: 'prior_run'], file(prior_run)])
        : Channel.value([[],[]])
        prior_run_ch.view { "PRIOR_RUN_CH: $it" }

    ORTHOFINDER(fasta_ch, prior_run_ch)

    emit:
    orthofinder = ORTHOFINDER.out.orthofinder
    working = ORTHOFINDER.out.working
    versions = ORTHOFINDER.out.versions
    
}
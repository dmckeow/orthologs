#!/usr/bin/env nextflow

include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main'
include { BROCCOLI } from '../modules/local/broccoli/main'
include { FILTER_ORTHOGROUPS } from '../modules/local/filter_orthogroups'
include { SEARCH as SEARCH_ORTHOFINDER } from '../modules/local/search/search'
include { SEARCH as SEARCH_BROCCOLI } from '../modules/local/search/search'

workflow ORTHOLOGY {
    take:
    fasta_dir
    orthofinder_prior_run
    orthofinder_min_sequences
    broccoli_args
    broccoli_min_sequences
    search_gene_family_info
    search_gene_family_name
    search_hmm_dir
    outdir
    run_orthofinder
    run_broccoli
    run_search

    main:
    ch_versions = Channel.empty()
    ch_orthofinder_results = Channel.empty()
    ch_broccoli_results = Channel.empty()
    ch_search_results = Channel.empty()

    if (run_orthofinder) {
        // Create a channel for the FASTA files
        fasta_ch = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
            .collect()
            .map { files -> [[id: "orthofinder"], files] }

        // Create a channel for the optional prior run
        prior_run_ch = orthofinder_prior_run
            ? Channel.value([[id: 'prior_run'], file(orthofinder_prior_run)])
            : Channel.value([[],[]])

        ORTHOFINDER(fasta_ch, prior_run_ch)
        ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

        // Create a channel for the Orthogroup_Sequences folder
        orthogroup_sequences = ORTHOFINDER.out.orthofinder
            .map { meta, path -> 
                def orthogroup_dir = file("${path}/Orthogroup_Sequences")
                return [meta, orthogroup_dir]
            }
        
        // Filter orthogroups by min number sequences
        FILTER_ORTHOGROUPS(orthogroup_sequences, orthofinder_min_sequences)
        ch_orthofinder_results = FILTER_ORTHOGROUPS.out.filtered_orthogroups

        ch_orthofinder_fastas = ch_orthofinder_results
            .map { meta, dir -> file("${dir}/*.fa") }
            .flatten()
            .map { file -> [ [id: file.baseName], file ] }

        if (run_search) {
            SEARCH_ORTHOFINDER(
                ch_orthofinder_fastas,
                search_gene_family_info,
                search_gene_family_name,
                file(search_hmm_dir),
                "orthofinder",
                outdir
            )
            ch_versions = ch_versions.mix(SEARCH_ORTHOFINDER.out.versions)
            ch_search_results = ch_search_results.mix(SEARCH_ORTHOFINDER.out.domfasta)
        } else {
            ch_search_results = ch_search_results.mix(ch_orthofinder_fastas)
        }
    }

    if (run_broccoli) {
        // Create a channel for the FASTA files
        fasta_ch = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
            .collect()
            .map { files -> [[id: "broccoli"], files] }

        BROCCOLI(
            fasta_ch,
            broccoli_args
        )
        ch_versions = ch_versions.mix(BROCCOLI.out.versions)
        
        // Add metadata to the orthologous_groups_sequences output
        orthogroups_with_meta = BROCCOLI.out.orthologous_groups_sequences
            .map { path -> 
                [[id: "broccoli_orthogroups"], path]
            }

        // Filter orthogroups by min number sequences
        FILTER_ORTHOGROUPS(orthogroups_with_meta, broccoli_min_sequences)
        ch_broccoli_results = FILTER_ORTHOGROUPS.out.filtered_orthogroups

        ch_broccoli_fastas = ch_broccoli_results
            .map { meta, dir -> file("${dir}/*.fa") }
            .flatten()
            .map { file -> [ [id: file.baseName], file ] }

        if (run_search) {
            SEARCH_BROCCOLI(
                ch_broccoli_fastas,
                search_gene_family_info,
                search_gene_family_name,
                file(search_hmm_dir),
                "broccoli",
                outdir
            )
            ch_versions = ch_versions.mix(SEARCH_BROCCOLI.out.versions)
            ch_search_results = ch_search_results.mix(SEARCH_BROCCOLI.out.domfasta)
        } else {
            ch_search_results = ch_search_results.mix(ch_broccoli_fastas)
        }
    }

    emit:
    versions = ch_versions
    orthofinder_results = ch_orthofinder_results
    broccoli_results = ch_broccoli_results
    search_results = ch_search_results
}
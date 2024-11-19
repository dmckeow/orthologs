#!/usr/bin/env nextflow

// Include modules
include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main'
include { BROCCOLI } from '../modules/local/broccoli/main'
include { FILTER_ORTHOGROUPS} from '../modules/local/filter_orthogroups'

workflow WF_ORTHOFINDER {
    take:
    fasta_dir
    prior_run
    min_sequences

    main:
    // Create a channel for the FASTA files
    fasta_ch = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
        .collect()
        .map { files -> [[id: "orthofinder"], files] }

    // Create a channel for the optional prior run
    prior_run_ch = prior_run
        ? Channel.value([[id: 'prior_run'], file(prior_run)])
        : Channel.value([[],[]])

    ORTHOFINDER(fasta_ch, prior_run_ch)

    // Create a channel for the Orthogroup_Sequences folder
    orthogroup_sequences = ORTHOFINDER.out.orthofinder
        .map { meta, path -> 
            def orthogroup_dir = file("${path}/Orthogroup_Sequences")
            return [meta, orthogroup_dir]
        }
    
    // Filter orthogroups by min number sequences
    FILTER_ORTHOGROUPS(orthogroup_sequences, min_sequences)

    emit:
    orthofinder = ORTHOFINDER.out.orthofinder
    working = ORTHOFINDER.out.working
    versions = ORTHOFINDER.out.versions
    orthogroup_sequences = FILTER_ORTHOGROUPS.out.filtered_orthogroups
    
}

workflow WF_BROCCOLI {
    take:
    fasta_dir
    broccoli_args
    min_sequences

    main:
    // Create a channel for the FASTA files
    fasta_ch = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
        .collect()
        .map { files -> [[id: "broccoli"], files] }

    BROCCOLI (
        fasta_ch,
        broccoli_args
    )
    
    // Add metadata to the orthologous_groups_sequences output
    orthogroups_with_meta = BROCCOLI.out.orthologous_groups_sequences
        .map { path -> 
            [[id: "broccoli_orthogroups"], path]
        }

    // Filter orthogroups by min number sequences
    FILTER_ORTHOGROUPS(orthogroups_with_meta, min_sequences)

    emit:
    orthologous_groups_sequences = FILTER_ORTHOGROUPS.out.filtered_orthogroups
    broccoli = BROCCOLI.out.broccoli

}
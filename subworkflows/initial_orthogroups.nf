#!/usr/bin/env nextflow

include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main'
include { BROCCOLI } from '../modules/local/broccoli/main'
include { FILTER_ORTHOGROUPS as FILTER_ORTHOGROUPS_ORTHOFINDER} from '../modules/local/filter_orthogroups/filter_orthogroups'
include { FILTER_ORTHOGROUPS as FILTER_ORTHOGROUPS_BROCCOLI} from '../modules/local/filter_orthogroups/filter_orthogroups'
include { SEARCH } from '../modules/local/search/search'

// DMND MMSEQS

include { CLUSTER_DMND_MCL } from '../modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER } from '../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { PARSE_MMSEQS_TO_FASTA } from '../modules/local/cluster_mmseqs/parse_mmseqs'

// SEARCH FIRST VERSION

process CONCATENATE_FASTAS {
    input:
    tuple val(meta), path(fasta_files)

    output:
    tuple val(meta), path("combined_sequences.fasta"), emit: combined_fasta

    script:
    """
    cat ${fasta_files} > combined_sequences.fasta
    """
}


workflow INITIAL_ORTHOGROUPS {
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
    cluster_dmnd_args // DMND MMSEQS
    cluster_mcl_args
    cluster_mcl_inflation
    run_cluster_dmnd_mcl
    run_cluster_mmseqs

    main:
    ch_versions = Channel.empty()
    ch_orthofinder_results = Channel.empty()
    ch_broccoli_results = Channel.empty()
    ch_input_fastas = Channel.empty()
    ch_dmnd_mcl_results = Channel.empty() // DMND MMSEQS
    ch_mmseqs_results = Channel.empty()

    // Create a channel for the FASTA files
    ch_fasta_files = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
        .map { file -> [ [id: file.baseName], file ] }

    if (run_search) {
        SEARCH(
            ch_fasta_files,
            search_gene_family_info,
            search_gene_family_name,
            file(search_hmm_dir),
            "initial",
            outdir
        )
        //ch_versions = ch_versions.mix(SEARCH.out.versions)
        ch_input_fastas = SEARCH.out.domfasta
    } else {
        ch_input_fastas = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
            //.collect()
            //.map { files -> [[id: "initial"], files] }
            .map { file -> [ [id: file.baseName], file ] }
    }

    //ch_input_fastas.view { it -> "ch_input_fastas: $it" }

    if (run_orthofinder) {
        // Collect all input fastas for OrthoFinder
        ch_orthofinder_input = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> 
                [[id: "orthofinder"], files]  // Create the structure OrthoFinder expects
            }

        //ch_orthofinder_input.view { it -> "ch_orthofinder_input: $it" }

        // Create a channel for the optional prior run
        prior_run_ch = orthofinder_prior_run
            ? Channel.value([[id: 'prior_run'], file(orthofinder_prior_run)])
            : Channel.value([[],[]])

        ORTHOFINDER(ch_orthofinder_input, prior_run_ch)
        ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

        // Create a channel for the Orthogroup_Sequences folder
        orthogroup_sequences = ORTHOFINDER.out.orthofinder
            .map { meta, path -> 
                def orthogroup_dir = file("${path}/Orthogroup_Sequences")
                return [meta, orthogroup_dir]
            }
        
        // Filter orthogroups by min number sequences
        FILTER_ORTHOGROUPS_ORTHOFINDER(orthogroup_sequences, orthofinder_min_sequences)
        ch_orthofinder_results = FILTER_ORTHOGROUPS_ORTHOFINDER.out.filtered_orthogroups

        //ch_orthofinder_results.view { it -> "ch_orthofinder_results: $it" }
    }

    if (run_broccoli) {
        // Collect all input fastas for Broccoli
        ch_broccoli_input = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> 
                [[id: "broccoli"], files]  // Create the structure OrthoFinder expects
            }

        BROCCOLI(
            ch_broccoli_input,
            broccoli_args
        )
        ch_versions = ch_versions.mix(BROCCOLI.out.versions)
        
        // Add metadata to the orthologous_groups_sequences output
        orthogroups_with_meta = BROCCOLI.out.orthologous_groups_sequences
            .map { path -> 
                [[id: "broccoli_orthogroups"], path]
            }

        // Filter orthogroups by min number sequences
        FILTER_ORTHOGROUPS_BROCCOLI(orthogroups_with_meta, broccoli_min_sequences)
        ch_broccoli_results = FILTER_ORTHOGROUPS_BROCCOLI.out.filtered_orthogroups

        //ch_broccoli_results.view { it -> "ch_broccoli_results: $it" }
    }

    if (run_cluster_dmnd_mcl) {
        // Group all input fastas for CLUSTER_DMND_MCL
        ch_dmnd_mcl_input = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> [[id: "combined"], files] }  // Add metadata

        CLUSTER_DMND_MCL(
            ch_dmnd_mcl_input, // multi input
            //ch_input_fastas, // to test single fasta input
            cluster_dmnd_args,
            cluster_mcl_args,
            cluster_mcl_inflation,
            "dmnd_mcl"
        )
        ch_dmnd_mcl_results = CLUSTER_DMND_MCL.out.dmnd_mcl_fastas
    }

// DMND MMSEQS
    if (run_cluster_mmseqs) {
        
        // Collect all input fastas for concatenation
        ch_fasta_for_concat = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> 
                [[id: "combined"], files]
            }
        //ch_fasta_for_concat.view { it -> "ch_fasta_for_concat: $it" }

        // Concatenate all FASTA files
        CONCATENATE_FASTAS(ch_fasta_for_concat)

        CONCATENATE_FASTAS.out.combined_fasta.view { it -> "CONCATENATE_FASTAS.out.combined_fasta: $it" }

        // Create MMseqs2 database
        MMSEQS_CREATEDB(CONCATENATE_FASTAS.out.combined_fasta)

        // Cluster sequences
        MMSEQS_CLUSTER(MMSEQS_CREATEDB.out.db)

        MMSEQS_CLUSTER.out.db_cluster.view { "MMSEQS_CLUSTER output: $it" }
        MMSEQS_CREATEDB.out.db.view { "MMSEQS_CREATEDB output: $it" }

        MMSEQS_CREATETSV(
            MMSEQS_CLUSTER.out.db_cluster,
            [[:],[]], //MMSEQS_CREATEDB.out.db,
            MMSEQS_CREATEDB.out.db
        )

        PARSE_MMSEQS_TO_FASTA(
            MMSEQS_CREATETSV.out.tsv,
            CONCATENATE_FASTAS.out.combined_fasta
        )
        
        

        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)
        ch_versions = ch_versions.mix(MMSEQS_CLUSTER.out.versions)
        
        //ch_mmseqs_results = PARSE_MMSEQS_TO_FASTA.out.fasta
    }

    emit:
    versions = ch_versions
    orthofinder_results = ch_orthofinder_results
    broccoli_results = ch_broccoli_results
    input_fastas = ch_input_fastas
    dmnd_mcl_results = ch_dmnd_mcl_results // DMND MMSEQS
    //mmseqs_results = ch_mmseqs_results
}
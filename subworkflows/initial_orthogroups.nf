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

        ch_orthofinder_input.view { it -> "ch_orthofinder_input: $it" }

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

        ch_orthofinder_results.view { it -> "ch_orthofinder_results: $it" }
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

        ch_broccoli_results.view { it -> "ch_broccoli_results: $it" }
    }

    if (run_cluster_dmnd_mcl) {
        CLUSTER_DMND_MCL(
            ch_input_fastas,
            cluster_dmnd_args,
            cluster_mcl_args,
            cluster_mcl_inflation,
            "initial/dmnd_mcl"
        )
        //ch_versions = ch_versions.mix(CLUSTER_DMND_MCL.out.versions)
        ch_dmnd_mcl_results = CLUSTER_DMND_MCL.out.dmnd_mcl_fastas
    }

// DMND MMSEQS
    if (run_cluster_mmseqs) {
        // Add the common prefix and .db suffix for CREATEDB, keeping original_id
        ch_input_for_createdb = ch_input_fastas.map { meta, fasta -> 
            [ meta + [id: "${meta.id}.db", original_id: meta.id], fasta ]
        }

        // Create MMseqs2 database
        MMSEQS_CREATEDB(ch_input_for_createdb)

        // Add .cluster suffix for CLUSTER, keeping original_id
        ch_input_for_cluster = MMSEQS_CREATEDB.out.db.map { meta, db ->
            [ meta + [id: "${meta.id.replace('.db', '.cluster')}"], db ]
        }

        // Cluster sequences
        MMSEQS_CLUSTER(ch_input_for_cluster)

        ch_parse_input = MMSEQS_CLUSTER.out.db_cluster
            .map { meta, db -> [meta.original_id, [meta, db]] }
            .join(MMSEQS_CREATEDB.out.db.map { meta, db -> [meta.original_id, [meta, db]] })
            .join(ch_input_fastas.map { meta, fasta -> [meta.id, [meta, fasta]] })
            .map { original_id, cluster_data, createdb_data, fasta_data ->
                def cluster_meta = cluster_data[0]
                def cluster_db = cluster_data[1]
                def query_db = createdb_data[1]
                def fasta = fasta_data[1]
                def result_meta = [id: "${original_id}"]
            
            [
                [meta: [id: original_id], querydb: query_db],
                [meta: [id: original_id], targetdb: query_db],
                [meta: result_meta, resultdb: cluster_db],
                [meta: [id: original_id], fasta: fasta]  // Include the original FASTA file
            ]
        }

        
        // Run PARSE_MMSEQS_TO_FASTA
        PARSE_MMSEQS_TO_FASTA(
            ch_parse_input.map { it[0] },  // querydb
            ch_parse_input.map { it[1] },  // targetdb
            ch_parse_input.map { it[2] },  // resultdb
            ch_parse_input.map { it[3] },   // fasta
            "initial"
        )

        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)
        ch_versions = ch_versions.mix(MMSEQS_CLUSTER.out.versions)
        //ch_versions = ch_versions.mix(PARSE_MMSEQS_TO_FASTA.out.versions)
        ch_mmseqs_results = PARSE_MMSEQS_TO_FASTA.out.fasta
    }

    emit:
    versions = ch_versions
    orthofinder_results = ch_orthofinder_results
    broccoli_results = ch_broccoli_results
    input_fastas = ch_input_fastas
    dmnd_mcl_results = ch_dmnd_mcl_results // DMND MMSEQS
    mmseqs_results = ch_mmseqs_results
}
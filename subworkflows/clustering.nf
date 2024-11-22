#!/usr/bin/env nextflow

include { CLUSTER_DMND_MCL } from '../modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER } from '../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { PARSE_MMSEQS_TO_FASTA } from '../modules/local/cluster_mmseqs/parse_mmseqs'

workflow WF_CLUSTERING {
    take:
    ch_input_fasta
    cluster_dmnd_args
    cluster_mcl_args
    cluster_mcl_inflation
    run_cluster_dmnd_mcl
    run_cluster_mmseqs
    input_source

    main:
    ch_versions = Channel.empty()
    ch_dmnd_mcl_results = Channel.empty()
    ch_mmseqs_results = Channel.empty()

    if (run_cluster_dmnd_mcl) {
        CLUSTER_DMND_MCL(
            ch_input_fasta,
            cluster_dmnd_args,
            cluster_mcl_args,
            cluster_mcl_inflation,
            "${input_source}/dmnd_mcl"
        )
        ch_versions = ch_versions.mix(CLUSTER_DMND_MCL.out.versions)
        ch_dmnd_mcl_results = CLUSTER_DMND_MCL.out.clusters
    }

    if (run_cluster_mmseqs) {
        // Add the common prefix and .db suffix for CREATEDB, keeping original_id
        ch_input_for_createdb = ch_input_fasta.map { meta, fasta -> 
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
            .join(ch_input_fasta.map { meta, fasta -> [meta.id, [meta, fasta]] })
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
            input_source
        )

        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)
        ch_versions = ch_versions.mix(MMSEQS_CLUSTER.out.versions)
        ch_versions = ch_versions.mix(PARSE_MMSEQS_TO_FASTA.out.versions)
        ch_mmseqs_results = PARSE_MMSEQS_TO_FASTA.out.fasta
    }

    emit:
    versions = ch_versions
    dmnd_mcl_results = ch_dmnd_mcl_results
    mmseqs_results = ch_mmseqs_results
}
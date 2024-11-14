include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER } from '../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { PARSE_MMSEQS_TO_FASTA } from '../modules/local/cluster_mmseqs/parse_mmseqs'


workflow WF_CLUSTER_MMSEQSX {
    take:
    ch_input_fasta // Channel with input FASTA files: [ meta, fasta ]
    val_prefix     // Single value for the prefix to be applied to all runs

    main:
    // Add the common prefix and .db suffix for CREATEDB
    ch_input_for_createdb = ch_input_fasta.map { meta, fasta -> 
        [ meta + [id: "${val_prefix}_${meta.id}.db"], fasta ]
    }

    // Create MMseqs2 database
    MMSEQS_CREATEDB (
        ch_input_for_createdb
    )

    // Add .cluster suffix for CLUSTER
    ch_input_for_cluster = MMSEQS_CREATEDB.out.db.map { meta, db ->
        [ meta + [id: "${meta.id.replace('.db', '.cluster')}"], db ]
    }

    // Cluster sequences
    MMSEQS_CLUSTER ( 
        ch_input_for_cluster
    )

    emit:
    //tsv = MMSEQS_CREATETSV.out.tsv
    db  = MMSEQS_CREATEDB.out.db
    cluster_db = MMSEQS_CLUSTER.out.db_cluster
}

workflow WF_CLUSTER_MMSEQS {
    take:
    ch_input_fasta // Channel with input FASTA files: [ meta, fasta ]
    val_prefix     // Single value for the prefix to be applied to all runs

    main:
    // Add the common prefix and .db suffix for CREATEDB, keeping original_id
    ch_input_for_createdb = ch_input_fasta.map { meta, fasta -> 
        [ meta + [id: "${val_prefix}_${meta.id}.db", original_id: meta.id], fasta ]
    }

    // Create MMseqs2 database
    MMSEQS_CREATEDB (
        ch_input_for_createdb
    )

    // Add .cluster suffix for CLUSTER, keeping original_id
    ch_input_for_cluster = MMSEQS_CREATEDB.out.db.map { meta, db ->
        [ meta + [id: "${meta.id.replace('.db', '.cluster')}"], db ]
    }

    // Cluster sequences
    MMSEQS_CLUSTER ( 
        ch_input_for_cluster
    )

    // Prepare input for PARSE_MMSEQS_TO_FASTA
    ch_parse_input = MMSEQS_CLUSTER.out.db_cluster
        .join(MMSEQS_CREATEDB.out.db, by: 0)  // Join by the entire meta map
        .map { cluster_meta, cluster_db, query_db ->
            def resultdb_meta = cluster_meta + [id: "${cluster_meta.id}_result"]
            log.info "PARSE input: cluster_meta=${cluster_meta}, query_meta=${cluster_meta}, resultdb_meta=${resultdb_meta}"
            [
                [cluster_meta, cluster_db],
                [cluster_meta, query_db],
                [resultdb_meta, cluster_db]  // Using cluster_db as resultdb
            ]
        }

    // Run PARSE_MMSEQS_TO_FASTA
    PARSE_MMSEQS_TO_FASTA (
        ch_parse_input.map { it[0] },  // querydb
        ch_parse_input.map { it[1] },  // targetdb
        ch_parse_input.map { it[2] }   // resultdb
    )

    emit:
    db  = MMSEQS_CREATEDB.out.db
    cluster_db = MMSEQS_CLUSTER.out.db_cluster
    abc = PARSE_MMSEQS_TO_FASTA.out.abc
    tsv = PARSE_MMSEQS_TO_FASTA.out.tsv
}
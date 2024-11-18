include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER } from '../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { PARSE_MMSEQS_TO_FASTA } from '../modules/local/cluster_mmseqs/parse_mmseqs'

workflow WF_CLUSTER_MMSEQS {
    take:
    ch_input_fasta // Channel with input FASTA files: [ meta, fasta ]
    input_source

    main:
    // Add the common prefix and .db suffix for CREATEDB, keeping original_id
    ch_input_for_createdb = ch_input_fasta.map { meta, fasta -> 
        [ meta + [id: "${meta.id}.db", original_id: meta.id], fasta ]
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

    // View the contents of ch_parse_input
    def logFile = file("${params.outdir}/ch_parse_input_logs.txt")
    ch_parse_input
    .map { querydb_tuple, targetdb_tuple, resultdb_tuple, fasta_tuple ->
        "\n\nch_parse_input:\n" +
        "  Query DB - Meta: ${querydb_tuple.meta}, Path: ${querydb_tuple.querydb}\n" +
        "  Target DB - Meta: ${targetdb_tuple.meta}, Path: ${targetdb_tuple.targetdb}\n" +
        "  Result DB - Meta: ${resultdb_tuple.meta}, Path: ${resultdb_tuple.resultdb}\n" +
        "  FASTA - Meta: ${fasta_tuple.meta}, Path: ${fasta_tuple.fasta}"
    }
    .collectFile(name: logFile.getName(), storeDir: logFile.getParent(), newLine: false)

    // Run PARSE_MMSEQS_TO_FASTA
    PARSE_MMSEQS_TO_FASTA (
        ch_parse_input.map { it[0] },  // querydb
        ch_parse_input.map { it[1] },  // targetdb
        ch_parse_input.map { it[2] },  // resultdb
        ch_parse_input.map { it[3] },   // fasta
        input_source
    )

    emit:
    db  = MMSEQS_CREATEDB.out.db
    cluster_db = MMSEQS_CLUSTER.out.db_cluster
    abc = PARSE_MMSEQS_TO_FASTA.out.abc
    tsv = PARSE_MMSEQS_TO_FASTA.out.tsv
    fasta = PARSE_MMSEQS_TO_FASTA.out.fasta
}
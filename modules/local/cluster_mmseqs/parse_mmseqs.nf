process PARSE_MMSEQS_TO_FASTA {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(querydb)
    tuple val(meta2), path(targetdb)
    tuple val(meta3), path(resultdb)

    output:
    tuple val(meta), path("*.abc"), emit: abc
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    """
    python {projectDir}/bin/parse_mmseqs.py
        --querydb ${querydb} \
        --targetdb ${targetdb} \
        --resultdb ${resultdb}
    """
}
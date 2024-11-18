process PARSE_MMSEQS_TO_FASTA {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(querydb)
    tuple val(meta2), path(targetdb, stageAs: 'targetdb/*')
    tuple val(meta3), path(resultdb)
    tuple val(meta4), path(fasta_file)
    val(input_source)

    output:
    tuple val(meta), path("${input_source}/${resultdb}.abc"), emit: abc
    tuple val(meta), path("${input_source}/${resultdb}.tsv"), emit: tsv
    tuple val(meta), path("${input_source}/*.*.mmseqs.fa"), emit: fasta

    script:
    def targetdb_name = targetdb.getBaseName()
    def targetdb_ext = targetdb.getExtension()
    """
    python ${projectDir}/bin/parse_mmseqs.py \
        -q ${querydb}/${querydb} \
        -t targetdb/${targetdb_name}.${targetdb_ext}/${targetdb_name}.${targetdb_ext} \
        -r ${resultdb}/${resultdb} \
        -o ${input_source}/
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m ${input_source}/${resultdb}.abc \
        -f ${fasta_file} \
        -o ${input_source}/
    """
}
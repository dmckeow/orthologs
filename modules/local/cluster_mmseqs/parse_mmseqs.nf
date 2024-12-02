process PARSE_MMSEQS_TO_FASTA {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(sequenceDB)
    tuple val(meta), path(resultDB)
    tuple val(meta), path(fasta_file)
    val(source)
    
    output:
    tuple val(meta), path("*.fa"),                      emit: fasta
    tuple val(meta), path("mmseqs.orthogroups.txt"),    emit: orthogroups

    script:
    def createseqfiledb_args = task.ext.createseqfiledb_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mmseqs createtsv \
    ${sequenceDB}/${sequenceDB} \
    ${sequenceDB}/${sequenceDB} \
    ${resultDB}/${resultDB} \
    clusters.tsv \
    --threads ${task.cpus}

    mmseqs createseqfiledb \
    ${sequenceDB}/${sequenceDB} \
    ${resultDB}/${resultDB} \
    sequences \
    --threads ${task.cpus} \
    $createseqfiledb_args

    python3 ${projectDir}/bin/parse_mmseqs.py \
        -i clusters.tsv
    
    python3 ${projectDir}/bin/parse_fastas_mcl.py \
        -m clusters.tsv.abc \
        -f ${fasta_file} \
        -o ./ \
        --source ${source}
    
    mv clusters.abc.tmp mmseqs.orthogroups.txt

    """
}


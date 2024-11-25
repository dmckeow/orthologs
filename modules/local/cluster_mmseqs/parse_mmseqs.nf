process PARSE_MMSEQS_TO_FASTA_1 {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(mmseqstsv)
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path("${resultdb}.abc"), emit: abc
    tuple val(meta), path("*.*.mmseqs.fa"), emit: fasta

    script:
    def mmseqstsv_basename = mmseqstsv.getBaseName()
    
    """
    python ${projectDir}/bin/parse_mmseqs.py \
        -i ${mmseqstsv}
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m ${mmseqstsv_basename}.abc \
        -f ${fasta_file} \
        -o ./
    """
}


process PARSE_MMSEQS_TO_FASTA {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(sequenceDB)
    tuple val(meta), path(resultDB)
    tuple val(meta), path(fasta_file)
    val(source)
    
    output:
    tuple val(meta), path("*.cluster.fa"), emit: fasta

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

    python ${projectDir}/bin/parse_mmseqs.py \
        -i clusters.tsv
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m clusters.tsv.abc \
        -f ${fasta_file} \
        -o ./ \
        --source ${source}

    """
}


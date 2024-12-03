process PARSE_MMSEQS_TO_FASTA {

    conda "${moduleDir}/environment.yml"
    
    //container 'oras://community.wave.seqera.io/library/mmseqs2_biopython_python:362838f7a003ac69'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f5/f5898f2f43274c4fbcff3e891df7295540f40cece35428b5762b0ac35f6ea955/data'

    input:
    tuple val(meta), path(sequenceDB)
    tuple val(meta), path(resultDB)
    tuple val(meta), path(fasta_file)
    val(source)
    
    output:
    tuple val(meta), path("cluster_sequences/*.fa"),     optional: true,           emit: fasta_files
    tuple val(meta), path("cluster_sequences"), optional: true,                    emit: fasta_dir
    tuple val(meta), path("mmseqs.orthogroups.txt"),     optional: true,           emit: orthogroup_list

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
        -o cluster_sequences/ \
        --source ${source}
    
    mv clusters.abc.tmp mmseqs.orthogroups.txt

    """
}


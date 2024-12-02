process SEARCH {
    conda "${moduleDir}/environment.yml"
    //container 'community.wave.seqera.io/library/bedtools_biopython_hmmer_samtools_python:c2a5749688500b49'
    //container 'oras://community.wave.seqera.io/library/bedtools_biopython_hmmer_samtools_python:6ff18068cb3833ec'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d06f9c1b929618ccc24024d4d13768bd9196a05416f9dcad078853ecbf40efa9/data'

    input:
    tuple val(meta), path(fasta)
    val(gene_family_info)
    val(gene_family_name)
    path(hmm_dir)
    val(input_source)
    val(output_dir)

    output:
    tuple val(meta), path("searches/${input_source}/${meta.id}.txt"), emit: txt
    tuple val(meta), path("searches/${input_source}/*.domains.fasta"), optional: true, emit: domfasta
    tuple val(meta), path("searches/${input_source}/*.domtable"), optional: true, emit: domtable
    
    script:
    """
    python3 ${projectDir}/bin/search.py \
        -f ${fasta} \
        -g ${gene_family_info} \
        -t ${task.cpus} \
        ${gene_family_name} \
        -H ${hmm_dir} \
        -i ${input_source} \
        -o searches

    touch searches/${input_source}/${meta.id}.txt
    """
}
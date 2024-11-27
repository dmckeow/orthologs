process SEARCH_OG {
    conda "${projectDir}/modules/local/search/environment.yml"
    
    input:
    tuple val(meta), path(fasta)
    val(gene_family_info)
    val(gene_family_name)
    path(hmm_dir)
    val(input_source)
    val(output_dir)

    output:
    tuple val(meta), path("searches/${input_source}/${fasta.baseName}.txt"), emit: txt

    tuple val(meta), path("searches/${input_source}/*.domains.fasta"), optional: true, emit: domfasta
    tuple val(meta), path("searches/${input_source}/*.domtable"), optional: true, emit: domtable
    
    script:
    """
    python ${projectDir}/bin/search.py \
        -f ${fasta} \
        -g ${gene_family_info} \
        -t ${task.cpus} \
        ${gene_family_name} \
        -H ${hmm_dir} \
        -i ${input_source} \
        -o searches

    touch searches/${input_source}/${fasta.baseName}.txt
    """
}
// SEARCH process with SASH:
process SEARCH {
    conda "${projectDir}/modules/local/search/environment.yml"
    
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
    python ${projectDir}/bin/search.py \
        -f ${fasta} \
        -g ${gene_family_info} \
        -t ${task.cpus} \
        ${gene_family_name} \
        -H ${hmm_dir} \
        -i ${input_source} \
        -o searches \
        --sample-name ${meta.id}

    touch searches/${input_source}/${meta.id}.txt
    """
}
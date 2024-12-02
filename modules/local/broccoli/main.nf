process BROCCOLI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    //container 'community.wave.seqera.io/library/biopython_diamond_fasttree_ete3_python:c345017048f04a9c'
    //container 'oras://community.wave.seqera.io/library/biopython_diamond_fasttree_ete3_python:50a4e9ceccbb9484'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')
    val(args)
    
    output:
    path("dir_step3/orthologous_groups.txt"), emit: orthologous_groups
    path("dir_step3/table_OGs_protein_names.txt"), emit: table_OGs_protein_names
    path("**"), emit: broccoli
    path("dir_step3/orthologous_groups_sequences"), emit: orthologous_groups_sequences
    path("dir_step4/orthologous_pairs.txt"), emit: orthologous_pairs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run Broccoli
    python3 ${projectDir}/broccoli/broccoli.py \\
        -dir input \\
        -threads ${task.cpus} \\
        $args
    
    python3 ${projectDir}/bin/parse_fastas_broccoli.py \\
        -b dir_step3/orthologous_groups.txt \\
        -f input \\
        -o dir_step3/orthologous_groups_sequences

    rm -fr input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        broccoli: \$(python ${projectDir}/broccoli/broccoli.py --version 2>&1 | sed 's/^.*Broccoli //; s/ .*\$//')
    END_VERSIONS
    """
}
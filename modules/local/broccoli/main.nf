process BROCCOLI {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/broccoli/environment.yml"

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
    python ${projectDir}/broccoli/broccoli.py \\
        -dir input \\
        -threads ${task.cpus} \\
        $args
    
    python ${projectDir}/bin/parse_fastas_broccoli.py \\
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
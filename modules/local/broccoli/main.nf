process BROCCOLI {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/broccoli/environment.yml"

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')
    val(args)
    
    output:
    path("broccoli/**"), emit: broccoli
    path("broccoli/dir_step3/orthologous_groups.txt"), emit: orthologous_groups
    path("broccoli/dir_step4/orthologous_pairs.txt"), emit: orthologous_pairs
    path "broccoli/versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run Broccoli
    rm -fr broccoli
    mkdir broccoli
    cd broccoli
    python ${projectDir}/broccoli/broccoli.py \\
        -dir ../input \\
        -threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        broccoli: \$(python ${projectDir}/broccoli/broccoli.py --version 2>&1 | sed 's/^.*Broccoli //; s/ .*\$//')
    END_VERSIONS
    """
}
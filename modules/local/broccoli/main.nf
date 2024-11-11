process BROCCOLI {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/broccoli/environment.yml"

    input:
    tuple val(meta), path(proteomes)
    
    output:
    tuple val(meta), path("*.orthologs.txt"), emit: orthologs
    tuple val(meta), path("*.tree"), emit: tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Clone Broccoli repository
    git clone https://github.com/rderelle/Broccoli.git
    cd broccoli

    # Run Broccoli
    python broccoli.py \\
        -i $proteomes \\
        -o ${prefix} \\
        $args

    # Move output files to the current directory
    mv ${prefix}* ..
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        broccoli: \$(python broccoli/broccoli.py --version 2>&1 | sed 's/^.*Broccoli //; s/ .*\$//')
    END_VERSIONS
    """
}
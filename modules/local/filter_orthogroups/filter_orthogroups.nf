process FILTER_ORTHOGROUPS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::biopython"

    input:
    tuple val(meta), path(orthogroup_dir)
    val min_sequences

    output:
    tuple val(meta), path("filtered_orthogroups"), emit: filtered_orthogroups
    tuple val(meta), path("removed_orthogroups"), emit: removed_orthogroups
    path "versions.yml", emit: versions

    script:
    """
    python ${projectDir}/bin/filter_orthogroups.py ${orthogroup_dir} ${min_sequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
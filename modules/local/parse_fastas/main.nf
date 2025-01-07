process PARSE_FASTAS {
    conda "bioconda::biopython"

    input:
    tuple val(meta), path(fastas, stageAs: 'input_fasta_dir/')
    tuple val(meta), path(og_file)
    val out_dir_name
    
    output:
    path "${out_dir_name}/*", emit: parsed_fastas
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${out_dir_name}

    python3 ${projectDir}/bin/parse_fastas.py \\
        -i ${og_file} \\
        -f input_fasta_dir \\
        -o ${out_dir_name} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${out_dir_name}
    touch ${out_dir_name}/dummy_output.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
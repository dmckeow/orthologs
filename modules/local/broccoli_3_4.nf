process BROCCOLI {
    tag "Calling initial orthogroups with Broccoli"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/broccoli_j",
        mode: 'copy',
        saveAs: { fn -> fn }
    )

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    path(fastas, stageAs: 'input/')
    val publish_subdir
    
    output:
    path("dir_step3/orthologous_groups.txt"), emit: orthologous_groups
    path("dir_step3/table_OGs_protein_names.txt"), emit: table_OGs_protein_names
    path("dir_step3/Orthogroup_Sequences"), emit: orthologous_groups_sequences
    path("dir_step4/orthologous_pairs.txt"), emit: orthologous_pairs
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Run Broccoli
    python3 ${projectDir}/broccoli/broccoli.py \\
        -dir input \\
        -threads ${task.cpus} \\
        -steps 3,4 \\
        $args
    
    python3 ${projectDir}/bin/parse_fastas.py \\
        -i dir_step3/orthologous_groups.txt \\
        -f input \\
        -o dir_step3/Orthogroup_Sequences

    
    """
}
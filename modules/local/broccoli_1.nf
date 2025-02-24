process BROCCOLI {
    tag "Calling initial orthogroups with Broccoli"
    label 'process_broccoli'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/broccoli",
        mode: 'copy',
        saveAs: { fn -> fn }
    )

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    path(fastas, stageAs: 'input/')
    val publish_subdir
    
    output:
    path("dir_step1/**"), emit: dir_step1
    path("dir_step1/species_index.pic"), emit: species_index
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    
    python3 ${projectDir}/broccoli/broccoli.py \\
        -dir input \\
        -threads ${task.cpus} \\
        -steps 1 \\
        $args
    

    
    """
}
process BROCCOLI {
    tag "Calling initial orthogroups with Broccoli"
    label 'process_broccoli_array'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/broccoli",
        mode: 'copy',
        saveAs: { fn -> fn }
    )

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    path(fastas, stageAs: 'input/')
    path(broccoli_results_dir, stageAs: 'dir_step2/') // from previous step
    path(files_start) // This will be a single line from step2 part1
    val publish_subdir
    
    output:
    path("dir_step2/**"), emit: dir_step2
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    

    python3 ${projectDir}/broccoli/broccoli.py \\
        -dir input \\
        -threads ${task.cpus} \\
        -steps 2 \\
        -sub_step 2 \\
        -sub_step2_input ${files_start} \\
        $args

    
    """
}
process BROCCOLI {
    tag "Building phylome: $files_start"
    label 'process_broccoli_array'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/broccoli",
        mode: 'copy',
        saveAs: { fn -> fn }
    )

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    path(fastas, stageAs: 'input/')
    path(dir_step1, stageAs: 'dir_step1/')
    path(prot_str_2_species, stageAs: 'dir_step2/')
    path(prot_int_2_species, stageAs: 'dir_step2/')
    path(files_start, stageAs: 'dir_step2/') // This will be a single line from step2 part1
    path(databases, stageAs: 'dir_step2/databases/')
    val publish_subdir
    
    output:
    path("dir_step2/dict_output/*.pic"), emit: dict_output
    path("dir_step2/dict_similarity_ortho/*.pic"), emit: dict_similarity_ortho
    path("dir_step2/dict_trees/*.pic"), emit: dict_trees
    
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
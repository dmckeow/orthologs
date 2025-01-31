process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_low'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/orthofinder",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        workflow.containerEngine == 'apptainer' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        '' }"

    input:
    path(fastas, stageAs: 'input/')
    val publish_subdir

    output:
    path "input/OrthoFinder/**.dmnd"           , emit: diamonds
    path "input/OrthoFinder/**.fa"             , emit: fastas
    path "input/OrthoFinder/**SequenceIDs.txt" , emit: seqIDs
    path "input/OrthoFinder/**SpeciesIDs.txt"  , emit: sppIDs

    script:
    """
    
    orthofinder \\
        -f input/ \\
        -t ${task.cpus} \\
        -op > tmp

    """
}

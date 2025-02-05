process FASTTREE {
    tag "$meta.og"
    label 'process_fasttree'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/fasttree_2.1.11:1.0.0':
        workflow.containerEngine == 'apptainer' ? 'arcadiascience/fasttree_2.1.11:1.0.0':
        '' }"

    publishDir(
        path: "${params.outdir}/${publish_subdir}/fasttree_gene_trees",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), file(alignment)
    val model // not used
    val publish_subdir

    output:
    tuple val(meta), path("*.treefile") , emit: phylogeny

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def og   = "${meta.og}"
    """
    # Make sure the number of threads are being specified properly
    export OMP_NUM_THREADS=${task.cpus}

    # Efficiently infer a gene family tree using FastTree!
    FastTreeDblMP \\
        $args \\
        $alignment > ${og}_ft.treefile

    # prevent zero-length branches (sometimes inferred with fasttree)
    ${projectDir}/bin/resolve_polytomies.R ${og}_ft.treefile resolved.tree
    mv resolved.tree ${og}_ft.treefile

    """
}

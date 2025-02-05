process POSSVM {
    tag "Parsing final orthgroups with POSSVM"
    label 'process_possvm'

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c73f545df55c871d54d273a80538bfd03d9f056a5444a556fc63edc0455486a/data'

    publishDir(
        path: "${params.outdir}/${publish_subdir}/${publish_subdir2}/possvm",
        mode: 'copy',
        saveAs: { fn -> fn }
    )

    input:
    tuple val(meta), path(generax_gft)
    val publish_subdir
    val publish_subdir2
    
    output:
    path("**.ortholog_groups.csv"), emit: ogs_csv
    path("**.ortholog_groups.newick"), emit: ogs_nwk
    path("**.ortholog_groups.newick.pdf"), optional: true, emit: ogs_nwk_pdf
    path("**.ortholog_pairs.csv"), optional: true, emit: ogs_pairs_csv
    path("**.pairs_all.csv"), optional: true, emit: ogs_pairs_all_csv

    script:
    def args = task.ext.args ?: ''
    def og   = "${meta.og}"
    """
    python3 ${projectDir}/possvm/possvm.py \\
        -i ${generax_gft} \\
        -ogprefix OG \\
        -p ${og} \\
        ${args}
    """
}
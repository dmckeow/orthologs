process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_high'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        workflow.containerEngine == 'apptainer' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        '' }"

    publishDir (
        path: "${params.outdir}/orthofinder_mcl",
        mode: 'copy',
        pattern: "${output_directory}/Results_Inflation*",
        saveAs: { fn -> fn }
    )

    input:
    each mcl_inflation
    file(blast)
    file(fasta)
    file(db)
    file(sppIDs)
    file(seqIDs)
    val output_directory

    output:
    path("${output_directory}/Results_Inflation_${mcl_inflation}"), emit: inflation_dir
    path("${output_directory}/Results_Inflation_${mcl_inflation}/Orthogroup_Sequences"), emit: initial_orthogroups_fa_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    

    orthofinder \\
        -b ./ \\
        -n "Inflation_${mcl_inflation}" \\
        -I $mcl_inflation \\
        -M msa -X -os -z \\
        -a ${task.cpus} \\
        $args

    
    dir=\$(pwd)
    cd \$(ls -d OrthoFinder/*/WorkingDirectory)
    tar -czvf Sequences_ids.tar.gz Sequences_ids
    rm -r Sequences_ids
    cd \$dir


    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mkdir ${output_directory}
    mv OrthoFinder/Results_Inflation_${mcl_inflation}/ ${output_directory}/
    rm -r OrthoFinder/
    """
}

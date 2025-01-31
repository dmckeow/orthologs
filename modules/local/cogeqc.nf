process COGEQC {
    tag "Orthogroup Summary"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/cogeqc_1.2.1:1.0.0':
        workflow.containerEngine == 'apptainer' ? 'arcadiascience/cogeqc_1.2.1:1.0.0':
        '' }"

    publishDir(
        path: "${params.outdir}/${publish_subdir}/orthogroup_summaries",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file orthofinder_outdir
    val min_spp             // Minimum number of species for orthogroup retention
    file prot_annotations   // Base filepath to where protein annotations are stored
    val publish_subdir

    output:
    path "*_cogeqc_summary.tsv" , emit: cogeqc_summary
    
    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    """
    # Assess orthogroups inferred using each inflation parameter, summarizing
    # how well they group proteins with the same domains together, as well as
    # other summary stats like number of ogs with >= the minimum # species,
    # per-species gene count per-og, etc.
    ${projectDir}/bin/cogeqc_summarize_ogs.R ${orthofinder_outdir} ${min_spp}
    
    """
}

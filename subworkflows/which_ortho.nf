include { GET_ORTHOGROUP_INFO } from '../modules/local/report/get_og_info'

workflow WHICH_ORTHO {
    take:
    combined_deflines
    orthofinder_orthogroups
    broccoli_orthogroups
    dmnd_mcl_clusters
    mmseqs_clusters
    search_status

    main:

    GET_ORTHOGROUP_INFO (
            combined_deflines,
            orthofinder_orthogroups,
            broccoli_orthogroups,
            dmnd_mcl_clusters,
            mmseqs_clusters,
            search_status
        )
    
    
    // Define channel for input FASTA files
    fasta_files = Channel.fromPath(params.input_fasta_dir + "/*.faa") // change to suit the outputs from init_ortho

    // Run orthogroup prediction tools
    ORTHOFINDER(fasta_files.collect())
    BROCCOLI(fasta_files.collect())
    DIAMOND_MCL(fasta_files.collect())
    MMSEQS(fasta_files.collect())

    // Calculate Jaccard indexes for all pairs of tools
    CALCULATE_JACCARD(ORTHOFINDER.out.orthogroups, BROCCOLI.out.orthogroups, 'OrthoFinder', 'Broccoli')
    CALCULATE_JACCARD(ORTHOFINDER.out.orthogroups, DIAMOND_MCL.out.orthogroups, 'OrthoFinder', 'DIAMOND_MCL')
    CALCULATE_JACCARD(ORTHOFINDER.out.orthogroups, MMSEQS.out.orthogroups, 'OrthoFinder', 'MMseqs2')
    CALCULATE_JACCARD(BROCCOLI.out.orthogroups, DIAMOND_MCL.out.orthogroups, 'Broccoli', 'DIAMOND_MCL')
    CALCULATE_JACCARD(BROCCOLI.out.orthogroups, MMSEQS.out.orthogroups, 'Broccoli', 'MMseqs2')
    CALCULATE_JACCARD(DIAMOND_MCL.out.orthogroups, MMSEQS.out.orthogroups, 'DIAMOND_MCL', 'MMseqs2')

   

    
    emit:
    GET_ORTHOGROUP_INFO.out.combined_orthogroups


}

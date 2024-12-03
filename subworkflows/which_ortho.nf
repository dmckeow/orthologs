include { GET_ORTHOGROUP_INFO } from '../modules/local/report/get_og_info'
include { CALCULATE_JACCARD } from '../modules/local/report/calc_jaccard'

workflow WHICH_ORTHO {
    take:
    combined_deflines
    orthofinder_orthogroups
    broccoli_orthogroups
    dmnd_mcl_clusters
    mmseqs_clusters
    search_status
    orthofinder_og_fa_dir
    broccoli_og_fa_dir
    dmnd_mcl_og_fa_dir
    mmseqs_og_fa_dir

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
    
    
    orthofinder_og_fa_dir
        .map { meta, path -> path.listFiles().findAll { it.name.endsWith('.fa') } }
        .flatten()
        .collect()
        .set { ch_collected_fa_files }
    
    broccoli_og_fa_dir
        .map { meta, path -> path.listFiles().findAll { it.name.endsWith('.fa') } }
        .flatten()
        .collect()
        .set { ch_collected_fa_files }

    // Calculate Jaccard indexes for all pairs of tools
    CALCULATE_JACCARD(orthofinder_og_fa_dir, broccoli_og_fa_dir, 'OrthoFinder', 'Broccoli')
    //CALCULATE_JACCARD(ORTHOFINDER.out.orthogroups, DIAMOND_MCL.out.orthogroups, 'OrthoFinder', 'DIAMOND_MCL')
    //CALCULATE_JACCARD(ORTHOFINDER.out.orthogroups, MMSEQS.out.orthogroups, 'OrthoFinder', 'MMseqs2')
    //CALCULATE_JACCARD(BROCCOLI.out.orthogroups, DIAMOND_MCL.out.orthogroups, 'Broccoli', 'DIAMOND_MCL')
    //CALCULATE_JACCARD(BROCCOLI.out.orthogroups, MMSEQS.out.orthogroups, 'Broccoli', 'MMseqs2')
    //CALCULATE_JACCARD(DIAMOND_MCL.out.orthogroups, MMSEQS.out.orthogroups, 'DIAMOND_MCL', 'MMseqs2')

   

    
    emit:
    GET_ORTHOGROUP_INFO.out.combined_orthogroups


}

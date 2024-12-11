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
    all_input_fastas

    main:

    GET_ORTHOGROUP_INFO (
        combined_deflines,
        orthofinder_orthogroups,
        broccoli_orthogroups,
        dmnd_mcl_clusters,
        mmseqs_clusters,
        search_status
    )

    
    
    
    orthogroups_file = GET_ORTHOGROUP_INFO.out.orthogroups
    // Calculate overall average jaccard score between tools for og membership by deflines
    CALCULATE_JACCARD (
        orthogroups_file
        )

    

    
    //emit:
    


}

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
    
    emit:
    GET_ORTHOGROUP_INFO.out.combined_orthogroups


}

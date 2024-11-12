#!/usr/bin/env nextflow

// Include modules
include { WF_ORTHOFINDER } from './subworkflows/initial_orthology'
include { WF_BROCCOLI } from './subworkflows/initial_orthology.nf'
include { SEARCH as SEARCH_ORTHOFINDER } from './modules/local/search/search'
include { SEARCH as SEARCH_BROCCOLI } from './modules/local/search/search'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_ORTHOFINDER } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_BROCCOLI } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'


workflow {
// Orthofinder and search
    if (params.run_orthofinder) {
        WF_ORTHOFINDER(
            params.fasta_dir,
            params.prior_run
        )
        // Set the channel from the orthofinder output
        ch_orthofinder_fastas = WF_ORTHOFINDER.out.orthogroup_sequences
            .map { meta, dir -> file("${dir}/*.fa") }
            .flatten()
            .map { file -> [ [id: file.baseName], file ] }

        if (params.run_search) {
            log.info "Running search as params.run_search is set to true - running clustering on proteins that had search hits and were in orthogroups created by OrthoFinder"
            // Run SEARCH_ORTHOFINDER process
            SEARCH_ORTHOFINDER(
                ch_orthofinder_fastas,
                params.search.gene_family_info,
                params.search.gene_family_name,
                file(params.search.hmm_dir),
                "orthofinder",
                params.outdir
            )
        if (params.run_cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_ORTHOFINDER(
                SEARCH_ORTHOFINDER.out.domfasta,
                params.cluster_dmnd_mcl.dmnd_params,
                params.cluster_dmnd_mcl.mcl_params,
                params.cluster_dmnd_mcl.mcl_inflation,
                "searches/orthofinder"
                )
            }
        } else {
            log.info "Skipping search as params.run_search is set to false - running clustering on all proteins in orthogroups created by OrthoFinder"
            if (params.run_cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_ORTHOFINDER(
                ch_orthofinder_fastas,
                params.cluster_dmnd_mcl.dmnd_params,
                params.cluster_dmnd_mcl.mcl_params,
                params.cluster_dmnd_mcl.mcl_inflation,
                "all/orthofinder"
                )
            }
        }
    }

// Broccoli and search
    if (params.run_broccoli) {
        log.info "Running Broccoli"
        WF_BROCCOLI(
            params.fasta_dir,
            params.broccoli.args
        )
    
        // Set the channel from the orthofinder output
        ch_broccoli_fastas = WF_BROCCOLI.out.orthologous_groups_sequences
        .map { dir -> 
            def files = file("${dir}/*.fa")
            files.collect { file -> [ [id: file.baseName], file ] }
        }
        .flatten()
        .collate(2)

        if (params.run_search) {
            log.info "Running search as params.run_search is set to true - running clustering on proteins that had search hits and were in orthogroups created by Broccoli"
            // Run SEARCH process
            SEARCH_BROCCOLI(
                ch_broccoli_fastas,
                params.search.gene_family_info,
                params.search.gene_family_name,
                file(params.search.hmm_dir),
                "broccoli",
                params.outdir
            )
            if (params.run_cluster_dmnd_mcl) {
                log.info "Running DIAMOND and MCL based clustering within orthogroups identified by broccoli and then search"
                CLUSTER_DMND_MCL_BROCCOLI(
                    SEARCH_BROCCOLI.out.domfasta,
                    params.cluster_dmnd_mcl.dmnd_params,
                    params.cluster_dmnd_mcl.mcl_params,
                    params.cluster_dmnd_mcl.mcl_inflation,
                    "searches/broccoli"
                    )
                }
        } else {
            log.info "Skipping search as params.run_search is set to false - running clustering on all proteins in orthogroups created by Broccoli"
            if (params.run_cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_BROCCOLI(
                ch_broccoli_fastas,
                params.cluster_dmnd_mcl.dmnd_params,
                params.cluster_dmnd_mcl.mcl_params,
                params.cluster_dmnd_mcl.mcl_inflation,
                "all/broccoli"
                )
            }
         }
    }
}

#!/usr/bin/env nextflow


// Include modules
include { BAM_INDEX } from './subworkflows/bam_index'
include { WF_ORTHOFINDER } from './subworkflows/initial_orthology'
include { SEARCH } from './modules/local/search/search'
include { CLUSTER_DMND_MCL } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'

workflow {
    if (params.run_bam_index) {
        BAM_INDEX(params.sample_sheet)
    }
// Orthofinder and search
    if (params.run_orthofinder) {
        WF_ORTHOFINDER(
            params.fasta_dir,
            params.prior_run
        )

        if (params.run_search) {
            // Wait for INITIAL_ORTHOLOGY to finish and use its output
            // Define channels
            ch_fasta = WF_ORTHOFINDER.out.orthogroup_sequences
                .map { meta, dir -> file("${dir}/*.fa") }
                //.view { "Files in Orthogroup_Sequences: $it" }
                .flatten()
                .map { file -> [ [id: file.baseName], file ] }
            
            // Run SEARCH_ORTHOFINDER process
            SEARCH(
                ch_fasta,
                params.search.gene_family_info,
                params.search.gene_family_name,
                file(params.search.hmm_dir),
                "orthofinder",
                params.outdir
            )
        if (params.run_cluster_dmnd_mcl) {
            // Run dmnd mcl clustering on search results only
            //ch_fasta = SEARCH.out.domfasta
              //  .map { meta, dir -> file("${dir}/*.domains.fasta") }
              //  .view { "Files found from search for cluster input: $it" }
               // .flatten()
               // .map { file -> [ [id: file.baseName], file ] }
            
            CLUSTER_DMND_MCL(
                SEARCH.out.domfasta,
                params.cluster_dmnd_mcl.dmnd_params,
                params.cluster_dmnd_mcl.mcl_params,
                params.cluster_dmnd_mcl.mcl_inflation,
                "orthofinder"
                )
        }
        
        }
    }

}
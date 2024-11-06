#!/usr/bin/env nextflow


// Include modules
include { BAM_INDEX } from './subworkflows/bam_index'
include { INITIAL_ORTHOLOGY } from './subworkflows/initial_orthology'
include { SEARCH } from './modules/local/search/search'

workflow {
    if (params.run_bam_index) {
        BAM_INDEX(params.sample_sheet)
    }
// Orthofinder and search
    if (params.run_orthofinder) {
        INITIAL_ORTHOLOGY(
            params.fasta_dir,
            params.prior_run
        )

        if (params.run_search) {
            // Wait for INITIAL_ORTHOLOGY to finish and use its output
            // Define channels
            ch_fasta = INITIAL_ORTHOLOGY.out.orthogroup_sequences
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
        }
    }

}
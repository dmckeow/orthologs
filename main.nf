#!/usr/bin/env nextflow

// Include modules
include { WF_ORTHOFINDER } from './subworkflows/initial_orthology'
include { WF_BROCCOLI } from './subworkflows/initial_orthology.nf'
include { SEARCH as SEARCH_ORTHOFINDER } from './modules/local/search/search'
include { SEARCH as SEARCH_BROCCOLI } from './modules/local/search/search'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_ORTHOFINDER } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_BROCCOLI } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'


workflow {
    // Initialize variables
def orthofinder_msg = "🔍 ❌ Will skip OrthoFinder"
def broccoli_msg = "🥦 ❌ Will skip Broccoli"
def search_msg = "🔨 ❌ Will skip HMMSEARCH on orthogroups"
def cluster_msg = "💎🍇 ❌ Will skip DIAMOND and MCL clustering of orthogroups"

if (!params.run.orthofinder && !params.run.broccoli) {
    log.error """
    ============================================================
    Error: Neither OrthoFinder nor Broccoli is set to run!
    
    You must set at least one of the following to true:
    - params.run.orthofinder
    - params.run.broccoli
    
    Exiting the pipeline.
    ============================================================
    """
    exit 1
}

// Orthofinder and search
if (params.run.orthofinder) {
    orthofinder_msg = "🔍 ✅ Will run OrthoFinder, then:"
}

if (params.run.broccoli) {
    broccoli_msg = "🥦 ✅ Will run Broccoli, then:"
}

    
if (params.run.search) {
    search_msg = "🔨 ✅ Will run HMMSEARCH on orthogroups, then:"
    
    if (params.run.cluster_dmnd_mcl) {
        cluster_msg = "💎🍇 ✅ Will run DIAMOND and MCL clustering of orthogroups, ❗ but ONLY on the orthogroups with HMMSEARCH hits ❗"
    }
} else if (params.run.cluster_dmnd_mcl) {
    cluster_msg = "💎🍇 ✅ Will run DIAMOND and MCL clustering of orthogroups, ❗ but on ALL orthogroups identified ❗"
}

// Print all messages as a single unit
log.info """
Pipeline workflow that will be executed:
---------------------------
    ${broccoli_msg}
        ${search_msg}
            ${cluster_msg}

    ${orthofinder_msg}
        ${search_msg}
            ${cluster_msg}
    """

// Orthofinder and search
    if (params.run.orthofinder) {
        WF_ORTHOFINDER(
            params.fasta_dir,
            params.orthofinder.prior_run
        )
        // Set the channel from the orthofinder output
        ch_orthofinder_fastas = WF_ORTHOFINDER.out.orthogroup_sequences
            .map { meta, dir -> file("${dir}/*.fa") }
            .flatten()
            .map { file -> [ [id: file.baseName], file ] }

        if (params.run.search) {
            // Run SEARCH_ORTHOFINDER process
            SEARCH_ORTHOFINDER(
                ch_orthofinder_fastas,
                params.search.gene_family_info,
                params.search.gene_family_name,
                file(params.search.hmm_dir),
                "orthofinder",
                params.outdir
            )
        if (params.run.cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_ORTHOFINDER(
                SEARCH_ORTHOFINDER.out.domfasta,
                params.cluster.dmnd.args,
                params.cluster.mcl.args,
                params.cluster.mcl.inflation,
                "searches/orthofinder"
                )
            }
        } else {
            if (params.run.cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_ORTHOFINDER(
                ch_orthofinder_fastas,
                params.cluster.dmnd.args,
                params.cluster.mcl.args,
                params.cluster.mcl.inflation,
                "all/orthofinder"
                )
            }
        }
    }

// Broccoli and search
    if (params.run.broccoli) {
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

        if (params.run.search) {
            // Run SEARCH process
            SEARCH_BROCCOLI(
                ch_broccoli_fastas,
                params.search.gene_family_info,
                params.search.gene_family_name,
                file(params.search.hmm_dir),
                "broccoli",
                params.outdir
            )
            if (params.run.cluster_dmnd_mcl) {
                CLUSTER_DMND_MCL_BROCCOLI(
                    SEARCH_BROCCOLI.out.domfasta,
                    params.cluster.dmnd.args,
                    params.cluster.mcl.args,
                    params.cluster.mcl.inflation,
                    "searches/broccoli"
                    )
                }
        } else {
            if (params.run.cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_BROCCOLI(
                ch_broccoli_fastas,
                params.cluster.dmnd.args,
                params.cluster.mcl.args,
                params.cluster.mcl.inflation,
                "all/broccoli"
                )
            }
         }
    }
}

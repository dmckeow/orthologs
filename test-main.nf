#!/usr/bin/env nextflow

include { INITIAL_ORTHOGROUPS } from './subworkflows/initial_orthogroups'
//include { WF_CLUSTERING } from './subworkflows/wf_clustering'

workflow {
    // Initialize variables for messages
    def orthofinder_msg      = "🔍 OrthoFinder: SKIP ❌"
    def broccoli_msg         = "🥦 Broccoli:    SKIP ❌"
    def search_msg           = "🔨 HMMSEARCH:   SKIP ❌"
    def cluster_dmnd_mcl_msg = "💎 DIAMOND:     SKIP ❌"
    def cluster_mmseqs_msg   = "🚀 MMseqs2:     SKIP ❌"
    def vs_msg               = "⚠️  Clustering ALL orthogroups  ⚠️"

    // Check if at least one of OrthoFinder or Broccoli is set to run
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

    // Update messages based on parameters
    if (params.run.orthofinder) orthofinder_msg     = "🔍 OrthoFinder: RUN ✅"
    if (params.run.broccoli)    broccoli_msg        = "🥦 Broccoli:    RUN ✅"
    if (params.run.search) {
        vs_msg = "⚠️  Clustering ONLY HMMSEARCH results  ⚠️"
        search_msg           = "🔨 HMMSEARCH:  RUN ✅"
    }
    if (params.run.cluster_dmnd_mcl) cluster_dmnd_mcl_msg = "💎 DIAMOND:    RUN ✅"
    if (params.run.cluster_mmseqs)   cluster_mmseqs_msg   = "🚀 MMseqs2:    RUN ✅"

    // Print workflow messages
    log.info """
    Pipeline workflow that will be executed:
    ---------------------------
    🥇 Initial orthology:
        ${broccoli_msg}
        ${orthofinder_msg}

    👪 Gene family search of orthologs (SEARCH):
        ${search_msg}

    🍇 Clustering:
        ${vs_msg}
        ${cluster_dmnd_mcl_msg}
        ${cluster_mmseqs_msg}
    """

    // Run the orthology workflow
    INITIAL_ORTHOGROUPS(
        params.fasta_dir,
        params.orthofinder.prior_run,
        params.orthofinder.min_sequences,
        params.broccoli.args,
        params.broccoli.min_sequences,
        params.search.gene_family_info,
        params.search.gene_family_name,
        params.search.hmm_dir,
        params.outdir,
        params.run.orthofinder,
        params.run.broccoli,
        params.run.search
    )

    // Run the clustering workflow
    //WF_CLUSTERING(
    //    WF_ORTHOLOGY.out.search_results,
    //    params.cluster.dmnd.args,
    //    params.cluster.mcl.args,
    //   params.cluster.mcl.inflation,
    //    params.run.cluster_dmnd_mcl,
    //    params.run.cluster_mmseqs,
    //    params.run.search ? "searches" : "all"
    //)
}
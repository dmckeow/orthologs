#!/usr/bin/env nextflow

include { INITIAL_ORTHOGROUPS } from './subworkflows/initial_orthogroups'
//include { WF_CLUSTERING } from './subworkflows/wf_clustering'

workflow {
    // Initialize variables for messages
    def search_msg           = "🔨 HMMSEARCH:     SKIP ❌ - ⚠️  Without this step, instead of identifying orthologs in gene families identified by HMMSEARCH, the pipeline will try to identify orthologs for all input sequences  ⚠️"
    def orthofinder_msg      = "🔍 OrthoFinder:   SKIP ❌"
    def broccoli_msg         = "🥦 Broccoli:      SKIP ❌"
    def cluster_dmnd_mcl_msg = "💎 DIAMOND + MCL: SKIP ❌"
    def cluster_mmseqs_msg   = "🚀 MMseqs2:       SKIP ❌"

    // Check if at least one of OrthoFinder or Broccoli is set to run
    if (!params.run.orthofinder && 
        !params.run.broccoli &&
        !params.run.cluster_dmnd_mcl &&
        !params.run.cluster_mmseqs
        ) {
        log.error """
        ============================================================
        Error: No initial orthogroup method is set to run! Interesting choice...
        
        You must set at least one of the following to true:
        - params.run.orthofinder
        - params.run.broccoli
        - params.run.cluster_dmnd_mcl
        - params.run.cluster_mmseqs
        
        Exiting the pipeline.
        ============================================================
        """
        exit 1
    }

    // Update messages based on parameters
    if (params.run.search) {
        search_msg                                        = "🔨 HMMSEARCH:     RUN ✅ - will only identify orthologs in gene families identified by HMMSEARCH"
    }

    if (params.run.orthofinder) orthofinder_msg           = "🔍 OrthoFinder:   RUN ✅"
    if (params.run.broccoli)    broccoli_msg              = "🥦 Broccoli:      RUN ✅"
    if (params.run.cluster_dmnd_mcl) cluster_dmnd_mcl_msg = "💎 DIAMOND + MCL: RUN ✅"
    if (params.run.cluster_mmseqs)   cluster_mmseqs_msg   = "🚀 MMseqs2:       RUN ✅"

    // Print workflow messages
    log.info """
    Pipeline workflow that will be executed:
    ---------------------------
    🔨 Gene family search of input sequences:
        ${search_msg}

    🍇 Processes to identify orthogroups:
        ${broccoli_msg}
        ${orthofinder_msg}
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
        params.run.search,
        params.cluster.dmnd.args,
        params.cluster.mcl.args,
        params.cluster.mcl.inflation,
        params.run.cluster_dmnd_mcl,
        params.run.cluster_mmseqs
    )
}
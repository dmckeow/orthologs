#!/usr/bin/env nextflow

include { INIT_ORTHO } from './subworkflows/init_ortho'
include { WHICH_ORTHO } from './subworkflows/which_ortho'
import groovy.json.JsonOutput

def printParams() {
        def output = [:]
        params.each { k, v ->
            if (v instanceof Map) {
                output[k] = v
            } else {
                output[k] = v.toString()
            }
        }
        log.info """
        =========================================
        ${workflow.manifest.name} v${workflow.manifest.version}
        =========================================
        Parameters:
        ${JsonOutput.prettyPrint(JsonOutput.toJson(output))}
        =========================================
        """
    }

workflow {

    // Initialize variables for messages
    def search_msg           = "🔨 HMMSEARCH:     SKIP ❌ - ⚠️  Without this step, instead of identifying orthologs in gene families identified by HMMSEARCH, the pipeline will try to identify orthologs for all input sequences  ⚠️"
    def orthofinder_msg      = "🔍 OrthoFinder:   SKIP ❌"
    def broccoli_msg         = "🥦 Broccoli:      SKIP ❌"
    def cluster_dmnd_mcl_msg = "💎 DIAMOND + MCL: SKIP ❌"
    def cluster_mmseqs_msg   = "🚀 MMseqs2:       SKIP ❌"

    // Print workflow paramters in startup
    printParams()

    // Check if at least one of OrthoFinder or Broccoli is set to run
    if (!params.run.orthofinder && 
        !params.run.broccoli &&
        !params.run.cluster_dmnd_mcl &&
        !params.run.cluster_mmseqs
        ) {
        log.error """
        ============================================================
        Error: No initial orthogroup method is set to run! I am afraid you must choose at least one...
        
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

    params.samplesheet = null

    if (params.samplesheet == null) {
        error "Please provide a samplesheet with the '--samplesheet' option."
    }

    // Run the orthology workflow
    if (params.run.init_ortho) {
        INIT_ORTHO (
            params.samplesheet,
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

    // Gather orthogroup info to generate reports and assess which OGs are the best for phylogeny
    if (params.run.which_ortho) {
        WHICH_ORTHO (
            INIT_ORTHO.out.combined_deflines,
            INIT_ORTHO.out.ch_orthofinder_out,
            INIT_ORTHO.out.ch_broccoli_og_table,
            INIT_ORTHO.out.ch_dmnd_mcl_og_table,
            INIT_ORTHO.out.ch_mmseqs_og_table,
            params.run.search,
            INIT_ORTHO.out.ch_orthofinder_og_fa_dir,
            INIT_ORTHO.out.ch_broccoli_og_fa_dir,
            INIT_ORTHO.out.ch_dmnd_mcl_og_fa_dir,
            INIT_ORTHO.out.ch_mmseqs_og_fa_dir,
            INIT_ORTHO.out.combined_fasta,
            params.interproscan.db
        )
    }
}
#!/usr/bin/env nextflow

// Include modules
include { WF_ORTHOFINDER } from './subworkflows/initial_orthology'
include { WF_BROCCOLI } from './subworkflows/initial_orthology.nf'
include { SEARCH as SEARCH_ORTHOFINDER } from './modules/local/search/search'
include { SEARCH as SEARCH_BROCCOLI } from './modules/local/search/search'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_ORTHOFINDER } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { CLUSTER_DMND_MCL as CLUSTER_DMND_MCL_BROCCOLI } from './modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { WF_CLUSTER_MMSEQS as WF_CLUSTER_MMSEQS_BROCCOLI } from './subworkflows/cluster_mmseqs'
include { WF_CLUSTER_MMSEQS as WF_CLUSTER_MMSEQS_ORTHOFINDER } from './subworkflows/cluster_mmseqs'


workflow {
    // Initialize variables
def orthofinder_msg      = "🔍 OrthoFinder: SKIP ❌"
def broccoli_msg         = "🥦 Broccoli:    SKIP ❌"
def search_msg           = "🔨 HMMSEARCH:   SKIP ❌"
def cluster_dmnd_mcl_msg = "💎 DIAMOND:     SKIP ❌"
def cluster_mmseqs_msg   = "🚀 MMseqs2:     SKIP ❌"
def vs_msg   = "⚠️  Clustering ALL orthogroups  ⚠️"

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
    orthofinder_msg     = "🔍 OrthoFinder: RUN ✅"
}

if (params.run.broccoli) {
    broccoli_msg        = "🥦 Broccoli:    RUN ✅"
}
    
if (params.run.search) {
    vs_msg = "⚠️  Clustering ONLY HMMSEARCH results  ⚠️"
    search_msg           = "🔨 HMMSEARCH:  RUN ✅"
}
if (params.run.cluster_dmnd_mcl) {
    cluster_dmnd_mcl_msg = "💎 DIAMOND:    RUN ✅"
}
if (params.run.cluster_mmseqs) {
    cluster_mmseqs_msg   = "🚀 MMseqs2:    RUN ✅"
}


// Print all messages as a single unit
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

// Orthofinder and search
    if (params.run.orthofinder) {
        WF_ORTHOFINDER(
            params.fasta_dir,
            params.orthofinder.prior_run,
            params.orthofinder.min_sequences
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
        SEARCH_ORTHOFINDER.out.domfasta.view()
        if (params.run.cluster_dmnd_mcl) {
            CLUSTER_DMND_MCL_ORTHOFINDER(
                SEARCH_ORTHOFINDER.out.domfasta,
                params.cluster.dmnd.args,
                params.cluster.mcl.args,
                params.cluster.mcl.inflation,
                "searches/orthofinder"
                )
            }
        if (params.run.cluster_mmseqs) {
            WF_CLUSTER_MMSEQS_ORTHOFINDER(
                SEARCH_ORTHOFINDER.out.domfasta,
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
            if (params.run.cluster_mmseqs) {
            WF_CLUSTER_MMSEQS_ORTHOFINDER(
                ch_orthofinder_fastas,
                "all/orthofinder"
                )
            }
        }
    }

// Broccoli and search
    if (params.run.broccoli) {
        WF_BROCCOLI(
            params.fasta_dir,
            params.broccoli.args,
            params.broccoli.min_sequences
        )

        // Set the channel from the orthofinder output
        //ch_broccoli_fastas = WF_BROCCOLI.out.orthologous_groups_sequences
        //.map { dir -> 
        //    def files = file("${dir}/*.fa")
        //    files.collect { file -> [ [id: file.baseName], file ] }
        //}
        //.flatten()
        //.collate(2)

        // Set the channel from the broccoli output
        ch_broccoli_fastas = WF_BROCCOLI.out.orthologous_groups_sequences
            .map { meta, dir -> file("${dir}/*.fa") }
            .flatten()
            .map { file -> [ [id: file.baseName], file ] }

        //ch_broccoli_fastas.view()

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
            SEARCH_BROCCOLI.out.domfasta.view()
            ch_search_broccoli_results = SEARCH_BROCCOLI.out.domfasta

            if (params.run.cluster_dmnd_mcl) {
                ch_search_broccoli_results
                    .branch {
                        empty: it.isEmpty()
                        data: true
                    }
                    .set { branched_results }

                branched_results.empty
                    .ifEmpty { /* do nothing */ }
                    .subscribe { log.warn "No SEARCH hits vs Broccoli orthogroups. Skipping CLUSTER_DMND_MCL_BROCCOLI." }


                CLUSTER_DMND_MCL_BROCCOLI(
                    branched_results.data,
                    params.cluster.dmnd.args,
                    params.cluster.mcl.args,
                    params.cluster.mcl.inflation,
                    "searches/broccoli"
                    )
                }
            if (params.run.cluster_mmseqs) {
                ch_search_broccoli_results
                    .branch {
                        empty: it.isEmpty()
                        data: true
                    }
                    .set { branched_results_mmseqs }

                branched_results_mmseqs.empty
                    .ifEmpty { /* do nothing */ }
                    .subscribe { log.warn "No SEARCH hits vs Broccoli orthogroups. Skipping WF_CLUSTER_MMSEQS_BROCCOLI and PARSE_MMSEQS_TO_FASTA." }

                WF_CLUSTER_MMSEQS_BROCCOLI(
                    branched_results_mmseqs.data,
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

            if (params.run.cluster_mmseqs) {
                WF_CLUSTER_MMSEQS_BROCCOLI(
                    ch_broccoli_fastas,
                    "all/broccoli"
                    )
                }
         }
    }
}

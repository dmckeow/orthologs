#!/usr/bin/env nextflow

include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main' //
include { BROCCOLI } from '../modules/local/broccoli/main'
include { SEARCH } from '../modules/local/search/search'


// DMND MMSEQS

include { CLUSTER_DMND_MCL } from '../modules/local/cluster_dmnd_mcl/cluster_dmnd_mcl'
include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER } from '../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { PARSE_MMSEQS_TO_FASTA } from '../modules/local/cluster_mmseqs/parse_mmseqs'

include { ORTHOFINDER as ORTHOFINDER_INITIAL } from '../modules/nf-core/orthofinder/main' // OVLP

process CONCATENATE_FASTAS {
    input:
    tuple val(meta), path(fasta_files)

    output:
    tuple val(meta), path("combined_sequences.fasta"), emit: combined_fasta

    script:
    """
    cat ${fasta_files} | sed 's/\\*/X/g' > combined_sequences.fasta
    # asterisks removed because of interproscan later
    """
}

// Processes to gather fasta names, sample name per fasta defline for later mapping
process EXTRACT_DEFLINES {
    
    input:
    tuple val(meta), path(fasta)

    output:
    path "${meta.id}_deflines.txt"

    script:
    """
    #!/bin/bash
    while read -r line; do
        if [[ \$line == ">"* ]]; then
            defline=\${line#>}
            echo "${meta.id}\t${fasta.getName()}\t\$defline" >> ${meta.id}_deflines.txt
        fi
    done < $fasta
    """
}

process COMBINE_DEFLINES {
    publishDir "${params.outdir}/deflines", mode: 'copy'

    input:
    path '*_deflines.txt'

    output:
    path 'deflines_combined.txt', emit: combined_deflines

    script:
    """
    #!/bin/bash
    echo "sample\tfasta\tdefline" > deflines_combined.txt
    cat *_deflines.txt >> deflines_combined.txt
    """
}


// OVLP

process ID_NO_OVERLAPS {
    input:
    path overlap_matrix
    
    output:
    path "no_overlaps.txt", emit: no_overlaps
    
    script:
    """
    python ${projectDir}/bin/id_orthofinder_no_overlap.py -i $overlap_matrix
    """
}

process FILTER_INPUTS_ORTHOFINDER {
    input:
    tuple val(meta), path(files)
    path no_overlaps

    output:
    tuple val(meta), path('filtered_files/*'), emit: filtered_input
    path 'filter_summary.txt', emit: summary

    script:
    """
    echo "FILTER_INPUTS_ORTHOFINDER:\nThe input fastas were assessed by their DIAMOND all-vs-all overlaps," > filter_summary.txt
    echo "to identify genomes sharing 0 overlaps with any other genome," >> filter_summary.txt
    echo "When two genomes shared an overlap of 0, the genomes with less total overlaps with all genomes was marked for removal" >> filter_summary.txt
    echo "Number of input files before filtering: \$(ls ${files} | wc -l)" >> filter_summary.txt
    echo "Genomes filtered out from OrthoFinder input:" >> filter_summary.txt

    mkdir filtered_files
    filtered_out=0
    for file in ${files}; do
        basename=\$(basename \$file .fasta)
        basename=\${basename%.domains}
        if grep -q "\$basename" $no_overlaps; then
            echo "\$file" >> filter_summary.txt
            filtered_out=\$((filtered_out+1))
        else
            cp \$file filtered_files/
        fi
    done

    echo "" >> filter_summary.txt
    echo "Number of files filtered out: \$filtered_out" >> filter_summary.txt
    echo "Number of files after filtering: \$(ls filtered_files | wc -l)" >> filter_summary.txt

    """
}


workflow INIT_ORTHO {
    take:
    samplesheet
    orthofinder_prior_run
    orthofinder_min_sequences
    broccoli_args
    broccoli_min_sequences
    search_gene_family_info
    search_gene_family_name
    search_hmm_dir
    outdir
    run_orthofinder
    run_broccoli
    run_search
    cluster_dmnd_args
    cluster_mcl_args
    cluster_mcl_inflation
    run_cluster_dmnd_mcl
    run_cluster_mmseqs

    main:
    ch_versions = Channel.empty()
    ch_orthofinder_out = Channel.empty() // Parent directory for output
    ch_orthofinder_og_fa_dir = Channel.empty() // Parent directory containing orthogroup fasta files
    ch_broccoli_og_fa_dir = Channel.empty()
    ch_input_fastas = Channel.empty()
    ch_dmnd_mcl_og_fa_list = Channel.empty()
    ch_mmseqs_og_fa_dir = Channel.empty()
    ch_broccoli_og_table = Channel.empty()
    ch_dmnd_mcl_og_table = Channel.empty()
    ch_mmseqs_og_table = Channel.empty()

    // Create a channel for the FASTA files
    //ch_fasta_files = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
        //.map { file -> [ [id: file.baseName], file ] }

    //ch_fasta_files.view { it -> "ch_fasta_files: $it" }


    // SASH: (replaced part above)
    ch_fasta_files = Channel.fromPath(samplesheet)
        .splitCsv(header:true, sep:',')
        .map { row -> 
            def meta = [id: row.id]
            def fasta = file(row.fasta)
            return [meta, fasta]
        }

    // Search

    if (run_search) {

        SEARCH(
            ch_fasta_files,
            search_gene_family_info,
            search_gene_family_name,
            file(search_hmm_dir),
            "initial",
            outdir
        )
        ch_input_fastas = SEARCH.out.domfasta

    } else {
        //ch_input_fastas = Channel.fromPath("${fasta_dir}/*.{fa,faa,fasta,fas,pep}")
            //.map { file -> [ [id: file.baseName], file ] }
        ch_input_fastas = ch_fasta_files // SASH, replaces part above
    }

    //ch_input_fastas.view { it -> "ch_input_fastas: $it\nExpected example: [[id:sample], /path/fasta]" }

    EXTRACT_DEFLINES(ch_input_fastas)
    COMBINE_DEFLINES(EXTRACT_DEFLINES.out.collect())

// Basic OrthoFINDER:
//    if (run_orthofinder) {
        // Collect all input fastas for OrthoFinder
//        ch_orthofinder_input = ch_input_fastas
//            .map { meta, file -> file }  // Extract just the file paths
//            .collect()  // Collect all files into a single list
//            .map { files -> 
//                [[id: "orthofinder"], files]  // Create the structure OrthoFinder expects
//            }

        // Create a channel for the optional prior run
//        prior_run_ch = orthofinder_prior_run
//            ? Channel.value([[id: 'prior_run'], file(orthofinder_prior_run)])
//            : Channel.value([[],[]])

//        ORTHOFINDER(ch_orthofinder_input, prior_run_ch)
//        ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)
//        ch_orthofinder_out = ORTHOFINDER.out.orthofinder

        // Create a channel for the Orthogroup_Sequences folder
//        ch_orthofinder_og_fa_dir = ORTHOFINDER.out.orthofinder
//            .map { meta, path -> 
//                def orthogroup_dir = file("${path}/Orthogroup_Sequences")
//                return [meta, orthogroup_dir]
//            }
//    }

// Orthofinder OVLP:
    if (run_orthofinder) {
        // Run initial OrthoFinder step

        // Collect all input fastas for OrthoFinder
                ch_orthofinder_input = ch_input_fastas
                    .map { meta, file -> file }  // Extract just the file paths
                    .collect()  // Collect all files into a single list
                    .map { files -> 
                        [[id: "orthofinder"], files]  // Create the structure OrthoFinder expects
                    }

                // Create a channel for the optional prior run
                prior_run_ch = orthofinder_prior_run
                    ? Channel.value([[id: 'prior_run'], file(orthofinder_prior_run)])
                    : Channel.value([[],[]])

        ORTHOFINDER_INITIAL(ch_orthofinder_input, prior_run_ch)
        
        // Extract the SpeciesOverlap.tsv file
        ch_overlap_matrix = ORTHOFINDER_INITIAL.out.orthofinder
            .map { meta, path -> 
                def overlap_file = file("${path}/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv")
                return overlap_file
            }

        ch_overlap_matrix.view { it -> "OrthoFinder's DIAMOND all-vs-all ortholog overlap matrix is here: $it" }
        
ID_NO_OVERLAPS(ch_overlap_matrix)

// Debug: View the contents of ch_orthofinder_input
//ch_orthofinder_input.view { it -> "ch_orthofinder_input before filtering: $it" }

// Run the filtering process
    FILTER_INPUTS_ORTHOFINDER(ch_orthofinder_input, ID_NO_OVERLAPS.out.no_overlaps)

    FILTER_INPUTS_ORTHOFINDER.out.summary.view()

    // Use the filtered output
    ch_filtered_input = FILTER_INPUTS_ORTHOFINDER.out.filtered_input

    // Debug: View the contents of ch_filtered_input
    ch_filtered_input.view { meta, files -> 
        "ch_filtered_input after filtering: meta=$meta, files=${files.size()}"
    }

// Run final OrthoFinder with filtered input
ORTHOFINDER(ch_filtered_input, prior_run_ch)
        
        ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)
        ch_orthofinder_out = ORTHOFINDER.out.orthofinder
        
        // Create a channel for the Orthogroup_Sequences folder
        ch_orthofinder_og_fa_dir = ORTHOFINDER.out.orthofinder
            .map { meta, path -> 
                def orthogroup_dir = file("${path}/Orthogroup_Sequences")
                return [meta, orthogroup_dir]
            }
    }


    if (run_broccoli) {
        // Collect all input fastas for Broccoli
        ch_broccoli_input = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> 
                [[id: "broccoli"], files]  // Create the structure OrthoFinder expects
            }

        BROCCOLI(
            ch_broccoli_input,
            broccoli_args
        )
        ch_versions = ch_versions.mix(BROCCOLI.out.versions)
        ch_broccoli_og_table = BROCCOLI.out.table_OGs_protein_names
        
        // Add metadata to the orthologous_groups_sequences output
        ch_broccoli_og_fa_dir = BROCCOLI.out.orthologous_groups_sequences
            .map { path -> 
                [[id: "broccoli_orthogroups"], path]
            }

        // Filter orthogroups by min number sequences
        //FILTER_ORTHOGROUPS_BROCCOLI(ch_broccoli_og_fa_dir, broccoli_min_sequences)
        //ch_broccoli_results = FILTER_ORTHOGROUPS_BROCCOLI.out.filtered_orthogroups

        //ch_broccoli_results.view { it -> "ch_broccoli_results: $it" }
    }

    if (run_cluster_dmnd_mcl) {
        // Group all input fastas for CLUSTER_DMND_MCL
        ch_dmnd_mcl_input = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> [[id: "combined"], files] }  // Add metadata

        CLUSTER_DMND_MCL(
            ch_dmnd_mcl_input, // multi input
            cluster_dmnd_args,
            cluster_mcl_args,
            cluster_mcl_inflation,
            "mcl"
        )

        ch_dmnd_mcl_og_fa_dir = CLUSTER_DMND_MCL.out.fasta_dir
        ch_dmnd_mcl_og_table = CLUSTER_DMND_MCL.out.orthogroup_list
    }

    if (run_cluster_mmseqs) {
        
        // Collect all input fastas for concatenation
        ch_fasta_for_concat = ch_input_fastas
            .map { meta, file -> file }  // Extract just the file paths
            .collect()  // Collect all files into a single list
            .map { files -> 
                [[id: "combined"], files]
            }
        //ch_fasta_for_concat.view { it -> "ch_fasta_for_concat: $it" }

        // Concatenate all FASTA files
        CONCATENATE_FASTAS(ch_fasta_for_concat)

        //CONCATENATE_FASTAS.out.combined_fasta.view { it -> "CONCATENATE_FASTAS.out.combined_fasta: $it" }

        // Create MMseqs2 database
        MMSEQS_CREATEDB(CONCATENATE_FASTAS.out.combined_fasta)

        // Cluster sequences
        MMSEQS_CLUSTER(MMSEQS_CREATEDB.out.db)

        //MMSEQS_CLUSTER.out.db_cluster.view { "MMSEQS_CLUSTER output: $it" }
        //MMSEQS_CREATEDB.out.db.view { "MMSEQS_CREATEDB output: $it" }

        PARSE_MMSEQS_TO_FASTA(
            MMSEQS_CREATEDB.out.db,
            MMSEQS_CLUSTER.out.db_cluster,
            CONCATENATE_FASTAS.out.combined_fasta,
            "mmseqs"
        )
        
        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)
        ch_versions = ch_versions.mix(MMSEQS_CLUSTER.out.versions)
        
        ch_mmseqs_og_fa_dir = PARSE_MMSEQS_TO_FASTA.out.fasta_dir
        ch_mmseqs_og_table = PARSE_MMSEQS_TO_FASTA.out.orthogroup_list

    }

        //ch_orthofinder_og_fa_dir.view { it -> "ch_orthofinder_og_fa_dir: $it" }
        //ch_broccoli_og_fa_dir.view { it -> "ch_broccoli_og_fa_dir: $it" }
        //ch_dmnd_mcl_og_fa_dir.view { it -> "ch_dmnd_mcl_og_fa_dir: $it" }
        //ch_mmseqs_og_fa_dir.view { it -> "ch_mmseqs_og_fa_dir: $it" }
        
    emit:
        versions = ch_versions
        input_fastas = ch_input_fastas
        combined_deflines = COMBINE_DEFLINES.out.combined_deflines
        combined_fasta = CONCATENATE_FASTAS.out.combined_fasta
        ch_orthofinder_og_fa_dir = ch_orthofinder_og_fa_dir
        ch_orthofinder_out = ch_orthofinder_out
        ch_broccoli_og_fa_dir = ch_broccoli_og_fa_dir
        ch_broccoli_og_table = ch_broccoli_og_table
        ch_dmnd_mcl_og_fa_dir = ch_dmnd_mcl_og_fa_dir
        ch_dmnd_mcl_og_table = ch_dmnd_mcl_og_table
        ch_mmseqs_og_fa_dir = ch_mmseqs_og_fa_dir
        ch_mmseqs_og_table = ch_mmseqs_og_table

}
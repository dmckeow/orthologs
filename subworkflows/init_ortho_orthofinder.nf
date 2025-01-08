#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// SINGLE USE MODULES FROM ORIGINAL INIT_ORTHO 
include { SEARCH } from '../modules/local/search/search'
include { ORTHOFINDER } from '../modules/nf-core/orthofinder/main'


// MULTI USE MODULES FROM ORIGINAL INIT_ORTHO 
include { ORTHOFINDER as ORTHOFINDER_INITIAL } from '../modules/nf-core/orthofinder/main'
include { PARSE_FASTAS as ORTHOFINDER_PARSE_FASTAS } from '../modules/local/parse_fastas/main'

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
    tuple val(meta), path(orthofinder_output, stageAs: "orthofinder_parent/*")
    
    output:
    path "no_overlaps.txt", emit: no_overlaps
    
    script:
    """
    python ${projectDir}/bin/id_orthofinder_no_overlap.py \
        -l orthofinder_parent \
        -o no_overlaps.txt
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
    echo "FILTER_INPUTS_ORTHOFINDER:\nWhenever OrthoFinder flagged pairs of genomes not sharing enough overlaps in DIAMOND all-vs-all, the problematic genomes were removed," > filter_summary.txt
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


workflow INIT_ORTHO_ORTHOFINDER {
    take:
    samplesheet
    orthofinder_prior_run
    orthofinder_min_sequences
    search_gene_family_info
    search_gene_family_name
    search_hmm_dir
    outdir
    run_search
    
    main:
    ch_versions = Channel.empty()
    ch_orthofinder_out = Channel.empty() // Parent directory for output
    ch_orthofinder_og_fa_dir = Channel.empty() // Parent directory containing orthogroup fasta files
    ch_input_fastas = Channel.empty()
    

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
        



ID_NO_OVERLAPS(ORTHOFINDER_INITIAL.out.orthofinder.map { meta, path -> [meta, path.parent] })

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
        // Create a channel for the Phylogenetic_Hierarchical_Orthogroups folder
        ch_orthofinder_hog_tsv = ORTHOFINDER.out.orthofinder
            .map { meta, path -> 
                def orthogroup_dir = file("${path}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
                return [meta, orthogroup_dir]
            }
        // Parse fastas for the hierarchical OGs of orthofinder
        ORTHOFINDER_PARSE_FASTAS(
            ch_filtered_input,
            ch_orthofinder_hog_tsv,
            'Phylogenetic_Hierarchical_Orthogroups_Sequences'
        )
    }


    

    

    
    emit:
        versions = ch_versions
        input_fastas = ch_input_fastas
        combined_deflines = COMBINE_DEFLINES.out.combined_deflines
        ch_orthofinder_og_fa_dir = ch_orthofinder_og_fa_dir
        ch_orthofinder_out = ch_orthofinder_out
}
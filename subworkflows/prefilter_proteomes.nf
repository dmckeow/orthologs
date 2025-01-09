#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// SINGLE USE MODULES FROM ORIGINAL INIT_ORTHO 
include { SEARCH } from '../modules/local/search/search'

workflow PREFILTER_PROTEOMES {
    take:
    samplesheet
    search_params
    outdir

    main:
    // Read the original sample sheet
    Channel
        .fromPath(samplesheet)
        .splitCsv(header:true)
        .map { row -> tuple([id: row.id], file(row.fasta)) }
        .set { ch_input_samples }

    // Read the new search parameters sheet
    Channel
        .fromPath(search_params)
        .splitCsv(header:true)
        .set { ch_search_params }

    // Combine input samples with each set of search parameters
    ch_search_inputs = ch_input_samples.combine(ch_search_params)

    // Run SEARCH for each combination
    SEARCH(
        ch_search_inputs.map { sample, fasta, params -> tuple(sample, fasta) },
        ch_search_inputs.map { sample, fasta, params -> params.gene_family_info },
        ch_search_inputs.map { sample, fasta, params -> params.gene_family_name },
        ch_search_inputs.map { sample, fasta, params -> file(params.hmm_dir) },
        Channel.value("initial"),
        Channel.value(outdir)
    )

    // Gather the meta info for the gene famiyl search results, and the fastas too
    defline_info_collected = SEARCH.out.defline_info
        .collectFile(name: 'all_defline_info.csv', keepHeader: false, skip: 1)

    // Process the collected information into a metamap
    defline_metamap = defline_info_collected
        .splitCsv(header: ['defline', 'fasta', 'gene_family_name'])
        .map { row -> 
            [
                defline: row.defline,
                fasta: row.fasta,
                gene_family_name: row.gene_family_name
            ]
        }
        .collect()
        .map { list -> [metamap: list] }
    
    domfastas_collected = SEARCH.out.domfasta.collect()

    defline_metamap.view { it -> "defline_metamap: $it" }
    domfastas_collected.view { it -> "domfastas_collected: $it" }

    emit:
    metamap = defline_metamap
    fastas  = domfastas_collected
}
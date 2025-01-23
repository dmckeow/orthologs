#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.preview.output = true

// SINGLE USE MODULES FROM ORIGINAL INIT_ORTHO 
include { PREFILTER_SEARCH } from '../modules/local/prefilter_search'

process CONCATENATE_FASTAS {
    tag "Concatentaing fastas from hmmsearch to keep one input per species or OrthoFinder etc"
    label 'process_low'
    publishDir "${params.outdir}/prefilter/initial/clean_fasta/merged", mode: 'copy'

    input:
    tuple val(id), path(fastas)

    output:
    tuple val(id), path("${id}.fasta"), emit: concatenated_fasta

    script:
    """
    cat ${fastas.join(' ')} > ${id}.fasta
    """
}




workflow PREFILTER {
    take:
    samplesheet
    search_params
    outdir
    run_prefilter_hmmsearch

    main:
        // Define the input samples channel
    ch_input_samples = Channel
        .fromPath(samplesheet)
        .splitCsv(header:true)
        .map { row -> tuple([id: row.id], file(row.fasta)) }

    // Define the search parameters channel, which will be empty if search_params is null
    ch_search_params = run_prefilter_hmmsearch ? Channel
        .fromPath(params.search_params)
        .splitCsv(header:true)
        : Channel.value([gene_family_info: 'no_hmmsearch', gene_family_name: 'no_hmmsearch', hmm_dir: 'no_hmmsearch'])

    // Combine input samples with search parameters, or just use input samples if search_params is null
    ch_search_inputs = ch_input_samples.combine(ch_search_params)

    // Run PREFILTER_SEARCH
    PREFILTER_SEARCH(
        ch_search_inputs.map { sample, fasta, params -> tuple(sample, fasta) },
        ch_search_inputs.map { sample, fasta, params -> params.gene_family_info },
        ch_search_inputs.map { sample, fasta, params -> params.gene_family_name },
        ch_search_inputs.map { sample, fasta, params -> params.hmm_dir ? file(params.hmm_dir) : null },
        Channel.value("initial"),
        Channel.value(outdir),
        run_prefilter_hmmsearch
    )

    // Gather the meta info for the gene family search results
    defline_info_collected = PREFILTER_SEARCH.out.defline_info
        .collectFile(name: 'all_defline_info.csv', keepHeader: false, skip: 1)

    // Process the collected information into a metamap
    defline_info_metamap = defline_info_collected
        .splitCsv(header: ['seq', 'parent_seq', 'clean_seq', 'clean_parent_seq', 'id', 'input_fasta_path', 'gene_family_name', 'preprocessed_fasta_path'])
        .map { row -> 
            def groupKey = [
                id: row.id
            ]
            def seqInfo = [
                seq: row.seq,
                parent_seq: row.parent_seq,
                clean_seq: row.clean_seq,
                clean_parent_seq: row.clean_parent_seq,
                input_fasta_path: row.input_fasta_path,
                gene_family_name: row.gene_family_name,
                preprocessed_fasta_path: row.preprocessed_fasta_path
            ]
            [groupKey, seqInfo]
        }
        .groupTuple(by: 0)

    
    
    fasta_info_collected = PREFILTER_SEARCH.out.fasta_info
        .collectFile(name: 'all_fasta_info.csv', keepHeader: false, skip: 1)
    
    fasta_info_metamap = fasta_info_collected
        .splitCsv(header: ['id', 'input_fasta_path', 'gene_family_name', 'preprocessed_fasta_path'])
        .map { row -> 
            def groupKey = [
                    id: row.id
                ]
            def seqInfo = [
                input_fasta_path: row.input_fasta_path,
                gene_family_name: row.gene_family_name,
                preprocessed_fasta_path: row.preprocessed_fasta_path
            ]
            [groupKey, seqInfo]
        }
        .groupTuple(by: 0)
    
    cleanfastas_collected = PREFILTER_SEARCH.out.cleanfasta.collect()
    //cleanfastas_collected.view { it -> "cleanfastas_collected initial: $it" }

    if (run_prefilter_hmmsearch) {
        // If hmmsearch was run, then we merge the outputs fastas for all gene families per genome
        preprocessed_fastas_by_id = fasta_info_metamap
        .map { entry ->
            def id = entry[0].id
            def fastas = entry[1].collect { it.preprocessed_fasta_path }
            return tuple(id, fastas)
        }

        CONCATENATE_FASTAS(preprocessed_fastas_by_id)

        cleanfastas_collected = CONCATENATE_FASTAS.out.concatenated_fasta
            .map { id, path -> [[id: id], path.toString()] }
            .collect()
    }


    //defline_info_metamap.view { it -> "defline_info_metamap: $it" }
    //fasta_info_metamap.view { it -> "fasta_info_metamap: $it" }
    //cleanfastas_collected.view { it -> "cleanfastas_collected final: $it" }

    emit:
    defline_info_metamap = defline_info_metamap // metainfo to later identify gene families
    fasta_info_metamap = fasta_info_metamap // metainfo to regroup fastas and later identify gene families
    cleanfastas_collected = cleanfastas_collected // this will be the set of fastas analysed downstream

}

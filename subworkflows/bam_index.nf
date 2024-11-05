#!/usr/bin/env nextflow


// Include modules
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

workflow BAM_INDEX {
    take:
    sample_sheet

    main:
    // Create input channel from a text file listing input file paths
    bam_ch = Channel.fromPath(sample_sheet)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.bam))
        }
    
    SAMTOOLS_INDEX(bam_ch)

    emit:
    bai      = SAMTOOLS_INDEX.out.bai
    versions = SAMTOOLS_INDEX.out.versions
}
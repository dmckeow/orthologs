process PARSE_MMSEQS_TO_FASTA {

    conda "${projectDir}/modules/local/cluster_mmseqs/environment.yml"

    input:
    tuple val(meta), path(mmseqstsv)
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path("${resultdb}.abc"), emit: abc
    tuple val(meta), path("*.*.mmseqs.fa"), emit: fasta

    script:
    def mmseqstsv_basename = mmseqstsv.getBaseName()
    
    """
    python ${projectDir}/bin/parse_mmseqs.py \
        -i ${mmseqstsv}
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m ${mmseqstsv_basename}.abc \
        -f ${fasta_file} \
        -o ./
    """
}
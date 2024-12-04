process CALCULATE_JACCARD {
    
    conda "${moduleDir}/calc_jaccard.yml"
    
    publishDir "${params.outdir}/report/jaccard_analysis", mode: 'copy'

    input:
    path orthogroups_file

    output:
    path "jaccard_similarity_matrix.csv"
    path "jaccard_similarity_heatmap.png"

    script:
    """
    python3 ${projectDir}/bin/calculate_jaccard.py \
        ${orthogroups_file} \
        jaccard_similarity_matrix.csv \
        jaccard_similarity_heatmap.png
    """
}
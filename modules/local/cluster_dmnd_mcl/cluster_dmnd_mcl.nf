process CLUSTER_DMND_MCL {
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_diamond_mcl_python:babf799e599f4f74'

    input:
    tuple val(meta), path(fasta_files)
    val(diamond_params)
    val(mcl_params)
    val(mcl_inflation)
    val(source)
    
    output:
    tuple val(meta), path("*.dmnd.csv"), emit: dmnd_csv
    tuple val(meta), path("*.dmnd.csv.abc"), emit: mcl_abc
    tuple val(meta), path("*.fa"), optional: true, emit: dmnd_mcl_fastas
    tuple val(meta), path("dmnd_mcl.orthogroups.txt"), optional: true, emit: dmnd_mcl_orthogroups

    script:
    def input_fasta = fasta_files instanceof List ? "combined_input.fasta" : fasta_files
    def output_prefix = fasta_files instanceof List ? "combined_input" : fasta_files.baseName
    """
    if [ -f "${input_fasta}" ]; then
        echo "single fasta input for dmnd mcl clustering"
        #ln -s ${fasta_files} ${input_fasta}
    else
        cat ${fasta_files} > ${input_fasta}
    fi
    
    diamond blastp ${diamond_params} -d ${input_fasta} -q ${input_fasta} -o ${output_prefix}.dmnd.csv --threads ${task.cpus}
    
    awk '{ print \$1,\$2,\$12 }' ${output_prefix}.dmnd.csv > ${output_prefix}.dmnd.csv.tmp
    mv ${output_prefix}.dmnd.csv.tmp ${output_prefix}.dmnd.csv
    
    mcl ${output_prefix}.dmnd.csv ${mcl_params} -I ${mcl_inflation} -o ${output_prefix}.dmnd.csv.abc --abc
    
    python3 ${projectDir}/bin/parse_fastas_mcl.py \
        -m ${output_prefix}.dmnd.csv.abc \
        -f ${input_fasta} \
        -o ./ \
        --source ${source}
    
    mv clusters.abc.tmp dmnd_mcl.orthogroups.txt

    """
}
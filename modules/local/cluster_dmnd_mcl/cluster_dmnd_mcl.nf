//process CLUSTER_DMND_MCL {
//    conda "${projectDir}/modules/local/cluster_dmnd_mcl/environment.yml"
//
//    input:
//    tuple val(meta), path(fasta_file)
//    val(diamond_params)
//    val(mcl_params)
//    val(mcl_inflation)
//    val(input_source)
//    
//    output:
//    tuple val(meta), path("${input_source}/${fasta_file.baseName}.dmnd.csv"), emit: dmnd_csv
//    tuple val(meta), path("${input_source}/${fasta_file.baseName}.dmnd.csv.abc"), emit: mcl_abc
//    tuple val(meta), path("${input_source}/*.mcl.fa"), optional: true, emit: dmnd_mcl_fastas
//
//    script:
//    """
//    mkdir -p ${input_source}
//    
//    diamond blastp ${diamond_params} -d ${fasta_file} -q ${fasta_file} -o ${input_source}/${fasta_file.baseName}.dmnd.csv --threads ${task.cpus}
//    
//    awk '{ print \$1,\$2,\$12 }' ${input_source}/${fasta_file.baseName}.dmnd.csv > ${input_source}/${fasta_file.baseName}.dmnd.csv.tmp
//    mv ${input_source}/${fasta_file.baseName}.dmnd.csv.tmp ${input_source}/${fasta_file.baseName}.dmnd.csv
//    
//    mcl ${input_source}/${fasta_file.baseName}.dmnd.csv ${mcl_params} -I ${mcl_inflation} -o ${input_source}/${fasta_file.baseName}.dmnd.csv.abc --abc
//    
//    python ${projectDir}/bin/parse_fastas_mcl.py \
//        -m ${input_source}/${fasta_file.baseName}.dmnd.csv.abc \
//        -f ${fasta_file} \
//        -o ${input_source}/
//    """
//}

process CLUSTER_DMND_MCL {
    conda "${projectDir}/modules/local/cluster_dmnd_mcl/environment.yml"

    input:
    tuple val(meta), path(fasta_files)
    val(diamond_params)
    val(mcl_params)
    val(mcl_inflation)
    val(input_source)
    
    output:
    tuple val(meta), path("${input_source}/*.dmnd.csv"), emit: dmnd_csv
    tuple val(meta), path("${input_source}/*.dmnd.csv.abc"), emit: mcl_abc
    tuple val(meta), path("${input_source}/*.mcl.fa"), optional: true, emit: dmnd_mcl_fastas

    script:
    def input_fasta = fasta_files instanceof List ? "combined_input.fasta" : fasta_files
    def output_prefix = fasta_files instanceof List ? "combined_input" : fasta_files.baseName
    """
    mkdir -p ${input_source}
    
    if [ -f "${input_fasta}" ]; then
        echo "single fasta input for dmnd mcl clustering"
        #ln -s ${fasta_files} ${input_fasta}
    else
        cat ${fasta_files} > ${input_fasta}
    fi
    
    diamond blastp ${diamond_params} -d ${input_fasta} -q ${input_fasta} -o ${input_source}/${output_prefix}.dmnd.csv --threads ${task.cpus}
    
    awk '{ print \$1,\$2,\$12 }' ${input_source}/${output_prefix}.dmnd.csv > ${input_source}/${output_prefix}.dmnd.csv.tmp
    mv ${input_source}/${output_prefix}.dmnd.csv.tmp ${input_source}/${output_prefix}.dmnd.csv
    
    mcl ${input_source}/${output_prefix}.dmnd.csv ${mcl_params} -I ${mcl_inflation} -o ${input_source}/${output_prefix}.dmnd.csv.abc --abc
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m ${input_source}/${output_prefix}.dmnd.csv.abc \
        -f ${input_fasta} \
        -o ${input_source}/
    """
}
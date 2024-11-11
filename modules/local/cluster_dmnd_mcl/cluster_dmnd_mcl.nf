process CLUSTER_DMND_MCL {
    conda "${projectDir}/modules/local/cluster_dmnd_mcl/environment.yml"

    input:
    tuple val(meta), path(fasta_file)
    val(diamond_params)
    val(mcl_params)
    val(mcl_inflation)
    val(input_source)
    
    output:
    tuple val(meta), path("clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv"), emit: dmnd_csv
    tuple val(meta), path("clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv.abc"), emit: mcl_abc
    tuple val(meta), path("clusters/dmnd_mcl/${input_source}/*.mcl.fa"), optional: true, emit: dmnd_mcl_fastas

    script:
    """
    mkdir -p clusters/dmnd_mcl/${input_source}
    
    diamond blastp ${diamond_params} -d ${fasta_file} -q ${fasta_file} -o clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv --threads ${task.cpus}
    
    awk '{ print \$1,\$2,\$12 }' clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv > clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv.tmp
    mv clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv.tmp clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv
    
    mcl clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv ${mcl_params} -I ${mcl_inflation} -o clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv.abc --abc
    
    python ${projectDir}/bin/parse_fastas_mcl.py \
        -m clusters/dmnd_mcl/${input_source}/${fasta_file.baseName}.dmnd.csv.abc \
        -f ${fasta_file} \
        -o clusters/dmnd_mcl/${input_source}/
    """
}
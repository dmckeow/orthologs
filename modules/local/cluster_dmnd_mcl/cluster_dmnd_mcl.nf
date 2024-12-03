process CLUSTER_DMND_MCL {
    conda "${moduleDir}/environment.yml"
    //container 'community.wave.seqera.io/library/biopython_diamond_mcl_python:babf799e599f4f74'
    //container 'oras://community.wave.seqera.io/library/biopython_diamond_mcl_python:7a0dd2d66fbb321d'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/314ed81fa1c23bc3db25f16df3f5ebf688deee3b47cf00459968b5b81f45c83c/data'

    input:
    tuple val(meta), path(fasta_files)
    val(diamond_params)
    val(mcl_params)
    val(mcl_inflation)
    val(source)
    
    output:
    tuple val(meta), path("diamond/*.dmnd.csv"), emit: dmnd_csv
    tuple val(meta), path("mcl/*.dmnd.csv.abc"), emit: mcl_abc
    tuple val(meta), path("cluster_sequences/*.fa"), optional: true, emit: fasta_files
    tuple val(meta), path("cluster_sequences"), optional: true, emit: fasta_dir
    tuple val(meta), path("dmnd_mcl.orthogroups.txt"), optional: true, emit: orthogroup_list

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
    
    mkdir -p diamond mcl

    diamond blastp ${diamond_params} -d ${input_fasta} -q ${input_fasta} -o diamond/${output_prefix}.dmnd.csv --threads ${task.cpus}
    
    awk '{ print \$1,\$2,\$12 }' diamond/${output_prefix}.dmnd.csv > diamond/${output_prefix}.dmnd.csv.tmp
    mv diamond/${output_prefix}.dmnd.csv.tmp diamond/${output_prefix}.dmnd.csv
    
    mcl diamond/${output_prefix}.dmnd.csv ${mcl_params} -I ${mcl_inflation} -o mcl/${output_prefix}.dmnd.csv.abc --abc
    
    python3 ${projectDir}/bin/parse_fastas_mcl.py \
        -m mcl/${output_prefix}.dmnd.csv.abc \
        -f ${input_fasta} \
        -o cluster_sequences/ \
        --source ${source}
    
    mv clusters.abc.tmp dmnd_mcl.orthogroups.txt

    """
}
process SEARCH {
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    //container 'community.wave.seqera.io/library/bedtools_biopython_hmmer_samtools_python:c2a5749688500b49'
    //container 'oras://community.wave.seqera.io/library/bedtools_biopython_hmmer_samtools_python:6ff18068cb3833ec'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d06f9c1b929618ccc24024d4d13768bd9196a05416f9dcad078853ecbf40efa9/data'

    input:
    tuple val(meta), path(fasta)
    val(gene_family_info)
    val(gene_family_name)
    path(hmm_dir)
    val(input_source)
    val(output_dir)

    output:
    tuple val(meta), path("searches/${input_source}/${meta.id}.txt"), emit: txt
    tuple val(meta), path("searches/${input_source}/*.domains.fasta"), optional: true, emit: domfasta
    tuple val(meta), path("searches/${input_source}/*.domtable"), optional: true, emit: domtable
    tuple val(meta), path("defline_info.csv"), optional: true, emit: defline_info
    
    script:
    
    """
    python3 ${projectDir}/bin/search.py \
        -f ${fasta} \
        -g ${gene_family_info} \
        -t ${task.cpus} \
        ${gene_family_name} \
        -H ${hmm_dir} \
        -i ${input_source} \
        -o searches

    touch searches/${input_source}/${meta.id}.txt

    # Extract deflines and create CSV with full paths
    echo "seq,parent_seq,fa_name,gene_family_name,id,input_fasta_path,domain_fasta_path" > defline_info.csv
    for domain_fasta in searches/${input_source}/*.domains.fasta; do
        grep ">" \$domain_fasta | sed 's/>//g' | awk -v OFS=',' \
        -v id="${meta.id}" \
        -v fa_name="${fasta}" \
        -v gf="${gene_family_name}" \
        -v input_path="\$(readlink -f ${fasta})" \
        -v domain_path="\$(readlink -f \$domain_fasta)" \
        '{
            original=\$0;
            gsub(/:.*/, "", \$0);
            print original, \$0, fa_name, gf, id, input_path, domain_path
        }' >> defline_info.csv
    done
    """
}
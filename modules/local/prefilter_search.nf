process PREFILTER_SEARCH {
    label 'process_low'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d06f9c1b929618ccc24024d4d13768bd9196a05416f9dcad078853ecbf40efa9/data'

    input:
    tuple val(meta), path(fasta)
    val(gene_family_info)
    val(gene_family_name)
    path(hmm_dir)
    val(input_source)
    val(output_dir)
    val(prefilter_hmmsearch)

    output:
    tuple val(meta), path("prefilter/${input_source}/*.domtable"), optional: true, emit: domtable
    tuple val(meta), path("prefilter/${input_source}/defline_info/*"), emit: defline_info
    tuple val(meta), path("prefilter/${input_source}/fasta_info/*"), emit: fasta_info
    tuple val(meta), path("prefilter/${input_source}/clean_fasta/*"), emit: cleanfasta
    
    script:
    def search_cmd = prefilter_hmmsearch ?
        """
        echo "RUNNING PREFILTER WITH SEARCH"
        python3 ${projectDir}/bin/search.py \
            -f ${fasta} \
            -g ${gene_family_info} \
            -t ${task.cpus} \
            ${gene_family_name} \
            -H ${hmm_dir} \
            -i ${input_source} \
            -o prefilter
        # make copy of fasta with cleaned deflines

        for f in prefilter/${input_source}/*.domains.fasta; do
            awk '/^>/ {gsub(/[^a-zA-Z0-9_>]/, "_", \$0)} {print}' \$f \
            > prefilter/${input_source}/clean_fasta/\$(basename \$f)
        done

        # Extract deflines and create CSV with full paths
        echo "seq,parent_seq,clean_seq,clean_parent_seq,id,input_fasta_path,gene_family_name,preprocessed_fasta_path" > \$out_defline_info
            for f in prefilter/${input_source}/*.domains.fasta; do
                grep ">" \$f | sed 's/>//g' | awk -v OFS=',' \
                -v id="${meta.id}" \
                -v gf="${gene_family_name}" \
                -v input_path="\$(readlink -f ${fasta})" \
                -v preprocessed_fasta_path="\$(realpath prefilter/${input_source}/clean_fasta/\$(basename \$f))" \
                '{
                    original=\$0;
                    parent=\$0;
                    gsub(/:.*/, "", parent);
                    clean_parent=parent;
                    clean_original=original;
                    gsub(/[^a-zA-Z0-9_]/, "_", clean_parent);
                    gsub(/[^a-zA-Z0-9_]/, "_", clean_original);
                    print original, parent, clean_original, clean_parent, id, input_path, gf, preprocessed_fasta_path
                }' >>  \$out_defline_info
            done

        sed '1,1d' \$out_defline_info | \
            cut -d "," --complement -f 1-4 | \
            sort -Vu | \
            sed -z 's/^/id,input_fasta_path,gene_family_name,preprocessed_fasta_path\\n/g' \
        > \$out_fasta_info
        
        """ :
        """
        echo "RUNNING PREFILTER WITHOUT SEARCH"
        # Make a copy of the fasta with clean deflines
        awk '/^>/ {gsub(/[^a-zA-Z0-9_>]/, "_", \$0)} {print}' ${fasta} \
            > prefilter/${input_source}/clean_fasta/\$(basename ${fasta})

        # Extract deflines and create CSV with full paths
        echo "seq,parent_seq,clean_seq,clean_parent_seq,id,input_fasta_path,gene_family_name,preprocessed_fasta_path" > \$out_defline_info
            for f in ${fasta}; do
                grep ">" \$f | sed 's/>//g' | awk -v OFS=',' \
                -v id="${meta.id}" \
                -v gf="${gene_family_name}" \
                -v input_path="\$(readlink -f ${fasta})" \
                -v preprocessed_fasta_path="\$(realpath prefilter/${input_source}/clean_fasta/\$(basename \$f))" \
                '{
                    original=\$0;
                    parent=\$0;
                    gsub(/:.*/, "", parent);
                    clean_parent=parent;
                    clean_original=original;
                    gsub(/[^a-zA-Z0-9_>]/, "_", clean_parent);
                    gsub(/[^a-zA-Z0-9_>]/, "_", clean_original);
                    print original, parent, clean_original, clean_parent, id, input_path, gf, preprocessed_fasta_path
                }' >> \$out_defline_info
            done

        sed '1,1d' \$out_defline_info | \
            cut -d "," --complement -f 1-4 | \
            sort -Vu | \
            sed -z 's/^/id,input_fasta_path,gene_family_name,preprocessed_fasta_path\\n/g' \
        > \$out_fasta_info

        """

    """
    #!/bin/bash
    out_defline_info="prefilter/${input_source}/defline_info/${meta.id}.${gene_family_name}.csv"
    out_fasta_info="prefilter/${input_source}/fasta_info/${meta.id}.${gene_family_name}.csv"
    mkdir -p prefilter/${input_source}
    mkdir -p prefilter/${input_source}/defline_info
    mkdir -p prefilter/${input_source}/fasta_info
    mkdir -p prefilter/${input_source}/clean_fasta

    echo ${search_cmd}
    ${search_cmd}
    """
}
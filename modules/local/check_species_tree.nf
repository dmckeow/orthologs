process CHECK_SPECIES_TREE {
    tag "Validating species tree"
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3dd76a137e4bfd2df8b00d8e07e5ffc355dc680adccf94bfbfbc2b5bbdef9efe/data'

    input:
    path samplesheet
    path species_tree

    output:
    path "species_tree_check.log"
    
    script:
    """
    # Run the Python script and capture its output
    python ${projectDir}/bin/check_species_tree.py ${samplesheet} ${species_tree}
    """
}
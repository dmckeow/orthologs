process FILTER_ORTHOGROUPS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::biopython"

    input:
    tuple val(meta), path(orthogroup_dir)
    val min_sequences

    output:
    tuple val(meta), path("filtered_orthogroups"), emit: filtered_orthogroups
    tuple val(meta), path("removed_orthogroups"), emit: removed_orthogroups
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import os
    import shutil

    os.makedirs("filtered_orthogroups", exist_ok=True)
    os.makedirs("removed_orthogroups", exist_ok=True)

    for filename in os.listdir("${orthogroup_dir}"):
        if filename.endswith(('.fa', '.faa', '.fasta', '.fas')):
            file_path = os.path.join("${orthogroup_dir}", filename)
            sequences = list(SeqIO.parse(file_path, "fasta"))
            if len(sequences) > ${min_sequences}:
                shutil.copy2(file_path, "filtered_orthogroups")
                print("Kept: ", filename)
            else:
                shutil.move(file_path, "removed_orthogroups")
                print("Removed: ", filename)

    with open("versions.yml", "w") as f:
        f.write("'${task.process}':\\n  python: \$( python --version | sed 's/Python //g' )\\n  biopython: \$( python -c 'import Bio; print(Bio.__version__)' )\\n")
    """
}
process CALCULATE_JACCARD {

    conda "${moduleDir}/calc_jaccard.yml"

    publishDir "${params.outdir}/jaccard_analysis", mode: 'copy'

    input:
    path orthogroups_tool1
    path orthogroups_tool2
    val tool1_name
    val tool2_name

    output:
    path "${tool1_name}_vs_${tool2_name}_jaccard.csv"
    path "${tool1_name}_vs_${tool2_name}_jaccard_heatmap.png"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.spatial.distance import jaccard

    def read_orthogroups(file_path):
        orthogroups = {}
        with open(file_path, 'r') as f:
            for line in f:
                og, proteins = line.strip().split(': ')
                orthogroups[og] = set(proteins.split())
        return orthogroups

    og1 = read_orthogroups('${orthogroups_tool1}')
    og2 = read_orthogroups('${orthogroups_tool2}')

    jaccard_matrix = np.zeros((len(og1), len(og2)))

    for i, (og1_name, og1_proteins) in enumerate(og1.items()):
        for j, (og2_name, og2_proteins) in enumerate(og2.items()):
            jaccard_matrix[i, j] = 1 - jaccard(og1_proteins, og2_proteins)

    df = pd.DataFrame(jaccard_matrix, index=og1.keys(), columns=og2.keys())
    df.to_csv('${tool1_name}_vs_${tool2_name}_jaccard.csv')

    plt.figure(figsize=(20, 20))
    sns.heatmap(df, cmap='YlOrRd', xticklabels=False, yticklabels=False)
    plt.title('Jaccard Similarity Heatmap: ${tool1_name} vs ${tool2_name}')
    plt.tight_layout()
    plt.savefig('${tool1_name}_vs_${tool2_name}_jaccard_heatmap.png', dpi=300)
    """
}
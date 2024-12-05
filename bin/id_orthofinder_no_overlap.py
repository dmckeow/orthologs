import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze genome overlap matrix to identify pairs with zero intersection and determine which genomes to keep.")
    parser.add_argument('-i', '--input', required=True, help="Input file containing the genome overlap matrix (tab-separated)")
    parser.add_argument('-o', '--output', default='no_overlaps.txt', help="Output file to write results (optional, default is stdout)")
    return parser.parse_args()

def count_total_overlaps(df, genome):
    return df.loc[genome].sum() - df.loc[genome, genome]

def analyze_matrix(df):
    zero_pairs = []
    for i in range(len(df.index)):
        for j in range(i+1, len(df.index)):
            if df.iloc[i, j] == 0:
                zero_pairs.append((df.index[i], df.index[j]))

    genomes_to_eliminate = set()
    for genome1, genome2 in zero_pairs:
        overlaps1 = count_total_overlaps(df, genome1)
        overlaps2 = count_total_overlaps(df, genome2)
        if overlaps1 >= overlaps2:
            genomes_to_eliminate.add(genome2)
        else:
            genomes_to_eliminate.add(genome1)

    genomes_to_keep = set(df.index) - genomes_to_eliminate
    return zero_pairs, genomes_to_eliminate, genomes_to_keep

def main():
    args = parse_arguments()

    # Read the matrix from the input file
    df = pd.read_csv(args.input, sep='\t', index_col=0)

    zero_pairs, genomes_to_eliminate, genomes_to_keep = analyze_matrix(df)

    # Prepare the output
    output = []
    output.append("Pairs with zero intersection:")
    for pair in zero_pairs:
        output.append(f"{pair[0]} - {pair[1]}")

    output.append("\nGenomes to eliminate:")
    for genome in genomes_to_eliminate:
        output.append(genome)

    output.append("\nGenomes to keep:")
    for genome in genomes_to_keep:
        output.append(genome)

    
    print('\n'.join(output)) # For log file

    # Save the genomes to eliminate for other process
    with open(args.output, 'w') as f:
        f.write('\n'.join(genomes_to_eliminate) + '\n')

if __name__ == "__main__":
    main()
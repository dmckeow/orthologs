#!/usr/bin/env python3
import sys
import random
from Bio import SeqIO
# not used in workflow - just used to reduce input data for testing purposes

def extract_random_sequences(input_file, output_file, percentage):
    with open(input_file, 'r') as f:
        sequences = list(SeqIO.parse(f, 'fasta'))
    
    num_sequences = len(sequences)
    num_extract = int(num_sequences * percentage / 100)
    print(f"Num seqs to extract ({percentage}% of {num_sequences}): {num_extract}")
    
    selected_sequences = random.sample(sequences, num_extract)
    
    with open(output_file, 'w') as f:
        SeqIO.write(selected_sequences, f, 'fasta')

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    percentage = float(sys.argv[3])
    extract_random_sequences(input_file, output_file, percentage)
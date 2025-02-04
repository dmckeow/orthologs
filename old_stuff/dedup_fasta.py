#!/usr/bin/env python3
import sys
from Bio import SeqIO

def deduplicate_fasta(input_file, output_file):
    # Dictionary to store unique sequences
    unique_sequences = {}

    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        # Convert sequence to string for comparison
        sequence = str(record.seq)
        
        # If the sequence is not in the dictionary, add it
        if sequence not in unique_sequences:
            unique_sequences[sequence] = record

    # Write unique sequences to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(unique_sequences.values(), output_handle, "fasta")

    print(f"Deduplicated sequences written to {output_file}")

if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_file.fasta output_file.fasta")
        sys.exit(1)

    # Get input and output file names from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Run the deduplication function
    deduplicate_fasta(input_file, output_file)
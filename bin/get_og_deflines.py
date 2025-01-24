import os
import re
import argparse
from collections import defaultdict

def process_fasta_files(input_directory, output_directory):
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Set default output file names
    output_file = os.path.join(output_directory, "Orthogroups.tsv")
    gene_count_file = os.path.join(output_directory, "Orthogroups.GeneCount.tsv")

    # Initialize a dictionary to hold the deflines for each sample
    data = defaultdict(lambda: defaultdict(list))

    # Regular expression to extract sample name from the fasta defline
    sample_regex = re.compile(r">([^_]+)")

    # Loop through all files in the input directory
    for filename in os.listdir(input_directory):
        # Only process .fa files
        if filename.endswith(".fa"):
            orthogroup = filename.split('.')[0]  # Extract orthogroup from filename
            
            # Open and process the fasta file
            with open(os.path.join(input_directory, filename), 'r') as file:
                current_sample = None
                # Iterate over each line in the file
                for line in file:
                    line = line.strip()
                    if line.startswith(">"):
                        # Extract the sample name from the defline
                        match = sample_regex.match(line)
                        if match:
                            current_sample = match.group(1)  # The sample name

                            # Add the defline to the appropriate sample in the dictionary
                            data[orthogroup][current_sample].append(line[1:]) # leading > is removed here
                    elif current_sample:
                        # skip if it is sequence
                        continue

    # Prepare the output
    header = ["Orthogroup"]  # Start with the Orthogroup column
    samples = sorted(set(sample for orthogroup in data.values() for sample in orthogroup.keys()))
    header.extend(samples)  # Add sample names as headers
    header.append("Total")  # Add Total column to header

    # Initialize the rows
    rows = []
    gene_count_rows = []

    # Loop through the orthogroups and create a row for each
    for orthogroup, samples_data in sorted(data.items()): # sorted here by orthogroup
        row = [orthogroup]
        gene_count_row = [orthogroup]
        total_count = 0  # Initialize total count for each orthogroup
        # For each sample, join its fasta deflines with a comma and space
        for sample in samples:
            if sample in samples_data:
                deflines = ", ".join(samples_data[sample])
                row.append(deflines)
                count = len(samples_data[sample])
                gene_count_row.append(str(count))  # Count of deflines
                total_count += count  # Add to total count
            else:
                row.append("")  # Empty if there are no deflines for that sample
                gene_count_row.append("0")  # Zero count if there are no deflines
        gene_count_row.append(str(total_count))  # Add total count to gene count row
        rows.append(row)
        gene_count_rows.append(gene_count_row)

    # Write the output to a tab-separated file
    with open(output_file, "w") as out_file:
        out_file.write("\t".join(header[:-1]) + "\n")  # Exclude 'Total' from this file
        for row in rows:
            out_file.write("\t".join(row) + "\n")

    # Write the gene count output to a tab-separated file
    with open(gene_count_file, "w") as out_file:
        out_file.write("\t".join(header) + "\n")  # Include 'Total' in this file
        for row in gene_count_rows:
            out_file.write("\t".join(row) + "\n")

def main():
    # Set up argparse for command-line arguments
    parser = argparse.ArgumentParser(description="Process fasta files to generate orthogroup summary.")
    parser.add_argument("input_directory", type=str, help="Path to the directory containing fasta files.")
    parser.add_argument("-o", "--output_directory", type=str, default=".", help="Path to the output directory. Default is the current directory.")
    
    args = parser.parse_args()
    
    # Call the function to process the fasta files in the provided directory
    process_fasta_files(args.input_directory, args.output_directory)

if __name__ == "__main__":
    main()
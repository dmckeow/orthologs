import os
import argparse
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial
import itertools

# OLD:
def mcl_input(input_file, output_file):
	# Extract the base file name without extension
	#base_filename = input_file.split('/')[-1].replace('.dmnd.csv.abc', '')
	base_filename = os.path.splitext(os.path.basename(input_file))[0]
	
	# Initialize an empty list to store the modified lines
	modified_lines = []
	
	# Add the header line
	modified_lines.append("#placeholder\theader")  # Modify headers as needed
	
	# Read the input file and modify the lines
	with open(input_file, 'r') as file:
		for index, line in enumerate(file):
			# Replace all tabs with spaces
			line = line.replace('\t', ' ')
			# Create the new identifier
			new_id = f"{base_filename}.{index + 1}.mcl"
			# Combine the new identifier with the modified line
			modified_line = f"{new_id}\t{line.strip()}"
			modified_lines.append(modified_line)
	
	# Write the modified lines to the output file
	with open(output_file, 'w') as file:
		for modified_line in modified_lines:
			file.write(modified_line + '\n')
def mmseqs_input(input_file, output_file):
	# Extract the base file name without extension
	#base_filename = input_file.split('/')[-1].replace('.cluster.abc', '')
	base_filename = os.path.splitext(os.path.basename(input_file))[0]
	
	# Initialize an empty list to store the modified lines
	modified_lines = []
	
	# Add the header line
	modified_lines.append("#placeholder\theader")  # Modify headers as needed
	
	# Read the input file and modify the lines
	with open(input_file, 'r') as file:
		for index, line in enumerate(file):
			# Replace all tabs with spaces
			line = line.replace('\t', ' ')
			# Create the new identifier
			new_id = f"{base_filename}.{index + 1}.mmseqs"
			# Combine the new identifier with the modified line
			modified_line = f"{new_id}\t{line.strip()}"
			modified_lines.append(modified_line)
	
	# Write the modified lines to the output file
	with open(output_file, 'w') as file:
		for modified_line in modified_lines:
			file.write(modified_line + '\n')
def extract_fasta_sequences(og_file, fasta_file, output_dir):
	# Read the OG to protein names mapping from the input file
	with open(og_file, "r") as og:
		for line in og.readlines()[1:]:  # Skip the header
			og_name, proteins = line.strip().split("\t")
			protein_list = proteins.split(" ")  # Get list of protein names

			# Create a new fasta file for the orthogroup
			output_file = os.path.join(output_dir, f"{og_name}.fa")
			
			with open(output_file, "w") as out_fasta:
				# Parse the single fasta file for matching records
				for record in SeqIO.parse(fasta_file, "fasta"):
					if record.id in protein_list:
						SeqIO.write(record, out_fasta, "fasta")

	print("parse_fastas_mcl.py: Fasta extraction complete.")

# NEW:

def process_input_file(input_file, output_file, source):
    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("#placeholder\theader\n")  # Write header
        
        for index, line in enumerate(infile, 1):
            new_id = f"{base_filename}.{index}.{source}"
            modified_line = f"{new_id}\t{line.replace('\t', ' ').strip()}\n"
            outfile.write(modified_line)

def extract_fasta_sequences(args):
    og_name, proteins, fasta_file, output_dir = args
    output_file = os.path.join(output_dir, f"{og_name}.cluster.fa")
    
    protein_set = set(proteins.split())
    
    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in protein_set:
                SeqIO.write(record, out_fasta, "fasta")

def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    output = os.path.splitext(os.path.basename(args.mcl_ogs))[0] + ".tmp"

    # Process input file
    process_input_file(args.mcl_ogs, output, args.source)

    # Read OG to protein names mapping
    with open(output, "r") as og_file:
        next(og_file)  # Skip header
        og_data = [line.strip().split("\t") for line in og_file]

    # Prepare arguments for multiprocessing
    extract_args = [(og_name, proteins, args.fasta, args.outdir) for og_name, proteins in og_data]

    # Use multiprocessing to extract sequences
    with Pool() as pool:
        pool.map(extract_fasta_sequences, extract_args)

    print(f"parse_fastas_{args.source}.py: Fasta extraction complete.")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-f','--fasta', required=True, help='Path to the input fasta file containing sequences')
	parser.add_argument('-m', '--mcl_ogs', required=True, help='Path to the .abc cluster file from MCL')
	parser.add_argument('-o', '--outdir', required=True, help='Path to the directory where mcl orthogroup fastas will be put')
	parser.add_argument(
		"--source",
		choices=["mcl", "mmseqs"],
		required=True,
		help="Choose either 'mcl' or 'mmseqs'."
	)
	
	args = parser.parse_args()
	main(args)

	# OLD:
	#os.makedirs(args.outdir, exist_ok=True)
	#output = os.path.splitext(os.path.basename(args.mcl_ogs))[0] + ".tmp"
	#if args.source == "mcl":
	#	mcl_input(args.mcl_ogs, output)
	#elif args.source == "mmseqs":
	#	mmseqs_input(args.mcl_ogs, output)
	#extract_fasta_sequences(output, args.fasta, args.outdir)

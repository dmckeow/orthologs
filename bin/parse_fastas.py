import os
import argparse
from Bio import SeqIO


def extract_fasta_sequences(og_file, fasta_dir, output_dir):
	# Read the OG to protein names mapping from the input file
	with open(og_file, "r") as og:
		for line in og.readlines()[1:]:  # Skip the header
			og_name, proteins = line.strip().split("\t")
			protein_list = proteins.split(" ")  # Get list of protein names

			# Create a new fasta file for the orthogroup
			output_file = os.path.join(output_dir, f"{og_name}.fa")
			
			with open(output_file, "w") as out_fasta:
				# Loop through each fasta file in the fasta_dir
				for fasta_file in os.listdir(fasta_dir):
					if fasta_file.endswith(".fa") or fasta_file.endswith(".fasta"):
						fasta_path = os.path.join(fasta_dir, fasta_file)
						# search for matching fastas between orthogroup list and input fastas in dir
						for record in SeqIO.parse(fasta_path, "fasta"):
							if record.id in protein_list:
								SeqIO.write(record, out_fasta, "fasta")

	print("parse_fastas.py: Fasta extraction complete.")

def process_of_hogs(input_file, output_file):
	with open(input_file, 'r') as infile, open (output_file, 'w') as outfile:
		header = infile.readline() # skip header by reading the first line to header
		for line in infile:
			columns = line.split('\t')
			columns = [column for i, column in enumerate(columns) if i not in [1, 2]] # remove columns 2 and 3
			first_column = columns[0]
			other_columns = ' '.join(columns[1:])
			result = first_column + '\t' + other_columns
			result = result.replace('N0.', '', 1)
			outfile.write(result + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-f','--fasta', required=True, help='Path to the input fasta file directory containing .fa or .fasta files')
	parser.add_argument('-i', '--input', required=True, help='Path to dir_step3/orthologous_groups.txt (Broccoli) OR orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv (OrthoFinder)')
	parser.add_argument('-o', '--outdir', required=True, help='Path to the directory where fastas will be put')
	
	args = parser.parse_args()

	os.makedirs(args.outdir, exist_ok=True)

	if os.path.basename(args.input) == "N0.tsv":
		processed_n0 = os.path.join(args.outdir, "N0.txt")
		process_of_hogs(args.input, processed_n0)
		extract_fasta_sequences(processed_n0, args.fasta, args.outdir)
	else:
		extract_fasta_sequences(args.input, args.fasta, args.outdir)

import os
import argparse
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial
import itertools

#def process_input_file(input_file, output_file, source):
def process_input_file(input_file, output_file):
	base_filename = os.path.splitext(os.path.basename(input_file))[0]
	
	with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
		outfile.write("#placeholder\theader\n")  # Write header
		
		for index, line in enumerate(infile, 1):
			#new_id = f"{base_filename}.{index}.{source}"
			#new_id = f"CL_{index}" # this is only valid in python 3.6 and over
			new_id = "CL_" + str(index)
			#modified_line = f"{new_id}\t{line.replace('\t', ' ').strip()}\n"
			modified_line = new_id + "\t" + line.replace('\t', ' ').strip() + "\n"
			outfile.write(modified_line)

def extract_fasta_sequences(args):
	og_name, proteins, fasta_file, output_dir = args
	output_file = os.path.join(output_dir, og_name + ".fa")
	
	protein_set = set(proteins.split())
	
	with open(output_file, "w") as out_fasta:
		for record in SeqIO.parse(fasta_file, "fasta"):
			if record.id in protein_set:
				SeqIO.write(record, out_fasta, "fasta")

def main(args):
	os.makedirs(args.outdir, exist_ok=True)
	#output = os.path.splitext(os.path.basename(args.mcl_ogs))[0] + ".abc.tmp"
	output = os.path.join(args.outdir, "clusters.abc.tmp")

	# Process input file
	#process_input_file(args.mcl_ogs, output, args.source)
	process_input_file(args.mcl_ogs, output)

	# Read OG to protein names mapping
	with open(output, "r") as og_file:
		next(og_file)  # Skip header
		og_data = [line.strip().split("\t") for line in og_file]

	# Prepare arguments for multiprocessing
	extract_args = [(og_name, proteins, args.fasta, args.outdir) for og_name, proteins in og_data]

	# Use multiprocessing to extract sequences
	with Pool() as pool:
		pool.map(extract_fasta_sequences, extract_args)

	print("parse_fastas_" + args.source + ".py: Fasta extraction complete.")

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

	

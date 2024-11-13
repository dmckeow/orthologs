import os
import argparse
import subprocess
import glob

def do_mmseqs_createdb(fasta, outdb):
	cmd = (
		f"mmseqs createdb {fasta} {outdb}"
	)
	subprocess.run(cmd, shell=True, check=True)
	
def do_mmseqs_createdb(fasta_files, outdb):
	# Concatenate all fasta files into a single temporary file
	with open('temp_combined.fasta', 'w') as outfile:
		for fasta_file in fasta_files:
			with open(fasta_file, 'r') as infile:
				outfile.write(infile.read())
	
	# Create the MMseqs2 database from the combined file
	cmd = f"mmseqs createdb temp_combined.fasta {outdb}"
	subprocess.run(cmd, shell=True, check=True)
	
	# Remove the temporary file
	os.remove('temp_combined.fasta')

def do_mmseqs_cluster(outdb):
	cmd = f"mmseqs cluster {outdb} {outdb}.cluster {outdb}.tmp"
	subprocess.run(cmd, shell=True, check=True)
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-i','--input_dir', required=True, help='Path to the input fasta file directory containing .fa or .fasta files')
	parser.add_argument('-o', '--outdb', required=True, help='Name for the output database')
	parser.add_argument('-p', '--pattern', default='*.fasta', help='File pattern to match (default: *.fasta)')
	parser.add_argument('-s', '--skip_createdb', action='store_true', help='Skip the createdb step if specified')
	
	args = parser.parse_args()
	
	# Fetch files matching the pattern in the input directory
	fasta_files = glob.glob(os.path.join(args.input_dir, args.pattern))

	if not fasta_files:
		print(f"No files matching the pattern '{args.pattern}' found in '{args.input_dir}'")
		exit(1)
		
	if not args.skip_createdb:
		print("Running createdb")
		do_mmseqs_createdb(fasta_files, args.outdb)
	
	print("Skipping createdb")
	do_mmseqs_cluster(args.outdb)
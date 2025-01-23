import os
import re
import argparse
from collections import defaultdict

def process_fasta_files(directory, output_file):
	# Initialize a dictionary to hold the deflines for each sample
	data = defaultdict(lambda: defaultdict(list))

	# Regular expression to extract sample name from the fasta defline
	sample_regex = re.compile(r">([^_]+)")

	# Loop through all files in the directory
	for filename in os.listdir(directory):
		# Only process .fa files
		if filename.endswith(".fa"):
			orthogroup = filename.split('.')[0]  # Extract orthogroup from filename
			
			# Open and process the fasta file
			with open(os.path.join(directory, filename), 'r') as file:
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

	# Initialize the rows
	rows = []

	# Loop through the orthogroups and create a row for each
	for orthogroup, samples_data in sorted(data.items()): # sorted here by orthogroup
		row = [orthogroup]
		# For each sample, join its fasta deflines with a comma and space
		for sample in samples:
			if sample in samples_data:
				deflines = ", ".join(samples_data[sample])
				row.append(deflines)
			else:
				row.append("")  # Empty if there are no deflines for that sample
		rows.append(row)

	# Write the output to a tab-separated file
	with open(output_file, "w") as out_file:
		out_file.write("\t".join(header) + "\n")
		for row in rows:
			out_file.write("\t".join(row) + "\n")

def main():
	# Set up argparse for command-line arguments
	parser = argparse.ArgumentParser(description="Process fasta files to generate orthogroup summary.")
	parser.add_argument("directory", type=str, help="Path to the directory containing fasta files.")
	parser.add_argument("-o", "--output", type=str, default="Orthogroups.tsv", help="Path to the output file. Default is 'Orthogroups.tsv' in the current directory.")
	
	args = parser.parse_args()
	
	# Call the function to process the fasta files in the provided directory
	process_fasta_files(args.directory, args.output)

if __name__ == "__main__":
	main()

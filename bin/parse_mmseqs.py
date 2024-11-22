import os
import argparse
import subprocess
import glob
import shutil
from collections import defaultdict

def flatten_clusters(input_file, output_file):
	# Dictionary to store clusters
	clusters = defaultdict(set)

	# Read the input file and populate the clusters dictionary
	with open(input_file, 'r') as f:
		for line in f:
			representative, member = line.strip().split('\t')
			clusters[representative].add(member)

	# Write the flattened clusters to the output file
	with open(output_file, 'w') as f:
		for representative, members in clusters.items():
			# Ensure the representative is first in the list
			cluster_list = [representative] + [m for m in members if m != representative]
			f.write('\t'.join(cluster_list) + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	# Search 
	parser.add_argument('-i','--input', required=True, help='Path to the tsv from mmseqs createtsv')

	args = parser.parse_args()

	output = os.path.basename(f"{args.input}.abc")

	flatten_clusters({args.input}, output)
	


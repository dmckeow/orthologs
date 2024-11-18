import os
import argparse
import subprocess
import glob
import shutil
from collections import defaultdict
	
def do_mmseqs_createtsv(querydb, targetdb, resultsdb):
	cmd = f"mmseqs createtsv {querydb} {targetdb} {resultsdb} {resultsdb}.tsv"
	subprocess.run(cmd, shell=True, check=True)

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
	parser.add_argument('-q','--querydb', required=True, help='Path to the query database used in mmseqs')

	parser.add_argument('-t','--targetdb', required=True, help='Path to the target database used in mmseqs - for cluster this is the same as querydb')

	parser.add_argument('-r','--resultdb', required=True, help='Path to the output database from mmseqs using the query and targets dbs provided')

	parser.add_argument('-o', '--outdir', required=True, help='Path to the directory where mcl orthogroup fastas will be put')
	
	args = parser.parse_args()

	os.makedirs(args.outdir, exist_ok=True)
	
	do_mmseqs_createtsv(args.querydb, args.targetdb, args.resultdb)

	flatten_clusters(f"{args.resultdb}.tsv", f"{args.resultdb}.abc")

	# Define paths for the result files
	tsv_file = f"{args.resultdb}.tsv"
	abc_file = f"{args.resultdb}.abc"

	# Move the files to the output directory
	shutil.move(tsv_file, os.path.join(args.outdir, os.path.basename(tsv_file)))
	shutil.move(abc_file, os.path.join(args.outdir, os.path.basename(abc_file)))
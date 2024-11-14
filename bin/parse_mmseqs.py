import os
import argparse
import subprocess
import glob
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
	
	args = parser.parse_args()
	
	do_mmseqs_createtsv(args.querydb, args.targetdb, args.resultdb)

	flatten_clusters(f"{args.resultdb}.tsv", f"{args.resultdb}.abc")
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse

def jaccard_similarity(set1, set2):
	intersection = len(set1.intersection(set2))
	union = len(set1.union(set2))
	return intersection / union if union > 0 else 0

def calculate_jaccard(orthogroups_file, output_matrix, output_heatmap):
	# Read the orthogroups file
	df = pd.read_csv(orthogroups_file)

	# Create a dictionary to store orthogroups for each source
	source_orthogroups = defaultdict(list)

	# Populate the dictionary
	for _, row in df.iterrows():
		source_orthogroups[row['og_source']].append(set(row['deflines'].split()))

	# Get unique sources
	sources = list(source_orthogroups.keys())

	# Initialize Jaccard similarity matrix
	jaccard_matrix = np.zeros((len(sources), len(sources)))

	# Calculate Jaccard similarity
        # For each pair of ortholog sources (source1 and source2), it iterates across the orthogroups in source1 (og1)
        # THEN for each og1 it calcualte the Jaccard similarity of it versus every orthogroup in source2 (og2)
        # it only keeps the highest similarity identified for each og1, adds up a total, and n comparisons (to calc average later)
        # finally, it averages the similarity of orthogroups between the two orthogroup sources
	for i, source1 in enumerate(sources):
		for j, source2 in enumerate(sources):
			if i <= j:  # We only need to calculate the upper triangle
				total_similarity = 0
				comparisons = 0
				for og1 in source_orthogroups[source1]:
					best_similarity = 0
					for og2 in source_orthogroups[source2]:
						similarity = jaccard_similarity(og1, og2)
						best_similarity = max(best_similarity, similarity)
					total_similarity += best_similarity
					comparisons += 1
				avg_similarity = total_similarity / comparisons if comparisons > 0 else 0
				jaccard_matrix[i, j] = avg_similarity
				jaccard_matrix[j, i] = avg_similarity  # Matrix is symmetric

	# Create a DataFrame from the Jaccard matrix
	jaccard_df = pd.DataFrame(jaccard_matrix, index=sources, columns=sources)

	# Save the Jaccard similarity matrix to a CSV file
	jaccard_df.to_csv(output_matrix)

	# Create a heatmap
	plt.figure(figsize=(10, 8))
	sns.heatmap(jaccard_df, annot=True, cmap="YlGnBu", vmin=0, vmax=1, fmt='.2f')
	plt.title("Jaccard Similarity between Orthogroup Sources")
	plt.tight_layout()
	plt.savefig(output_heatmap, dpi=300)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate Jaccard similarity between orthogroup sources')
	parser.add_argument('orthogroups_file', help='Path to the orthogroups CSV file')
	parser.add_argument('output_matrix', help='Path to save the Jaccard similarity matrix CSV')
	parser.add_argument('output_heatmap', help='Path to save the Jaccard similarity heatmap PNG')
	args = parser.parse_args()

	calculate_jaccard(args.orthogroups_file, args.output_matrix, args.output_heatmap)
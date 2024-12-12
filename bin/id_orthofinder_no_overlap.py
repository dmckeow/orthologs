import argparse
import re
from collections import defaultdict
import os
import glob

def parse_arguments():
	parser = argparse.ArgumentParser(description="Analyze OrthoFinder log to identify low overlap pairs and determine which species to keep.")
	parser.add_argument('-l', '--orthofinder_output', required=True, help="Path to the orthofinder output directory")
	parser.add_argument('-o', '--output', default='no_overlaps.txt', help="Output file to write results")
	return parser.parse_args()

def find_command_log(orthofinder_output):
	search_path = os.path.join(orthofinder_output, "*", ".command.log")
	file = glob.glob(search_path)
	if not file:
		raise FileNotFoundError(f"Cannot find {search_path}")
	return file[0]

def find_species_ids(orthofinder_output):
	search_path = os.path.join(orthofinder_output, "*", "orthofinder/WorkingDirectory/SpeciesIDs.txt")
	file = glob.glob(search_path)
	if not file:
		raise FileNotFoundError(f"Cannot find {search_path}")
	return file[0]

def find_species_overlaps(orthofinder_output):
	search_path = os.path.join(orthofinder_output, "*", "orthofinder/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv")
	file = glob.glob(search_path)
	if not file:
		raise FileNotFoundError(f"Cannot find {search_path}")
	return file[0]

def parse_log_file(log_file):
	low_overlap_pairs = []
	pattern = r"WARNING: Too few hits between species (\d+) and species (\d+)"
	
	with open(log_file, 'r') as f:
		for line in f:
			match = re.search(pattern, line)
			if match:
				species1, species2 = int(match.group(1)), int(match.group(2))
				low_overlap_pairs.append((species1, species2))
	
	return low_overlap_pairs

def count_low_overlaps(pairs):
	counts = defaultdict(int)
	for s1, s2 in pairs:
		counts[s1] += 1
		counts[s2] += 1
	return counts

def determine_species_to_remove(pairs, species_overlaps):
	species_in_pairs = set(s for pair in pairs for s in pair)
	to_remove = set()
	
	while pairs:
		# Find the species with the least total overlaps
		species_to_remove = min(species_in_pairs, key=lambda s: species_overlaps[s])
		to_remove.add(species_to_remove)
		
		# Remove all pairs containing the removed species
		pairs = [pair for pair in pairs if species_to_remove not in pair]
		species_in_pairs = set(s for pair in pairs for s in pair)
		
		if not pairs:
			break
	
	return to_remove

def parse_species_ids(species_ids_file):
	species_map = {}
	reverse_species_map = {}
	with open(species_ids_file, 'r') as f:
		for line in f:
			species_id, filename = line.strip().split(': ')
			species_id = int(species_id)
			species_map[species_id] = filename
			reverse_species_map[filename] = species_id
	return species_map, reverse_species_map

def parse_species_overlaps(overlap_matrix_file, species_map):
	species_overlaps = {}
	with open(overlap_matrix_file, 'r') as f:
		header = f.readline().strip().split('\t')
		for line in f:
			parts = line.strip().split('\t')
			species_name = parts[0]
			overlaps = sum(int(x) for x in parts[1:])
			# Find the species ID that corresponds to this species name
			for species_id, filename in species_map.items():
				if species_name in filename:
					species_overlaps[species_id] = overlaps
					break
	return species_overlaps

def main():
	args = parse_arguments()

	log_file = find_command_log(args.orthofinder_output)
	species_ids = find_species_ids(args.orthofinder_output)
	species_overlap_matrix = find_species_overlaps(args.orthofinder_output)
	print("log_file:", log_file)
	print("species_ids:", species_ids)
	print("species_overlap_matrix:", species_overlap_matrix)

	if not os.path.exists(log_file):
		raise FileNotFoundError(f"Cannot find .command.log file at {log_file}")

	low_overlap_pairs = parse_log_file(log_file)
	species_map, reverse_species_map = parse_species_ids(species_ids)
	species_overlaps = parse_species_overlaps(species_overlap_matrix, species_map)

	print("\nTotal species overlap for each species (sorted by overlap):")
	sorted_species = sorted(species_overlaps.items(), key=lambda x: x[1], reverse=True)
	for species_id, overlap in sorted_species:
		species_name = species_map[species_id]
		print(f"Species {species_id} ({species_name}): {overlap}")

	species_to_remove = determine_species_to_remove(low_overlap_pairs, species_overlaps)
	
	files_to_remove = [species_map[species_id] for species_id in species_to_remove if species_id in species_map]

	with open(args.output, 'w') as f:
		for filename in files_to_remove:
			f.write(f"{filename}\n")

	print(f"\nAnalysis complete. {len(files_to_remove)} files identified for removal.")
	print("Removed species:")
	for species_id in species_to_remove:
		if species_id in species_map:
			print(f"  Species {species_id} ({species_map[species_id]})")
	print("Species in low overlap pairs, but spared removal:")
	printed_species = set()
	for pair in low_overlap_pairs:
		for species_id in pair:
			if species_id in species_map and species_id not in species_to_remove and species_id not in printed_species:
				print(f"  Species {species_id} ({species_map[species_id]})")
				printed_species.add(species_id)
	print(f"Remaining low-overlap pairs: {len(low_overlap_pairs) - len(species_to_remove)}")

if __name__ == "__main__":
	main()
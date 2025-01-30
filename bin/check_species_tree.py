#!/usr/bin/env python
from Bio import Phylo
import csv
import argparse

def check_species_tree(samplesheet, tree_file):
    # Read the species tree
    tree = Phylo.read(tree_file, "newick")
    
    # Get all leaf node names from the tree
    tree_species = set(leaf.name for leaf in tree.get_terminals())
    
    # Read species IDs from the samplesheet
    with open(samplesheet, 'r') as f:
        reader = csv.DictReader(f)
        sample_species = set(row['id'] for row in reader)
    
    # Check for missing species
    missing_species = sample_species - tree_species
    
    # Write results to log file
    with open("species_tree_check.log", "w") as log:
        if missing_species:
            log.write("ERROR: The following species from the samplesheet are missing in the species tree:\n")
            for species in missing_species:
                log.write(f"- {species}\n")
            log.write("\nPipeline execution should be stopped due to missing species in the tree.\n")
        else:
            log.write("SUCCESS: All species from the samplesheet are present in the species tree.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check if all species in samplesheet are present in the species tree.")
    parser.add_argument("samplesheet", help="Path to the samplesheet CSV file")
    parser.add_argument("species_tree", help="Path to the species tree Newick file")
    args = parser.parse_args()

    check_species_tree(args.samplesheet, args.species_tree)
    
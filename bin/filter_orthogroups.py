import argparse
from Bio import SeqIO
import os
import shutil

def filter_orthogroups(orthogroup_dir, min_sequences):
    os.makedirs("filtered_orthogroups", exist_ok=True)
    os.makedirs("removed_orthogroups", exist_ok=True)

    for filename in os.listdir(orthogroup_dir):
        if filename.endswith(('.fa', '.faa', '.fasta', '.fas')):
            file_path = os.path.join(orthogroup_dir, filename)
            sequences = list(SeqIO.parse(file_path, "fasta"))
            if len(sequences) > min_sequences:
                shutil.copy2(file_path, "filtered_orthogroups")
                print(f"Kept: {filename}")
            else:
                shutil.move(file_path, "removed_orthogroups")
                print(f"Removed: {filename}")

def main():
    parser = argparse.ArgumentParser(description='Filter orthogroups based on minimum sequence count.')
    parser.add_argument('orthogroup_dir', help='Directory containing orthogroup files')
    parser.add_argument('min_sequences', type=int, help='Minimum number of sequences required to keep an orthogroup')
    args = parser.parse_args()

    filter_orthogroups(args.orthogroup_dir, args.min_sequences)

if __name__ == "__main__":
    main()
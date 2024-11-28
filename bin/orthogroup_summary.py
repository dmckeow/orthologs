import re
import csv
import argparse
from collections import defaultdict

def normalize_defline(defline):
    return re.sub(r'[^a-zA-Z0-9_-]', '_', defline)

def modify_fasta_field(fasta):
    fasta = fasta.replace('.domains.fasta', '')
    parts = fasta.split('.')
    if len(parts) >= 2:
        return '.'.join(parts[-2:])
    return fasta

def parse_input_file(file_path, modify_fasta=False):
    data = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            sample, fasta, defline = parts
            if modify_fasta:
                fasta = modify_fasta_field(fasta)
            normalized_defline = normalize_defline(defline)
            data[normalized_defline] = {'sample': sample, 'fasta': fasta, 'original': defline}
    return data

def parse_orthofinder_file(file_path):
    data = defaultdict(lambda: defaultdict(str))
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            og = parts[0]
            for genome in parts[1:]:
                deflines = genome.split(', ')
                for defline in deflines:
                    normalized_defline = normalize_defline(defline)
                    data[normalized_defline]['og'] = og
                    data[normalized_defline]['original'] = defline
    return data

def parse_broccoli_file(file_path):
    data = defaultdict(lambda: defaultdict(str))
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            og = parts[0]
            for genome in parts[1:]:
                deflines = genome.split()
                for defline in deflines:
                    normalized_defline = normalize_defline(defline)
                    data[normalized_defline]['og'] = og
                    data[normalized_defline]['original'] = defline
    return data

def parse_diamond_mcl_file(file_path):
    data = defaultdict(lambda: defaultdict(str))
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            cluster = parts[0]
            deflines = parts[1].split()
            for defline in deflines:
                normalized_defline = normalize_defline(defline)
                data[normalized_defline]['cluster'] = cluster
                data[normalized_defline]['original'] = defline
    return data

def parse_mmseqs_file(file_path):
    data = defaultdict(lambda: defaultdict(str))
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            cluster = parts[0]
            deflines = parts[1].split()
            for defline in deflines:
                normalized_defline = normalize_defline(defline)
                data[normalized_defline]['cluster'] = cluster
                data[normalized_defline]['original'] = defline
    return data

def main():
    parser = argparse.ArgumentParser(description='Combine orthogroup data from multiple sources.')
    parser.add_argument('--input', required=True, help='Path to the SEARCH deflines file OR the deflines file of original fasta if search skipped (deflines/deflines_combined.txt)')
    parser.add_argument('--orthofinder', help='Path to the OrthoFinder file')
    parser.add_argument('--broccoli', help='Path to the Broccoli file')
    parser.add_argument('--diamond_mcl', help='Path to the Diamond MCL file')
    parser.add_argument('--mmseqs', help='Path to the MMseqs2 file')
    parser.add_argument('--output', default='combined_orthogroups.csv', help='Path to the output CSV file')
    parser.add_argument('--search', action='store_true', help='If set, modifies the fasta field format')
    args = parser.parse_args()

    all_data = {}
    fieldnames = ['original_defline', 'sample', 'input']

    input_data = parse_input_file(args.input, modify_fasta=args.search)
    all_data['input'] = input_data

    if args.orthofinder:
        orthofinder_data = parse_orthofinder_file(args.orthofinder)
        all_data['orthofinder'] = orthofinder_data
        fieldnames.append('orthofinder_og')

    if args.broccoli:
        broccoli_data = parse_broccoli_file(args.broccoli)
        all_data['broccoli'] = broccoli_data
        fieldnames.append('broccoli_og')

    if args.diamond_mcl:
        diamond_mcl_data = parse_diamond_mcl_file(args.diamond_mcl)
        all_data['diamond_mcl'] = diamond_mcl_data
        fieldnames.append('diamond_mcl_cluster')

    if args.mmseqs:
        mmseqs_data = parse_mmseqs_file(args.mmseqs)
        all_data['mmseqs'] = mmseqs_data
        fieldnames.append('mmseqs_cluster')

    # Write the combined data to a CSV file
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for normalized_defline, input_data in all_data['input'].items():
            row = {
                'original_defline': input_data['original'],
                'sample': input_data['sample'],
                'input': input_data['fasta']
            }
            if 'orthofinder' in all_data and normalized_defline in all_data['orthofinder']:
                row['orthofinder_og'] = all_data['orthofinder'][normalized_defline]['og']
            if 'broccoli' in all_data and normalized_defline in all_data['broccoli']:
                row['broccoli_og'] = all_data['broccoli'][normalized_defline]['og']
            if 'diamond_mcl' in all_data and normalized_defline in all_data['diamond_mcl']:
                row['diamond_mcl_cluster'] = all_data['diamond_mcl'][normalized_defline]['cluster']
            if 'mmseqs' in all_data and normalized_defline in all_data['mmseqs']:
                row['mmseqs_cluster'] = all_data['mmseqs'][normalized_defline]['cluster']
            writer.writerow(row)

    print(f"Combined data has been written to '{args.output}'")

if __name__ == "__main__":
    main()
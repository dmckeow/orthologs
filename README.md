# Installation
First, install conda and nextflow on your system
```
git clone --recurse-submodules https://github.com/dmckeow/crg-bcaortho.git
```
# Usage
```
# Prepare your input fasta file directory by making symbolic soft links (ln -s) to your protein fastas (one file per genome) and place them all within a single directory
# Running locally with example data provided
nextflow run main.nf -profile local

# Running using SLURM with example data provided
sbatch submit_nf.sh main.nf -profile slurm
```
## Changing parameters, resources, etc
The following files control the paramters for the pipeline:
    * nextflow.config - parameters unlikely to be changed. If you change anything in this it will be the memory and cpu paramaters
        Recommendation - create new profiles in this file (see local, slurm for examples) and load them from the command line to change memory and cpus (see below)
    * params.config - parameters that are likely to be changed frequently, such as arguments for software, etc. It is loaded by the nextflow.config by default
        Recommendation - copy this file with a memorable name and change the parameters within, then provide to command line with -c
You can run the pipeline with a different config.file e.g.:
```
sbatch submit_nf.sh main.nf -profile slurm -c params-fullhmms.config

```

# Notes
### How Broccoli was incorporated in the pipeline

```
# Broccoli added as a git submodule
git submodule add https://github.com/rderelle/Broccoli.git broccoli

# Broccoli submodule committed
git add .gitmodules broccoli
git commit -m "Add Broccoli as a submodule"

# To update the Broccoli submodule
cd broccoli
git pull origin master
cd ..
git add broccoli
git commit -m "Update Broccoli submodule"


# If the broccoli module mysteriously disappears
git config --file=.gitmodules --get-regexp path # check if the submodule still setup
git submodule init
git submodule update
```


### How much resources?
## OrthoFinder

## Broccoli
8 threads & 0.1-0.5 Gb per eukaryotic genome (from Broccoli documentation)


# Running with test dataset

* Target gene families:
    * Homeodomains
    * Ets
    * TFs
    * zf-met
    * zf-c2h2

* Target genomes:
    * **List of genomes to include:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/species_list_annotation.txt
    * **Genomes:** /users/asebe/xgraubove/genomes/data/*long.pep.fasta

* HMM profies:
    * From Grisha (has an extra field for domain name):
        * **gene_family_info:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/results_phylogenies/gene_families_searchinfo.csv
        * **hmm_dir:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/results_phylogenies/hmms
    * From Xavi (larger and newer dataset):
        * **gene_family_info:** /users/asebe/xgraubove/climate-adaptation-corals-ops/results_annotation/data/gene_families_searchinfo.csv
        * **hmm_dir:** /users/asebe/xgraubove/climate-adaptation-corals-ops/results_annotation/data/hmms



# OrthoFinder problem continues...

.command.log

WARNING: Too few hits between species 21 and species 35 to normalise the scores, these hits will be ignored
WARNING: Too few hits between species 35 and species 10 to normalise the scores, these hits will be ignored
WARNING: Too few hits between species 35 and species 21 to normalise the scores, these hits will be ignored

input/OrthoFinder/Results_orthofinder/WorkingDirectory/SpeciesIDs.txt

10: Coeast_long.pep.None.Homeodomains.domains.fasta
21: Mertsp_long.pep.None.Homeodomains.domains.fasta
35: Sros_long.pep.None.Homeodomains.domains.fasta

then add file name to no_overlaps file
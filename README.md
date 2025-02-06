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
### How Broccoli, Possvm was incorporated in the pipeline

```
# Broccoli added as a git submodule
git submodule add https://github.com/rderelle/Broccoli.git broccoli

# Broccoli submodule committed
git add .gitmodules broccoli
git commit -m "Add Broccoli as a submodule"

# possvm added as a git submodule
git submodule add https://github.com/xgrau/possvm-orthology.git possvm

# possvm submodule committed
git add .gitmodules possvm
git commit -m "Add possvm as a submodule"

# To update a submodule
cd broccoli
git pull origin master
cd ..
git add broccoli
git commit -m "Update Broccoli submodule"


# If the a submodule mysteriously disappears
git config --file=.gitmodules --get-regexp path # check if the submodule still setup
git submodule init
git submodule update
```
### Broccoli

* Originally, we had the original Broccoli as a submodule
* However, its speed scales very poorly due to the multithreaded step of building phylomes (DIAMOND and fasttree)
* If these processes could be run as an array of independent jobs, then it would be so much faster
* This must be addressed if we are to use Broccoli in this project
First, I **forked** the repository of broccoli v1.2 https://github.com/rderelle/Broccoli.git to https://github.com/dmckeow/Broccoli.git then:
```
git submodule add https://github.com/dmckeow/Broccoli.git broccoli
git submodule sync
git submodule update --init --recursive --remote
```
See changelog within the submodule for changes I have made

### Novel tree
* Novel was forked to examine its implementation, which is somewhat similar to our goals
* In the end, the only substantial code from it that ended up being useful was their generax modules
First, I **forked** the repository of noveltree v1.0.2 https://github.com/Arcadia-Science/noveltree.git to https://github.com/dmckeow/noveltree.git then:
```
git submodule add https://github.com/dmckeow/noveltree.git noveltree
git submodule sync
git submodule update --init --recursive --remote
```
See changelog within the submodule for changes I have made


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



# New ioo pipeline
Many goals:
* Use modules inspired by noveltree
* restructure subworkflows
  * search is become a separate workflow that is part of an optional prefilter process
```
nextflow run ioo-main.nf -profile local -c ioo-nextflow.config -params-file ioo-params.json
```
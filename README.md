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
Everything is controlled from within the nextflow.config file.  
You can run the pipeline with a different config.file e.g.:
```
nextflow run main.nf -profile local -c path/to/other/nextflow.config

# You can also provide different parameters using a different params file:
nextflow run main.nf -profile local -params-file path/to/other/params.json
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

# To do
## High priority
* Examine the reportho pipeline html report for comparing ortholog calling tools
* Begin the whichortho workflow to generate criteria for classifying the runnability of orthogroups
* Run initortho on real data to assess runtimes and resource requirements
* Need to refine execution reports for initortho - report and timeline show no time and resource information
* Restructure pipeline
    * Main to subworkflow - main just runs the subworkflows - merge current main with init_orthology subworkflow

## Med priority

## Low priority
* Running a successful search that produces no results will cause downstream clustering to fail. Interactively this cancels the main process, but via sbatch this is not a problem
* Change logic to allow search and all to run simultaneously?
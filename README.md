# Installation
First, install conda and nextflow on your system
```
git clone --recurse-submodules https://github.com/dmckeow/crg-bcaortho.git
```
# Usage
```
# You should make symbolic soft links to your protein fastas (one file per genome) and place them all within a single directory as the initial input
# Running locally
nextflow run main.nf -profile local,toy

# Running using SLURM
sbatch submit_nf.sh main.nf -profile slurm
```
## Changing parameters, resources, etc
Everything is controlled from within the nextflow.config file.  
You can run the pipeline with a different config.file e.g.:
```
nextflow run main.nf -profile local -c path/to/other/nextflow.config
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
```
### How much resources?
## OrthoFinder

## Broccoli
8 threads & 0.1-0.5 Gb per eukaryotic genome (from Broccoli documentation)
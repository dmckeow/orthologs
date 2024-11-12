# Installation
```
git clone --recurse-submodules https://github.com/dmckeow/crg-bcaortho.git
```
# Usage
```
# Running locally
nextflow run main.nf -profile local

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
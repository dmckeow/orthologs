#!/bin/bash
#SBATCH --mem 32G
#SBATCH --cpus-per-task 8
#SBATCH --time=03:00:00
#SBATCH -p genoa64
#SBATCH --qos shorter
#SBATCH --output=/users/asebe/dmckeown/projects/crg-bcaortho/logs/slurm-nf.%j.out
#SBATCH --error=/users/asebe/dmckeown/projects/crg-bcaortho/logs/slurm-nf.%j.err

Rscript -e 'library(rmarkdown); rmarkdown::render("summary.Rmd", "html_document")'

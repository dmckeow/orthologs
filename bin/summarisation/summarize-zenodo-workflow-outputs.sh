#!/bin/bash
#SBATCH --mem 16G
#SBATCH --time=03:00:00
#SBATCH -p genoa64
#SBATCH --qos shorter
#SBATCH --output=/users/asebe/dmckeown/projects/crg-bcaortho/logs/slurm-nf.%j.out
#SBATCH --error=/users/asebe/dmckeown/projects/crg-bcaortho/logs/slurm-nf.%j.err

# Run the R Markdown script to produce the interactive HTML:
#Rscript -e 'library(rmarkdown); rmarkdown::render("noveltree-results-summary.Rmd", "html_document")'
# with renv:
Rscript -e 'renv::restore(); library(rmarkdown); rmarkdown::render("noveltree-results-summary.Rmd", "html_document")'

#!bin/bash

# Download workflow outputs:
#wget https://zenodo.org/record/8237421/files/noveltree-results-tsar-eukaryotes-06062023.tar.gz

# Decompress these outputs:
#tar -xzvf noveltree-results-tsar-eukaryotes-06062023.tar.gz

# Clean up:
#rm noveltree-results-tsar-eukaryotes-06062023.tar.gz

# Run the R Markdown script to produce the interactive HTML:
Rscript -e 'library(rmarkdown); rmarkdown::render("noveltree-results-summary.Rmd", "html_document")'
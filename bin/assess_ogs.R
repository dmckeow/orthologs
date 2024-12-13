suppressMessages(library(cogeqc))
suppressMessages(library(tidyverse))

# Inputs needed:
# orthogroup tsvs for all callers:
# orthogroups_deflines.csv, GET_ORTHOGROUP_INFO.orthogroups.deflines
# annotations from interproscan for all proteins

# Command line arguments
args <- if (interactive()) {
  c("results/WhichTest/report/orthogroups_deflines.csv",
    "example_data/samplesheet-WhichTestAnnot.csv")
} else {
  commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    print("Usage: Rscript assess_ogs.R 
    <orthogroups_file> <samplesheet_file> ")
    quit(status = 1)
  }
}

orthogroups_file <- args[1]
samplesheet_file <- args[2]


# Load orthogroups
og <- tryCatch(
  read.csv(orthogroups_file, header = TRUE, stringsAsFactors = FALSE),
  error = function(e) {
    print(paste("Error reading orthogroups file:", e))
    quit(status = 1)
  }
)

# DFunction for parsing differene orthogroup callers results
process_og <- function(data, variable_chosen) {
  # Rename variables
  data$Orthogroup <- data[[variable_chosen]]
  data$Species <- data$sample
  data$Gene <- data$cleaned_defline
  data <- subset(data, !is.na(Orthogroup) & Orthogroup != "")
  data <- data %>%
    distinct(Orthogroup, Species, Gene, .keep_all = TRUE)
  # Select the columns
  output_data <- data.frame(
    Orthogroup = data$Orthogroup,
    Species = data$Species,
    Gene = data$Gene
  )
  # Return the output data frame
  return(output_data)
}

# Load samplesheet
samplesheet <- tryCatch(
  read_csv(samplesheet_file),
  error = function(e) {
    print(paste("Error reading samplesheet file:", e))
    quit(status = 1)
  }
)

# Function to load annotation file
load_annotation <- function(file_path) {
  tryCatch(
    {
      # Read the file, skipping the header if necessary
      df <- read_tsv(file_path, skip = 1, col_names = c("Gene", "Annotation"))
      return(df)
    },
    error = function(e) {
      print(paste("Error reading annotation file:", file_path, "-", e))
      return(NULL)
    }
  )
}


# Load annotations and create a list of data frames
annotation <- list()
for (i in 1:nrow(samplesheet)) {
  sample_id <- samplesheet$id[i]
  annotation_file <- samplesheet$annotation[i]
  df <- load_annotation(annotation_file)
  if (!is.null(df)) {
    annotation[[sample_id]] <- df
  }
}


# Print summary of loaded data
print(paste("Loaded", length(annotation), "annotation files"))
print("Samples:")
print(names(annotation))

og$cleaned_defline <- sub(":[0-9]+-[0-9]+$", "", og$original_defline)

# Identify protein names that do not match between orthologs and annotations
# Remove those with no matches from og and list those that were removed

og_cleaned_and_unmatched <- lapply(names(annotation), function(species) {
  # Extract species-specific genes from annotation
  species_genes <- annotation[[species]]$Gene
  
  # Filter og for the given species
  species_og <- og %>%
    filter(sample == species)
  
  # Find rows in og where cleaned_defline does not match any Gene in annotation
  unmatched_rows <- species_og %>%
    filter(!cleaned_defline %in% species_genes)
  
  # If unmatched rows exist, print them
  if (nrow(unmatched_rows) > 0) {
    cat("Species:", species, "- Unmatched entries:", nrow(unmatched_rows), "\n")
    print(unmatched_rows)
  }
  
  # Remove the unmatched rows from the species_og data
  species_og_cleaned <- species_og %>%
    filter(cleaned_defline %in% species_genes)  # Keep only matched rows
  
  # Return both the cleaned data and unmatched rows
  list(cleaned_data = species_og_cleaned, unmatched_data = unmatched_rows)
})

# Combine all cleaned species data frames
og <- bind_rows(lapply(og_cleaned_and_unmatched, function(x) x$cleaned_data))

# Combine all unmatched rows into a single data frame for review
og_unmatched <- bind_rows(lapply(og_cleaned_and_unmatched, function(x) x$unmatched_data))





# Extract ortholog table for each source tool
og_orthofinder <- process_og(og, "orthofinder_og")
og_broccoli <- process_og(og, "broccoli_og")
og_dmndmcl <- process_og(og, "diamond_mcl_cluster")
og_mmseqs <- process_og(og, "mmseqs_cluster")

# Generate homogeneity by species
cat("Starting assessment for og_orthofinder...\n")
og_assessment_orthofinder <- assess_orthogroups(og_orthofinder, annotation)

cat("Starting assessment for og_broccoli...\n")
og_assessment_broccoli <- assess_orthogroups(og_broccoli, annotation)

cat("Starting assessment for og_dmndmcl...\n")
og_assessment_dmndmcl <- assess_orthogroups(og_dmndmcl, annotation)

cat("Starting assessment for og_mmseqs...\n")
og_assessment_mmseqs <- assess_orthogroups(og_mmseqs, annotation)


#mean(og_assessment_orthofinder$Mean_score)


process_annotations <- function(og_data, annotation) {
  # Process each species
  result_list <- lapply(names(annotation), function(species) {
    # Select the species-specific annotation data
    species_annotation <- annotation[[species]]
    
    # Merge with the og data for that species
    merged_data <- og_data %>%
      filter(Species == species) %>%
      left_join(species_annotation, by = "Gene")
    
    # Identify unmatched rows 
    unmatched_rows <- merged_data %>%
      filter(is.na(Orthogroups))
    
    # Print information only if there are unmatched rows
    if (nrow(unmatched_rows) > 0) {
      cat("Species:", species, "- Unmatched rows:", nrow(unmatched_rows), "\n")
      print(unmatched_rows)
    }
    
    # Select desired columns
    merged_data <- merged_data %>%
      filter(!is.na(Annotation)) %>%
      select(Gene, Orthogroup, Annotation)
    
    merged_data # Return species-specific merged data
  })
  
  # Combine all species data frames
  combined_df <- bind_rows(result_list)
  
  # Calculate homogeneity scores
  homogeneity_scores <- calculate_H(combined_df)
  
  return(homogeneity_scores)
}


# Generate homogeneity scores for each Orthogroup across all species
# not the same as mean score across species from assess_orthorgoups
hscores_allspecies_orthofinder <- process_annotations(og_orthofinder, annotation)
hscores_allspecies_broccoli <- process_annotations(og_broccoli, annotation)
hscores_allspecies_dmndmcl <- process_annotations(og_dmndmcl, annotation)
hscores_allspecies_mmseqs <- process_annotations(og_mmseqs, annotation)


# Visualisation

#    Specific OrthoFinder results import
#stats_dir <- system.file("extdata", package = "cogeqc")
#ortho_stats <- read_orthofinder_stats(stats_dir)

#    Species tree
#data(tree)
#    Four summary plots for OrthoFinder:
#plot_species_tree(tree, stats_list = ortho_stats)
#plot_duplications(ortho_stats)
#plot_genes_in_ogs(ortho_stats)
#plot_species_specific_ogs(ortho_stats)
#    To plot all these four together:
#plot_orthofinder_stats(
#  tree,
#  xlim = c(-0.1, 2),
#  stats_list = ortho_stats
#)

#plot_og_overlap(ortho_stats)

#plot_og_sizes(og)
#plot_og_sizes(og, log = TRUE) # Or with Natural log
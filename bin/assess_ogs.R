suppressMessages(library(cogeqc))
suppressMessages(library(tidyverse))

# Inputs needed:
# orthogroup tsvs for all callers:
# orthogroups_deflines.csv, GET_ORTHOGROUP_INFO.orthogroups.deflines
# annotations from interproscan for all proteins

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  print("Usage: Rscript assess_ogs.R 
  <orthogroups_file> <samplesheet_file> ")
  quit(status = 1)
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

head(og)

# Define the function
process_og <- function(data, variable_chosen) {
  # Rename variables
  data$Orthogroup <- data[[variable_chosen]]
  data$Species <- data$sample
  data$Gene <- data$original_defline
  # Select the columns
  output_data <- data.frame(
    Orthogroup = data$Orthogroup,
    Species = data$Species,
    Gene = data$Gene
  )
  # Return the output data frame
  return(output_data)
}

og_orthofinder <- process_og(og, "orthofinder_og")
head(og_orthofinder)

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


# Load test data
#data(og)
#>     Orthogroup Species      Gene
#> 1 HOM05D000001     Ath AT1G02310
#> 2 HOM05D000001     Ath AT1G03510

#str(og)
#'data.frame':   88934 obs. of  3 variables:
# $ Orthogroup: chr  "HOM05D000001" "HOM05D000001" "HOM05D000001" "HOM05D000001" ...
# $ Species   : chr  "Ath" "Ath" "Ath" "Ath" ...
# $ Gene      : chr  "AT1G02310" "AT1G03510" "AT1G03540" "AT1G04020" ...


##data(interpro_ath)
#>        Gene Annotation
#> 1 AT1G01010  IPR036093
#> 2 AT1G01010  IPR003441

#str(interpro_ath)
#data.frame':   131661 obs. of  2 variables:
# $ Gene      : chr  "AT1G01010" "AT1G01010" "AT1G01010" "AT1G01020" ...
# $ Annotation: chr  "IPR036093" "IPR003441" "IPR036093" "IPR007290" ...


##data(interpro_bol)
#>           Gene Annotation
#> 1 BolC1t00001H  IPR014710
#> 2 BolC1t00001H  IPR018490

# Assess orthogroups (homogeneity score) - compare species

# Create a list of annotation data frames
##annotation <- list(
##  Ath = interpro_ath,
##  Bol = interpro_bol
##)
#str(annotation)
#> List of 2
#>  $ Ath:'data.frame': 131661 obs. of  2 variables:
#Gene      : chr [1:131661] "AT1G01010" "AT1G01010" "AT1G01010" "AT1G01020" ...
#Annotation: chr [1:131661] "IPR036093" "IPR003441" "IPR036093" "IPR007290" ...
str(og_orthofinder)
str(annotation)
og_assessment_orthofinder <- assess_orthogroups(og_orthofinder, annotation)
#>    Orthogroups    Ath_score Bol_score Mean_score Median_score
#> 1 HOM05D000002  0.143273487 0.5167253  0.3299994    0.3299994
#> 2 HOM05D000003  1.006908255        NA  1.0069083    1.0069083

mean(og_assessment_orthofinder$Mean_score)
#> [1] 1.797598


# Assess orthogroups (homogeneity score) - without species
#     Merge og with annotations
result_list <- lapply(names(annotation), function(species) {
  # Select the species-specific annotation data
  species_annotation <- annotation[[species]]
  # Merge with the og data for that species
  merged_data <- og %>%
    filter(Species == species) %>%
    left_join(species_annotation, by = "Gene")

  merged_data <- merged_data %>%
    select(Gene, Orthogroup, Annotation)
  merged_data # This returns it
})

# Combine the individual species data frames into one
orthogroup_df <- bind_rows(result_list)
# DO the homoengeity scores
og_assessment <- calculate_H(orthogroup_df)

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
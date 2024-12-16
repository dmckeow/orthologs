suppressMessages(library(cogeqc))
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(coin))
suppressMessages(library(ggtext))
suppressMessages(library(treeio))
library(here)


source(here("bin", "utils.R"))
# Inputs needed:
# orthogroup tsvs for all callers:
# orthogroups_deflines.csv, GET_ORTHOGROUP_INFO.orthogroups.deflines
# annotations from interproscan for all proteins

# Command line arguments
args <- if (interactive()) {
  c("results/WhichTest/report/orthogroups_deflines.csv",
    "example_data/samplesheet-WhichTestAnnot.csv",
    "results/WhichTest/orthofinder/Comparative_Genomics_Statistics",
    "results/WhichTest/orthofinder/Species_Tree/SpeciesTree_rooted_node_labels.txt"
  )
} else {
  commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    print("Usage: Rscript assess_ogs.R 
    <orthogroups_file> <samplesheet_file> <orthofinder Comparative_Genomics_Statistics path> <OrthoFinder Species_Tree/SpeciesTree_rooted_node_labels.txt>")
    quit(status = 1)
  }
}

orthogroups_file <- args[1]
samplesheet_file <- args[2]
orthofinder_dir <- args[3]
orthofinder_tree <- args[4]


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


# To filter annotation list only by what is in the corresponding orthogroup calling
filter_annotation_by_genes <- function(annotation_list, gene_df, gene_column = "Gene") {
  # Extract unique genes from the gene_df
  unique_genes <- unique(gene_df[[gene_column]])

  # Iterate through each dataframe in the annotation_list
  filtered_annotations <- lapply(annotation_list, function(annotation_df) {
    # Filter the current annotation_df by matching genes
    filtered_annotation <- annotation_df %>%
      filter(Gene %in% unique_genes)
    
    # Return filtered dataframe (empty ones will have 0 rows)
    return(filtered_annotation)
  })
  
  # Remove empty dataframes (those with 0 rows)
  filtered_annotations <- Filter(function(x) nrow(x) > 0, filtered_annotations)

  # Return the filtered list of annotations
  return(filtered_annotations)
}


# Extract ortholog table for each source tool
og_orthofinder <- process_og(og, "orthofinder_og")
og_broccoli <- process_og(og, "broccoli_og")
og_dmndmcl <- process_og(og, "diamond_mcl_cluster")
og_mmseqs <- process_og(og, "mmseqs_cluster")

# Plot og sizes

OUT_sizes_ogs <- patchwork::wrap_plots(
    plot_og_sizes(og_orthofinder) + ggtitle("OrthoFinder"), 
    plot_og_sizes(og_broccoli) + ggtitle("Broccoli") + theme(axis.text.y = element_blank()),
    plot_og_sizes(og_dmndmcl) + ggtitle("DIAMOND + MCL") + theme(axis.text.y = element_blank()),
    plot_og_sizes(og_mmseqs) + ggtitle("MMseqs2") + theme(axis.text.y = element_blank()),
    nrow = 1, ncol = 4
)

OUT_sizes_ogs

####################################
ogs <- list(og_orthofinder, og_broccoli, og_dmndmcl, og_mmseqs)

# Calculate OG sizes for each run
og_sizes <- lapply(ogs, function(x) {
    sizes <- as.matrix(table(x$Orthogroup, x$Species))
    total <- rowSums(sizes)
    
    sizes_df <- data.frame(unclass(sizes))
    sizes_df$Total <- total
    return(sizes_df)
})

# What is the percentage of OGs with >=100 genes? And with >50 genes?
percentage_size <- function(size_df, n = 100) {
    return(sum(size_df$Total > n) / nrow(size_df) * 100)
}

names(og_sizes) <- c("orthofinder", "broccoli", "dmndmcl", "mmseqs")

percentages <- data.frame(
    Mode = names(og_sizes),
    P200 = unlist(lapply(og_sizes, percentage_size, n = 200)),
    P100 = unlist(lapply(og_sizes, percentage_size, n = 100)),
    P50 = unlist(lapply(og_sizes, percentage_size, n = 50)),
    OGs = unlist(lapply(og_sizes, nrow))
)

# Reorder rows from lowest to highest mcl inflation (not needed?)
#orders <- c("orthofinder", "broccoli", "dmndmcl", "mmseqs")
#percentages <- percentages[orders, ]

# Visual exploration
OUT_percent_sizes <- percentages %>%
    tidyr::pivot_longer(cols = !Mode) %>%
    mutate(name = str_replace_all(
        name,
        c(
            "OGs" = "Number of OGs",
            "P200" = "% OGs with >200 genes",
            "P100" = "% OGs with >100 genes",
            "P50" = "% OGs with >50 genes"
        )
    )) %>%
    ggplot(., aes(y = Mode, x = value)) +
    geom_col(aes(fill = Mode), show.legend = FALSE) +
    scale_fill_manual(
        values = c("orthofinder" = "#08519C",
                   "broccoli" = "#6BAED6",
                   "dmndmcl" = "#006D2C",
                   "mmseqs" = "#74C476")
    ) +
    facet_wrap(~name, ncol = 4, scales = "free_x") +
    theme_bw() +
    labs(
        x = "", y = "Tool",
        title = "Relationship between the number of orthogroups and orthogroup size per tool"
    )

OUT_percent_sizes








#################################

# filter annoation by gene list for each orthololog caller
annotation_orthofinder <- filter_annotation_by_genes(annotation, og_orthofinder)
annotation_broccoli <- filter_annotation_by_genes(annotation, og_broccoli)
annotation_dmndmcl <- filter_annotation_by_genes(annotation, og_dmndmcl)
annotation_mmseqs <- filter_annotation_by_genes(annotation, og_mmseqs)


# Generate homogeneity by species
cat("Starting assessment for og_orthofinder...\n")
og_assessment_orthofinder <- assess_orthogroups(og_orthofinder, annotation_orthofinder)

cat("Starting assessment for og_broccoli...\n")
og_assessment_broccoli <- assess_orthogroups(og_broccoli, annotation_broccoli)

cat("Starting assessment for og_dmndmcl...\n")
og_assessment_dmndmcl <- assess_orthogroups(og_dmndmcl, annotation_dmndmcl)

cat("Starting assessment for og_mmseqs...\n")
#og_assessment_mmseqs <- assess_orthogroups(og_mmseqs, annotation_mmseqs)


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
      filter(is.na(Orthogroup))
    
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
hscores_allspecies_orthofinder <- process_annotations(og_orthofinder, annotation_orthofinder)
hscores_allspecies_broccoli <- process_annotations(og_broccoli, annotation_broccoli)
hscores_allspecies_dmndmcl <- process_annotations(og_dmndmcl, annotation_dmndmcl)
#hscores_allspecies_mmseqs <- process_annotations(og_mmseqs, annotation)


H_combined <- bind_rows(
    hscores_allspecies_orthofinder %>% mutate(Source = "orthofinder"),
    hscores_allspecies_broccoli %>% mutate(Source = "broccoli"),
    hscores_allspecies_dmndmcl %>% mutate(Source = "dmndmcl"),
)

# Scale scores to maximum, so that they range from 0 to 1
H_combined$Score_scaled <- H_combined$Score / max(H_combined$Score)
# Force any values below 0 to be 0
H_combined$Score_scaled <- pmax(H_combined$Score_scaled, 0)


# Compare homogeneity scores - all vs all
db_wilcox <- compare(H_combined, "Score_scaled ~ Source")

OUT_db_wilcox_table <- db_wilcox %>%
    filter_comparison() %>%
    knitr::kable(
        caption = "Mann-Whitney U test for differences in orthogroup scores with Wilcoxon effect sizes.",
        digits = 10
    )


    # Comparisons to be made
comps <- list(
    c("orthofinder", "broccoli", "dmndmcl")
)

# Change order of levels according to comparison results
H_combined$Source <- factor(
    H_combined$Source, levels = rev(c(
        "orthofinder", "broccoli", "dmndmcl", "mmseqs"
    ))
)

# Visualize distributions with significant differences highlighted
OUT_distros <- ggviolin(
    H_combined, y = "Score_scaled", x = "Source", 
    orientation = "horiz", trim = TRUE, add = c("boxplot", "mean"), 
    fill = "Source", add.params = list(fill = "white"), palette = "jama"
) +
    ggpubr::stat_compare_means(
        comparisons = comps,
        label = "p.signif",
        method = "wilcox.test"
    ) +
    theme(legend.position = "none") +
    labs(y = "Scaled homogeneity scores", x = "Source of orthogroups",
         title = "Distribution of mean homogeneity scores for orthogroups",
         subtitle = "Scores were calculated across all genomes included per orthogroup") +
    theme(plot.subtitle = ggtext::element_markdown())

OUT_distros

# Other things

#    Specific OrthoFinder results import

ortho_stats <- read_orthofinder_stats(orthofinder_dir)

#    Species tree
#    Four summary plots for OrthoFinder:
tree <- treeio::read.tree(orthofinder_tree)
#plot_species_tree(tree, stats_list = ortho_stats)
#plot_duplications(ortho_stats)
#plot_genes_in_ogs(ortho_stats)
#plot_species_specific_ogs(ortho_stats)
#    To plot all these four together:
OUT_orthofinder_results <- plot_orthofinder_stats(
  tree,
  xlim = c(-0.1, 2),
  stats_list = ortho_stats
)

OUT_orthofinder_results

OUT_overlaps_ogs <- plot_og_overlap(ortho_stats)

OUT_overlaps_ogs


# Reference protein recall?

# OrthoBench?
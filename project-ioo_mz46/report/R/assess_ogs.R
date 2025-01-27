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





##### new version for Rmd


clean_pfam_scan_file <- function(file_path, header) {
  lines <- readLines(file_path)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nchar(lines) > 0] # remove empty lines
  lines <- trimws(lines) # remove trailing leading whitespace
  lines <- gsub("\\s+", ",", lines)
  data <- read.table(text = lines, header = FALSE, sep = ",", stringsAsFactors = FALSE)
  colnames(data) <- header
  return(data)
}

# Path to save the functional_annotations RDS file
orthogroups_data <- "data/orthogroups.rds"
defline_info_data <- "data/defline_info.rds"
functional_annotations_data <- "data/functional_annotations.rds"

# Check if all three RDS files exist
if (file.exists(orthogroups_data) && file.exists(defline_info_data) && file.exists(functional_annotations_data)) {
  cat("RDS files already exist - loading data from them...\n")
  # If all RDS files exist, load them
  orthogroups <- readRDS(orthogroups_data)
  defline_info <- readRDS(defline_info_data)
  functional_annotations <- readRDS(functional_annotations_data)
} else {
  # If any file is missing, proceed with processing (you can add your processing logic here)
  cat("One or more RDS files are missing. Proceeding with data processing...\n")
  
  # Load orthogroups
  orthogroups <- cogeqc::read_orthogroups("/users/asebe/dmckeown/projects/crg-bcaortho/project-ioo_mz46/results/orthofinder_results/orthofinder_mcl/Orthogroups.tsv")

  # Get all defline info 
  file_paths <- list.files("/users/asebe/dmckeown/projects/crg-bcaortho/project-ioo_mz46/results/prefilter/initial/defline_info/", 
                           pattern = "*.csv", full.names = TRUE)
  defline_info <- lapply(file_paths, read.csv)

  # Get functional annotations
  species_list <- unique(orthogroups$Species)
  file_paths <- list.files("/users/asebe/xgraubove/genomes/annotation_functional/", 
                           pattern = "*_long.pep.pfamscan.csv", full.names = TRUE)

  # Filter the file paths based on species in the orthogroups
  filtered_file_paths <- file_paths[sapply(file_paths, function(x) {
    any(sapply(species_list, function(species) grepl(species, x)))
  })]

  pfam_header <- c("seq_id", "alignment_start", "alignment_end", "envelope_start", "envelope_end", 
                    "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end", "hmm_length", 
                    "bit_score", "E_value", "significance", "clan")

  # Read the filtered CSV files
  functional_annotations <- lapply(filtered_file_paths, clean_pfam_scan_file, header = pfam_header)

  # Extract seq_id values from defline_info (assuming parent_seq contains the seq_ids)
  seq_ids_defline <- unique(unlist(lapply(defline_info, function(x) x$parent_seq)))

  # Filter functional_annotations based on matching seq_id in deflines
  functional_annotations <- lapply(functional_annotations, function(annotation_df) {
    filtered_df <- annotation_df[annotation_df$seq_id %in% seq_ids_defline, ]
    filtered_df <- filtered_df[, c("seq_id", "hmm_acc")]
    colnames(filtered_df) <- c("Gene", "Annotation")
    return(filtered_df)
  })

  functional_annotations <- do.call(rbind, functional_annotations)

  # Create a mapping between parent_seq and clean_seq from defline_info
parent_seq_to_clean_seq <- do.call(rbind, lapply(defline_info, function(x) {
  unique(x[, c("parent_seq", "clean_seq")])
}))

# Merge functional_annotations with the clean_seq mapping
functional_annotations <- merge(functional_annotations, parent_seq_to_clean_seq, 
                                by.x = "Gene", by.y = "parent_seq", 
                                all.x = TRUE)

# Replace 'Gene' column (seq_id) with 'clean_seq'
functional_annotations$Gene <- functional_annotations$clean_seq

# Drop the 'parent_seq' and 'clean_seq' columns, as we don't need them anymore
functional_annotations <- functional_annotations[, c("Gene", "Annotation")]


  # Replace Gene annotations with clean version

  # Save the processed functional_annotations to the RDS file
  saveRDS(orthogroups, orthogroups_data)
  saveRDS(defline_info, defline_info_data)
  saveRDS(functional_annotations, functional_annotations_data)
}


og_assessment_orthofinder <- assess_orthogroups(orthogroups, functional_annotations)






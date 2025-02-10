


##############################
# Plotting with cogeqc

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

prepare_all_data <- function(
  results_dir = NULL,
  samplesheet_fpath = NULL,
  orthogroups_path = NULL,
  functional_annotations = NULL,
  OG_source = NULL  # Add new parameter for OG_source
) {

  # Read the samplesheet
  samplesheet <- read.delim(samplesheet_fpath, sep = ",")
  
  # Read the orthogroups data
  orthogroups <- cogeqc::read_orthogroups(orthogroups_path)
  
  # Get all defline info 
  file_paths <- list.files(paste0(results_dir, "prefilter/initial/defline_info/"), 
                            pattern = "*.csv", full.names = TRUE)
  defline_info <- lapply(file_paths, read.csv)
  defline_info <- do.call(rbind, defline_info)

  # Link functional annotations to sequence ids
  merged_data <- defline_info %>%
    left_join(functional_annotations, by = c("parent_seq" = "seq_id"))
  
  # Merge with samplesheet and rename columns
  merged_data <- merged_data %>%
    left_join(samplesheet, by = c("id" = "id")) %>%
    rename(
      supergroup = taxonomy
    ) %>%
    select(-fasta)

  # Merge with orthogroups data
  merged_data <- merged_data %>%
    left_join(orthogroups, by = c("clean_seq" = "Gene")) %>%
    select(-Species)
  
  # Add the new OG_source column to the merged data
  merged_data <- merged_data %>%
    mutate(OG_source = OG_source)  # Fill the column with the passed value
  
  # Return the merged data with the new OG_source column
  return(merged_data)
}


prepare_cogeqc_inputs <- function(
    input_df = NULL,
    species_col = NULL,
    gene_col = NULL,
    annot_col = NULL,
    annot_filter_col = NULL,
    annot_filter_value = NULL
  ) {
    orthogroups <- input_df %>%
      filter(!is.na(Orthogroup)) %>%
      select(
          Orthogroup,
          {{species_col}},
          {{gene_col}}
        ) %>%
      rename(
        Species = {{species_col}},
        Gene = {{gene_col}}
      ) %>%
      unique()

      annotations <- input_df %>%
    filter(!is.na({{annot_col}})) %>%  # Filter out rows with NA in annot_col
    filter({{annot_filter_col}} == annot_filter_value) %>%  # Filter based on annot_filter_col and value
    select(
      {{species_col}},
      {{gene_col}},
      {{annot_col}}
    ) %>%
    rename(
      Species = {{species_col}},
      Gene = {{gene_col}},
      Annotation = {{annot_col}}
    ) %>%
    unique()
  # split into lists by species
  annotations_split <- annotations %>%
    split(annotations$Species) %>%
    lapply(function(df) {
      df %>%
        select(-Species)  # Remove the 'Species' column from each split dataframe
    })
  
  # Return both dataframes (orthogroups and annotations_list)
  return(list(orthogroups = orthogroups, annotations = annotations_split))
  }


get_cogeqc_hscores <- function(og_data, annotation) {
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


compare_homogeneity_scores <- function(..., comps = NULL) {
  # Collect all dataframes passed as arguments
  df_list <- list(...)
  
  # Check if there are source names passed, if not, use the names of the dataframes
  source_names_raw <- sapply(substitute(list(...))[-1], deparse)
  source_names <- gsub('.*\\"(.*?)\\".*', '\\1', source_names_raw)
  
  # Add 'Source' column to each dataframe and combine them
  H_combined <- bind_rows(
    lapply(seq_along(df_list), function(i) {
      df_list[[i]] %>% 
        mutate(Source = source_names[i])
    })
  )

  # Scale scores to maximum, so that they range from 0 to 1
  H_combined$Score_scaled <- H_combined$Score / max(H_combined$Score)
  
  # Force any values below 0 to be 0
  H_combined$Score_scaled <- pmax(H_combined$Score_scaled, 0)
  
  # Visualize distributions with significant differences highlighted
  OUT_distros <- ggviolin(
    H_combined, y = "Score_scaled", x = "Source", 
    orientation = "horiz", trim = TRUE, add = c("boxplot", "mean"), 
    fill = "Source", add.params = list(fill = "white"), palette = "jama"
  ) +
    theme(legend.position = "none") +
    labs(y = "Scaled homogeneity scores", x = "Source of orthogroups",
         title = "Distribution of mean homogeneity scores for orthogroups",
         subtitle = "Scores were calculated across all genomes included per orthogroup") +
    theme(plot.subtitle = ggtext::element_markdown())
  
  # Return the table and plot as a list
  return(list(
    og_hscores_scaled_all = H_combined,
    og_hscores_scaled_violin = OUT_distros
  ))
}


compare_wilcox <- function(data, form, ref = NULL) {
    # Wilcoxon test - greater and less alternatives
    wilcoxtest_greater <- tibble::as_tibble(data) %>%
        rstatix::wilcox_test(
            formula(form), p.adjust.method = "BH", ref.group = ref,
            alternative = "greater"
        )
    pg <- ifelse("p.adj" %in% names(wilcoxtest_greater), "p.adj", "p")
    wilcoxtest_greater <- wilcoxtest_greater %>% dplyr::select(
            group1, group2, n1, n2, padj_greater = all_of(pg)
        )
    
    wilcoxtest_less <- tibble::as_tibble(data) %>%
        rstatix::wilcox_test(
            formula(form), p.adjust.method = "BH", ref.group = ref,
            alternative = "less"
        )
    pl <- ifelse("p.adj" %in% names(wilcoxtest_less), "p.adj", "p")
    wilcoxtest_less <- wilcoxtest_less %>% dplyr::select(
        group1, group2, n1, n2, padj_less = all_of(pl)
    )
    
    wilcox_summary <- dplyr::inner_join(wilcoxtest_greater, wilcoxtest_less) %>%
        dplyr::mutate(padj_interpretation = dplyr::case_when(
            padj_less < 0.05 ~ "less",
            padj_greater < 0.05 ~ "greater",
            TRUE ~ "ns"
        ))
    
    # Effect sizes for Wilcoxon test - greater and less alternatives
    
    effsize <- tibble::as_tibble(data) %>%
        rstatix::wilcox_effsize(
            formula(form), ref.group = ref,
        ) %>%
        dplyr::select(
            group1, group2, effsize, magnitude
        )
    
    
    result <- as.data.frame(inner_join(wilcox_summary, effsize))
        
    return(result)
}


filter_comparison_wilcox <- function(compare_output) {
    
    filtered_df <- compare_output |>
        dplyr::filter(padj_interpretation != "ns") |>
        mutate(padj = case_when(
            padj_interpretation == "greater" ~ padj_greater,
            padj_interpretation == "less" ~ padj_less
        )) |>
        dplyr::select(group1, group2, n1, n2, padj, effsize, magnitude)
    
    return(filtered_df)
}


## Jaccard pairwise
  # Function to calculate Jaccard similarity between two sets of sequences
jaccard_similarity <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersection_size / union_size)
}

# Function to calculate pairwise Jaccard similarity between two orthogroup lists
calculate_jaccard_similarity <- function(ogs1, ogs2) {
  # Get the orthogroup IDs (names of the list)
  og_ids1 <- names(ogs1)
  og_ids2 <- names(ogs2)
  
  # Initialize the similarity matrix
  jaccard_matrix <- matrix(0, nrow = length(og_ids1), ncol = length(og_ids2))
  rownames(jaccard_matrix) <- og_ids1
  colnames(jaccard_matrix) <- og_ids2
  
  # Loop through all pairs of orthogroups and calculate Jaccard similarity
  for (i in 1:length(og_ids1)) {
    for (j in 1:length(og_ids2)) {
      # Get the sets of sequences
      set1 <- ogs1[[og_ids1[i]]]
      set2 <- ogs2[[og_ids2[j]]]
      
      # Calculate Jaccard similarity and store it in the matrix
      jaccard_matrix[i, j] <- jaccard_similarity(set1, set2)
    }
  }
  
  # Return the Jaccard similarity matrix
  return(jaccard_matrix)
}


calculate_jaccard_metrics <- function(data, tool1_name, tool2_name) {
  data_matrix <- as.matrix(data)
  all_values <- as.vector(data_matrix)
  total_ogs <- nrow(data_matrix) + ncol(data_matrix) - 1

  rows_with_one <- sum(apply(data_matrix, 1, function(row) any(row == 1, na.rm = TRUE)))
  cols_with_one <- sum(apply(data_matrix, 1, function(col) any(col == 1, na.rm = TRUE)))
  total_ogs_with_one = rows_with_one + cols_with_one - 1
  portion_ones <- (total_ogs_with_one / total_ogs)

  max_jaccard_per_column <- apply(data_matrix, 2, max, na.rm = TRUE)
  mean_best_jaccard_index <- mean(max_jaccard_per_column, na.rm = TRUE)

  # Create a data frame with the results
  result <- data.frame(
    tool1 = tool1_name,
    tool2 = tool2_name,
    Portion_of_OGs_identical_pw = portion_ones,
    Mean_Best_Jaccard = mean_best_jaccard_index
  )
  
  return(result)
}


create_presence_matrix <- function(df, group_col, value_col) {
  df %>%
    select({{ value_col }}, {{ group_col }}) %>%  # Select the specified columns
    distinct() %>%  # Remove duplicate combinations
    mutate(Presence = 1) %>%  # Mark presence with 1
    spread(key = {{ group_col }}, value = Presence, fill = 0) %>%  # Convert to wide format
    rename(Name = {{ value_col }}) %>%
    filter(!is.na(Name))
}



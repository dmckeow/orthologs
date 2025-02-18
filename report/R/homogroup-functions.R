


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


get_cogeqc_hscores <- function(og_data, annotation, max_og_size = 200) {
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
  homogeneity_scores <- calculate_H(
                          combined_df,
                          max_size = {{max_og_size}}
                          )

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
  H_combined$Score_scaled_all <- H_combined$Score / max(H_combined$Score)

  H_combined <- H_combined %>%
    group_by(Source) %>%
    mutate(Score_scaled_self = Score / max(Score)) %>%
    ungroup()
  
  # Force any values below 0 to be 0
  H_combined$Score_scaled_self <- pmax(H_combined$Score_scaled_self, 0)
  H_combined$Score_scaled_all <- pmax(H_combined$Score_scaled_all, 0)
  
  # Return the table and plot as a list
  return(
    og_hscores_scaled_all = H_combined
  )
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
  total_ogs <- nrow(data_matrix) + ncol(data_matrix) - 1

  rows_with_one <- sum(apply(data_matrix, 1, function(row) any(row == 1, na.rm = TRUE)))
  cols_with_one <- sum(apply(data_matrix, 1, function(col) any(col == 1, na.rm = TRUE)))
  total_ogs_with_one = rows_with_one + cols_with_one - 1
  portion_ones <- (total_ogs_with_one / total_ogs)

  max_jaccard_per_col <- apply(data_matrix, 2, max, na.rm = TRUE)

  mean_best_jaccard_index <- mean(max_jaccard_per_col, na.rm = TRUE)

  # Create a data frame with the results
  jaccard_by_og_source <- data.frame(
    tool1 = tool1_name,
    tool2 = tool2_name,
    Portion_of_OGs_identical_pw = portion_ones,
    Mean_Best_Jaccard = mean_best_jaccard_index
  )

  return(jaccard_by_og_source)
}

calculate_jaccard_metrics_per_og <- function(data, tool1_name, tool2_name) {
  data_matrix <- as.matrix(data)

  # Calculate the maximum Jaccard value for each orthogroup (column)
  max_jaccard_per_col <- apply(data_matrix, 2, max, na.rm = TRUE)
  
  # Find the row indices corresponding to the maximum Jaccard value for each column
  max_jaccard_row_indices <- apply(data_matrix, 2, which.max)
  
  # Get the orthogroup names corresponding to the maximum values
  max_jaccard_og_names <- rownames(data_matrix)[max_jaccard_row_indices]

  # Calculate the mean Jaccard for each orthogroup (column)
  mean_jaccard_per_col <- apply(data_matrix, 2, mean, na.rm = TRUE)
  
  # Create a data frame with the necessary information
  jaccard_by_og <- data.frame(
    Orthogroup = colnames(data_matrix),
    OG_source = tool2_name,
    vs_Orthogroup = max_jaccard_og_names,
    vs_OG_source = tool1_name,
    Best_Jaccard = max_jaccard_per_col,
    Mean_Jaccard = mean_jaccard_per_col
  )
  
  # Return the data frame
  return(jaccard_by_og)
}




create_presence_matrix <- function(df, group_col, value_col) {
  df %>%
    select({{ value_col }}, {{ group_col }}) %>%
    distinct() %>%
    mutate(Presence = TRUE) %>%
    spread(key = {{ group_col }}, value = Presence, fill = FALSE) %>%  # to wide
    filter(!is.na({{ value_col }}))
}


get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

count_mode <- function(x, mode_value) {
  sum(x == mode_value)
}




# This does the opposite of separate_longer_delim
collapse_col <- function(df, col_name, group_by_cols, separator = ",", unique_values = FALSE) {
  # Convert col_name and group_by_cols to symbols
  col_name <- rlang::ensym(col_name)  # Convert col_name to a symbol
  group_by_cols <- rlang::syms(group_by_cols)  # Convert group_by_cols to symbols

  # Group by multiple columns using `!!!` to splice the symbols correctly
  df %>%
    group_by(!!!group_by_cols) %>%  # Group by the specified columns
    mutate(
      !!col_name := str_c(
        if (unique_values) {
          unique(!!col_name)  # Apply unique if requested
        } else {
          !!col_name  # Otherwise, use the original column
        }, 
        collapse = separator
      )
    ) %>% 
    ungroup()
}



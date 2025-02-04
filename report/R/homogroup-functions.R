




####################################



###########################################################
# OrthoFinder summary functions
###########################################################
# Function to read in and get per-species gene copy number
# for each orthogroup
# @
get_per_spp_og_counts <- 
  function(results_dir = NULL, out_dir = NULL, gene_count_tsv = NULL){
    
    # read in per-species gene copy number per gene family
    og_counts <- 
      read.table(gene_count_tsv, 
                 header = T, check.names = F)
    colnames(og_counts) <- gsub("\\..*", "", colnames(og_counts))
    
    # And calculate the number of species in each gene family
    og_counts$NumSpecies <- 
      rowSums(og_counts[,-c(1,ncol(og_counts))] > 0)
    
    # Write out to file if an output directory is provided
    if(!is.null(out_dir)){
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file
      write.table(og_counts, file = paste0(out_dir, "og_copy_num_per_spp.tsv"), 
                  sep = "\t", quote = F, row.names = F, col.names = T)
    }
    return(og_counts)
  }


# Function to summarize gene family content, both with respect to number of 
# species included in each gene family, and the per-species mean copy number
# @
summ_genefam_composition <- 
  function(results_dir = NULL, show_plots = F, out_dir = NULL){
    # Read in the gene copy number and species count for each gene family
    res <- read.table(paste0(results_dir, "filtered_orthogroups/all_ogs_counts.csv"), sep = ",", header = T)
    
    # Plot the mean species gene copy number against the number of species 
    # in each gene family
    copynum_v_numspp <- 
      ggplot(data = res, aes(y = mean_copy_num, x = num_spp)) + 
      geom_point(position = position_jitter(width = 0.2),
                 alpha = 0.25, size = 0.75) + 
      theme_classic(base_size = 14) + 
      scale_y_log10() + 
      annotation_logticks(sides = 'l') + 
      xlab("# Species in gene family") + 
      ylab("Mean per-species\ngene copy #")
    
    # And plot the histogram of species counts across all gene families
    fam_sppcount_dist <- 
      ggplot(data = res, aes(x = num_spp)) + 
      geom_histogram(fill = 'lightgrey', color = 'black') + 
      theme_classic(base_size = 14) + 
      scale_y_log10() + 
      annotation_logticks(sides = 'l') + 
      xlab("# Species in gene family") + 
      ylab("# Gene families")
    
    # Combine them
    genefam_composition_plt <- 
      suppressMessages(plot_grid(copynum_v_numspp, fam_sppcount_dist, ncol = 2))
    
    # print the plots if that is requested
    if(show_plots == T){print(genefam_composition_plt)}
    
    # And save
    ggsave(genefam_composition_plt, filename = paste0(out_dir, "genefamily_compostion.png"), 
           height = 6, width = 12, dpi = 600)
    ggsave(genefam_composition_plt, filename = paste0(out_dir, "genefamily_compostion.pdf"), 
           height = 6, width = 12)
  }

# Summarize and plot orthogroup statistics per species
# @
get_per_spp_ofinder_stats <-
  function(tree = NULL, tree_fpath = NULL,
           results_dir = NULL, show_plots = T,
           samplesheet = NULL, samplesheet_fpath = NULL, 
           grp_name = "Group", tip_grp_cols = NULL, 
           count_cols = arcadia_cividis,
           prop_cols = arcadia_viridis, outgroup = NULL,
           out_dir = NULL, og_stats_perspp_fpath = NULL){
    
    if(is.null(tree)){tree <- read.tree(tree_fpath)}
    if(is.null(tree) & is.null(tree_fpath)){
      print("Error! Provide tree object or filepath to tree")}
    if(is.null(samplesheet)){samplesheet <- read.delim(samplesheet_fpath, sep = ",")}
    if(is.null(samplesheet) & is.null(samplesheet_fpath)){
      print("Error! Provide samplesheet object or filepath to samplesheet")}
    
    ###########################################
    # Orthogroup statistics per species       #
    ###########################################
    species <- tree$tip.label[order(tree$tip.label)]
    
    # First read in the overall stats per species and write out as a table
    og_stats_perspp <- 
      read_tsv(og_stats_perspp_fpath, skip = 1, n_max = 10,
               col_types = c("f", rep("n", length(species))),
               col_names=c("Statistic", species))
    
    # Transpose and then plot stats as a heatmap alongside the the species tree
    og_stats_perspp <- og_stats_perspp %>%
      gather(key = species, value = value, 2:ncol(og_stats_perspp)) %>% 
      spread(key = names(og_stats_perspp)[1],value = "value")
    # Check the above for use of "spread_"
    og_stats_perspp <- 
      data.frame(og_stats_perspp[order(match(og_stats_perspp$species, 
                                             tree$tip.label)),])
    rownames(og_stats_perspp) <- tree$tip.label
    
    # Pull out and rename relevant stats
    plot_dat_nums <- og_stats_perspp[c(2,3,4,7,9,10)]
    plot_dat_props <- og_stats_perspp[c(5,6,8,11)]
    colnames(plot_dat_nums) <- 
      c("# genes", "# genes in OGs", "# unassigned genes", "# OGs containing spp.",
        "# spp. specific OGs", "# genes in spp. specific OGs")
    colnames(plot_dat_props) <- 
      c("% genes in OGs", "% unassigned genes", "% OGs containing spp.", 
        "% genes in spp. specific OGs")
    
    # Now prepare to plot
    # Get metadata for plotting
    taxa_tibble <- tibble(
      tip = samplesheet$id,
      grp = samplesheet$taxonomy)
    
    # Root if using the asteroid tree and a specified outgroup
    if(!is.null(outgroup)){
      # find the most recent common ancestor of the outgroup
      mrca_out <- getMRCA(tree, tip = outgroup)
      # find the parent of this MRCA node
      parent_node <- which(tree$edge[, 2] == mrca_out)
      # the parent edge number is the second column value at the corresponding row
      parent_edge_number <- tree$edge[parent_node, 2]
      # root the tree at the parent node
      tree <- root(tree, edgelabel = T, node = parent_edge_number, position = 0.25)
      
      # Make the basic tree plot
      base_p <- suppressMessages(
        ggtree(tree, size = 0.75, branch.length = 'none') %<+% taxa_tibble+ 
          geom_tiplab(size = 7, as_ylab = T))
    }else{
      # Make the basic tree plot
      base_p <- suppressMessages(
        ggtree(tree, size = 0.75) %<+% taxa_tibble+ 
          geom_tiplab(size = 7, as_ylab = T))
    }
    
    # Now one with the tip labels
    tree_p <- suppressMessages(
      base_p + 
        geom_tippoint(aes(fill = grp), pch = 22, 
                      size = 3, color = "black") + 
        scale_fill_manual(values = SUPERGROUP_COLS) + 
        guides(fill = guide_legend(title = grp_name, nrow = 2,
                                   title.position="top", title.hjust = 0)) + 
        theme(legend.position = "right", legend.justification = "left") + 
        new_scale_fill()) 
    
    # Combine with heatmap of per-species orthogroup count stats
    count_p <- suppressMessages(
      gheatmap(base_p, plot_dat_nums, offset = 0.0, 
               font.size = 2.5, width = 0.6, hjust = 1,
               colnames_angle = 25, 
               colnames_offset_y = 0.25) + 
        scale_fill_gradientn(colors = count_cols$color_dict,
                             values = count_cols$values, 
                             trans = "log1p") + 
        theme(legend.position = "none"))
    
    # Get legends...
    # For the tree & tip labels
    leg1 <- 
      get_legend(tree_p + 
                   theme(legend.text = element_text(size = 8),
                         legend.title = element_text(size = 10),
                         legend.key.width = unit(1,"cm")))
    
    # For the count stats
    leg2 <- get_legend(
      count_p + 
        labs(fill = "Count") + 
        guides(fill = guide_colorbar(title.position="top", 
                                     title.hjust = 0)) +
        theme(legend.position = "right", legend.justification = "left",
              legend.key.width = unit(1,"cm"), 
              plot.margin = margin(10, 10, 15,10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10)))
    
    # And finally, proportional stats
    leg3 <- suppressMessages(get_legend(
      gheatmap(base_p, plot_dat_props, offset = 0.5, 
               font.size = 2.5, width = 0.6, hjust = 1,
               colnames_angle = 25, 
               colnames_offset_y = 0.25) + 
        scale_fill_gradientn(colors = prop_cols$color_dict,
                             values = prop_cols$values, 
                             trans = "log1p", 
                             labels = c(0, 5, 10, 25, 50, 75, 100),
                             breaks = c(0, 5, 10, 25, 50, 75, 100)) +
        theme(legend.position = "right", legend.justification = "left",
              legend.key.width = unit(1,"cm"), 
              plot.margin = margin(10, 10, 15,10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10)) + 
        guides(fill = guide_colorbar(title.position="top", 
                                     title.hjust = 0)) +
        labs(fill = "Percent")))
    
    og_spp_stats_p <- suppressMessages(
      gheatmap(tree_p, plot_dat_nums, offset = 0.0, 
               font.size = 2.5, width = 0.6, hjust = 1,
               colnames_angle = 25, 
               colnames_offset_y = 0.25) + 
        scale_fill_gradientn(colors = count_cols$color_dict,
                             values = count_cols$values, 
                             trans = "log1p") + 
        theme(legend.position = "none") +
        new_scale_fill())
    
    og_spp_stats_p <- suppressMessages(
      gheatmap(og_spp_stats_p, plot_dat_props, offset = 3, 
               font.size = 2.5, width = 0.4, hjust = 1,
               colnames_angle = 25, 
               colnames_offset_y = 0.25) + 
        scale_fill_gradientn(colors = prop_cols$color_dict,
                             values = prop_cols$values, 
                             trans = "log1p", 
                             labels = c(0, 5, 10, 25, 50, 75, 100),
                             breaks = c(0, 5, 10, 25, 50, 75, 100)) + 
        theme(legend.position = "none",
              plot.margin = margin(-10, 10, 35,10)) + 
        coord_cartesian(clip = "off"))
    
    legs <- 
      plot_grid(leg1, NULL, leg2, leg3, ncol = 2)
    
    og_spp_stats_p_final <- 
      plot_grid(legs, og_spp_stats_p, ncol = 1, nrow = 2,
                rel_heights = c(0.5, 1))
    
    # Plot out if requested
    if(show_plots == T){print(og_spp_stats_p_final)}
    
    # Write out to file if an output directory is provided
    if(!is.null(out_dir)){
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file c( both the tabular summary and the figure)
      write.table(og_stats_perspp, file = paste0(out_dir, "per_spp_og_stats.tsv"),
                  row.names = F, col.names = T, quote = F)
      
      # That was involved.... Save this to a pdf
      ggsave(og_spp_stats_p_final, 
             filename = paste0(out_dir, "species_tree_gene_og_stat_heatmaps.pdf"),
             height = 8, width = 8)
      ggsave(og_spp_stats_p_final, 
             filename = paste0(out_dir, "species_tree_gene_og_stat_heatmaps.png"),
             height = 8, width = 8, dpi = 600)
    }
    
    # Return both the dataframe and plot of per-species stats (in the inevitable 
    # case that plotting size needs to tweaked for the given dataset)
    return(list(og_stats_per_spp_table = og_stats_perspp,
                og_stats_per_spp_plt = og_spp_stats_p_final))
  }


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

prepare_all_data <- function(samplesheet_fpath = NULL) {
  samplesheet <- read.delim(samplesheet_fpath, sep = ",")
  # Path to save the functional_annotations RDS file
  df_main_rds <- "data/df_main.rds"

  # Check if all three RDS files exist
  if (file.exists(df_main_rds)) {
    cat(paste("File", df_main_rds, "exists - loading data from it...\n"))
    # If all RDS files exist, load them
    merged_data <- readRDS(df_main_rds)
    return(merged_data)

  } else {
    # If any file is missing, proceed with processing (you can add your processing logic here)
    cat(paste("File", df_main_rds, "does not exist - generating it...\n"))
    
    # Load orthogroups
    orthogroups_of <- cogeqc::read_orthogroups("../project-ioo_mz155/results/orthofinder_results/orthofinder_mcl/Orthogroups.tsv")
    orthogroups_br <- cogeqc::read_orthogroups("../project-ioo_mz155/results/broccoli_results/broccoli/Orthogroups.tsv")

    # Get all defline info 
    file_paths <- list.files("../project-ioo_mz155/results/prefilter/initial/defline_info/", 
                            pattern = "*.csv", full.names = TRUE)
    defline_info <- lapply(file_paths, read.csv)
    defline_info <- do.call(rbind, defline_info)

    # Get functional annotations
    species_list <- unique(defline_info$id)
    file_paths <- list.files("../../../../xgraubove/genomes/annotation_functional/", 
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
    functional_annotations <- do.call(rbind, functional_annotations)


  # We are linking functional annotations to seq ids, 
  # so there will be multiple duplicated seqid lines due to multiple Pfam hits, domain, family, alignment range, etc
  merged_data <- defline_info %>%
    left_join(functional_annotations, by = c("parent_seq" = "seq_id"))
  
  merged_data <- merged_data %>%
    left_join(samplesheet, by = c("id" = "id")) %>%
    rename(
      supergroup = taxonomy
    ) %>%
    select(-fasta)

  merged_data <- merged_data %>%
    left_join(orthogroups_of, by = c("clean_seq" = "Gene")) %>%
    rename(
      Orthogroup_Orthofinder = Orthogroup
    ) %>%
    select(-Species)

  merged_data <- merged_data %>%
    left_join(orthogroups_br, by = c("clean_seq" = "Gene")) %>%
    rename(
      Orthogroup_Broccoli = Orthogroup
    ) %>%
    select(-Species)

    saveRDS(merged_data, df_main_rds)

    return(merged_data)
  }

}

prepare_cogeqc_inputs <- function(
    input_df = NULL,
    orthogroup_col = NULL,
    species_col = NULL,
    gene_col = NULL,
    annot_col = NULL,
    annot_filter_col = NULL,
    annot_filter_value = NULL
  ) {
    orthogroups <- input_df %>%
      filter(!is.na({{orthogroup_col}})) %>%
      select(
          {{orthogroup_col}},
          {{species_col}},
          {{gene_col}}
        ) %>%
      rename(
        Orthogroup = {{orthogroup_col}},
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
  source_names <- sapply(substitute(list(...))[-1], deparse)
  
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
  
  # Perform Wilcoxon test to compare homogeneity scores between groups
  db_wilcox <- compare_wilcox(H_combined, "Score_scaled ~ Source")
  
  # Prepare table for output with significant comparisons
  OUT_db_wilcox_table <- db_wilcox %>%
    filter_comparison_wilcox() %>%
    knitr::kable(
      caption = "Mann-Whitney U test for differences in orthogroup scores with Wilcoxon effect sizes.",
      digits = 10
    )
  
  # Change order of levels in Source factor for plotting
  H_combined$Source <- factor(
    H_combined$Source, levels = rev(source_names)
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
  
  # Return the table and plot as a list
  return(list(
    og_hscores_scaled_all = H_combined,
    og_hscores_scaled_all_wilcox = OUT_db_wilcox_table,
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
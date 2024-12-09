library(cogeqc)
library(dplyr)

# Load test data
data(og)
#>     Orthogroup Species      Gene
#> 1 HOM05D000001     Ath AT1G02310
#> 2 HOM05D000001     Ath AT1G03510

data(interpro_ath)
#>        Gene Annotation
#> 1 AT1G01010  IPR036093
#> 2 AT1G01010  IPR003441

data(interpro_bol)
#>           Gene Annotation
#> 1 BolC1t00001H  IPR014710
#> 2 BolC1t00001H  IPR018490

# Assess orthogroups (homogeneity score) - compare species

# Create a list of annotation data frames
annotation <- list(
  Ath = interpro_ath,
  Bol = interpro_bol
)
#> List of 2
#>  $ Ath:'data.frame': 131661 obs. of  2 variables:
#Gene      : chr [1:131661] "AT1G01010" "AT1G01010" "AT1G01010" "AT1G01020" ...
#Annotation: chr [1:131661] "IPR036093" "IPR003441" "IPR036093" "IPR007290" ...

og_assessment <- assess_orthogroups(og, annotation)
#>    Orthogroups    Ath_score Bol_score Mean_score Median_score
#> 1 HOM05D000002  0.143273487 0.5167253  0.3299994    0.3299994
#> 2 HOM05D000003  1.006908255        NA  1.0069083    1.0069083

mean(og_assessment$Mean_score)
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
stats_dir <- system.file("extdata", package = "cogeqc")
ortho_stats <- read_orthofinder_stats(stats_dir)

#    Species tree
data(tree)
#    Four summary plots for OrthoFinder:
plot_species_tree(tree, stats_list = ortho_stats)
plot_duplications(ortho_stats)
plot_genes_in_ogs(ortho_stats)
plot_species_specific_ogs(ortho_stats)
#    To plot all these four together:
plot_orthofinder_stats(
  tree,
  xlim = c(-0.1, 2),
  stats_list = ortho_stats
)

plot_og_overlap(ortho_stats)

plot_og_sizes(og)
plot_og_sizes(og, log = TRUE) # Or with Natural log
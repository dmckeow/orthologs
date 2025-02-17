library(tidyverse)
library(ggtree)
library(ggnewscale)
library(parallel)
library(reshape2)
library(cowplot)
library(phytools)
library(ape)
library(ggforce)
library(plotly)
library(cogeqc)
library(DiagrammeR)
library(purrr)
library(ggpubr)
library(rstatix)
library(here)
library(coin)
library(ggtext)
library(treeio)
library(pheatmap)
library(ComplexUpset)
library(plotly)
library(viridis)
library(RColorBrewer)
library(patchwork)

suppressPackageStartupMessages(suppressWarnings(source('R/homogroup-functions.R')))

################################################################################
# INPUTS
################################################################################

################################################################################
# Pipeline run inputs
################################################################################
RESULT_DIR  <-  "/no_backup/asebe/dmckeown/test_br_array_vs_br/results/"
SAMPLESHEET <- "../inputs/samplesheet1-10.csv"
PIPELINE_RUN_NAME <- "test_br_array_vs_br" # we will use this later to identify runs if comparing them this way

#SPECIES_TREE <- paste0(RESULT_DIR, 'speciesrax/species_trees/inferred_species_tree.newick') # internal speciesrax tree will be here
SPECIES_TREE <- '../../../../gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/data_annotation/species_tree.newick'

################################################################################
# Orthogroup caller specific stuff (subworkflow i.e. orthofinder, broccoli)
################################################################################

OGS_TSV_PATH_OF <- paste0(RESULT_DIR, "orthofinder_results/orthofinder_mcl/Orthogroups.tsv")
OGS_TSV_PATH_BR <- paste0(RESULT_DIR, "broccoli_array_results/broccoli/Orthogroups.tsv")

LONG_PEP_PFAMSCAN_CSVS <- list.files("../../../../xgraubove/genomes/annotation_functional/", pattern = "*_long.pep.pfamscan.csv", full.names = TRUE)

GENE_FAMILIES_SEARCH_INFO <- "../../../../climate-adaptation-corals-ops/results_annotation/data/gene_families_searchinfo.csv"

################################################################################
# LOAD DATA
################################################################################
########################################
# data from pipeline run samplesheet
samplesheet <- read.delim(SAMPLESHEET, sep = ",")
species_list <- unique(samplesheet$id)
########################################
# Get all defline info 
file_paths <- list.files(paste0(RESULT_DIR, "prefilter/initial/defline_info/"), 
                          pattern = "*.csv", full.names = TRUE)
seqinf <- lapply(file_paths, read.csv)
seqinf <- do.call(rbind, seqinf)

########################################
# Prepare functional annotations for your genes
# Filter down to only csvs for ids in the pipeline run (filename dependent)
LONG_PEP_PFAMSCAN_CSVS <- LONG_PEP_PFAMSCAN_CSVS[sapply(LONG_PEP_PFAMSCAN_CSVS, function(x) {
  any(sapply(species_list, function(species) grepl(species, x)))
})]

# This is the specific header for these annotations
funann_head <- c("pfam_seq_id", "pfam_alignment_start", "pfam_alignment_end", "pfam_envelope_start", "pfam_envelope_end", 
                  "pfam_hmm_acc", "pfam_hmm_name", "pfam_type", "pfam_hmm_start", "pfam_hmm_end", "pfam_hmm_length", 
                  "pfam_bit_score", "pfam_E_value", "pfam_significance", "pfam_clan")

# Read the filtered CSV files
funann <- lapply(LONG_PEP_PFAMSCAN_CSVS, clean_pfam_scan_file, header = funann_head)
# Filter down to only genes within the OG dataset

funann <- lapply(funann, function(annotation_df) {
  annotation_df %>%
    filter(pfam_seq_id %in% seqinf$parent_seq)  # Filter where seq_id is in parent_seq
})

# Final merge of the annotations
funann <- do.call(rbind, funann)

# Join defline info and functional annotations into long format
def_fun <- seqinf %>%
  left_join(funann, by = c("parent_seq" = "pfam_seq_id"))

########################################
# Prep COGEQC inputs
# Reduce to annotations for COGEQC
# clean_parent_seq is used to use whole protein annotations for COGEQC

funann_cogeqc <- def_fun %>%
    filter(!is.na(pfam_hmm_acc)) %>%
    filter(pfam_type == "Domain") %>%
    select(
      id,
      clean_parent_seq,
      pfam_hmm_acc
    ) %>%
    rename(
      Species = id,
      Gene = clean_parent_seq,
      Annotation = pfam_hmm_acc
    ) %>%
    unique()
  
# split into lists by species
funann_cogeqc <- funann_cogeqc %>% 
  split(funann_cogeqc$Species) %>%
  lapply(function(df) {
    df %>%
      select(-Species)
  })


########################################
# Load Orthogroup info for each OG caller
orthgs <- list()
orthgs$orthofinder <- cogeqc::read_orthogroups(OGS_TSV_PATH_OF) %>%
  filter(!is.na(Orthogroup)) %>%
  mutate(OG_source = "orthofinder") %>%
  unique()
orthgs$broccoli <- cogeqc::read_orthogroups(OGS_TSV_PATH_BR) %>%
  filter(!is.na(Orthogroup)) %>%
  mutate(OG_source = "broccoli") %>%
  unique()
orthgs <- do.call(rbind, orthgs)

# reshape for COGEQC
# Here we need to change from clean_seq to clean_parent_seq (use genes not domains)
orthgs_cogeqc <- list()

name_replace_cogeqc <- seqinf %>%
  select(clean_seq, clean_parent_seq) %>%
  unique()


orthgs_cogeqc$orthofinder <- cogeqc::read_orthogroups(OGS_TSV_PATH_OF) %>%
  filter(!is.na(Orthogroup)) %>%
  left_join(name_replace_cogeqc, by = c("Gene" = "clean_seq")) %>%
  select(-Gene) %>%
  rename(Gene = clean_parent_seq) %>%
  unique()
orthgs_cogeqc$broccoli <- cogeqc::read_orthogroups(OGS_TSV_PATH_BR) %>%
  filter(!is.na(Orthogroup)) %>%
  left_join(name_replace_cogeqc, by = c("Gene" = "clean_seq")) %>%
  select(-Gene) %>%
  rename(Gene = clean_parent_seq) %>%
  unique()

# Calculate mean orthogroup homogeneity scores across species
# Some OGs have no scores despite having many seqs e.g. OG_12 broccoli
# It has scores from assess_orthogroups
# but NOT from get_cogeqc_hscores
# It has many functional annotations
# cogeqc should use whole protein as genes, not separate domains
MAX_OG_SIZE = 2000
og_hscores <- list()

og_hscores[["orthofinder"]] <- get_cogeqc_hscores(orthgs_cogeqc[["orthofinder"]], funann_cogeqc, max_og_size = MAX_OG_SIZE)

og_hscores[["broccoli"]] <- get_cogeqc_hscores(orthgs_cogeqc[["broccoli"]], funann_cogeqc, max_og_size = MAX_OG_SIZE)

# Generate  homogeneity scores between tools/runs, make scores scaled
og_hscores_scaled_all <- compare_homogeneity_scores(
  og_hscores[["orthofinder"]],
  og_hscores[["broccoli"]],
  comps = list(c("orthofinder", "broccoli"))
)
# Add COGEQC scores to OG info
orthgs <- orthgs %>%
  left_join(og_hscores_scaled_all, by = c("Orthogroup" = "Orthogroup")) %>%
  select(-Source) %>%
  rename(
      cogeqc_hscore = Score,
      cogeqc_hscore_scaled_self = Score_scaled_self,
      cogeqc_hscore_scaled_all = Score_scaled_all
  ) %>%
  unique()

########################################
# Merge useful information with the defline info

m <- def_fun %>%
    left_join(samplesheet, by = c("id" = "id")) %>%
    rename(
      supergroup = taxonomy
    ) %>%
    select(-fasta)

# Merge with orthogroups data
def_fun_sam_ogs <- m %>%
    left_join(orthgs, by = c("clean_seq" = "Gene")) %>%
    select(-Species)

# sum total of certain vars
total_seqs_all <- n_distinct(def_fun_sam_ogs$seq)
total_parent_seqs_all <- n_distinct(def_fun_sam_ogs$parent_seq)
total_ids_all <- n_distinct(def_fun_sam_ogs$id)
total_supergroups_all <- n_distinct(def_fun_sam_ogs$supergroup)

def_fun_sam_ogs <- def_fun_sam_ogs %>%
  group_by(Orthogroup) %>%
  mutate(
    total_seqs = length(unique(seq)),
    total_parent_seqs = length(unique(parent_seq)),
    total_ids = length(unique(id)),
    total_supergroups = length(unique(supergroup)),
    
    most_common_pfam_hmm_name = get_mode(pfam_hmm_name),
    total_most_common_pfam_hmm_name = count_mode(pfam_hmm_name, get_mode(pfam_hmm_name)),
    percent_most_common_pfam_hmm_name = (total_most_common_pfam_hmm_name / length(pfam_hmm_name)) * 100,

    # Percentage of the total across the entire dataset
    percent_total_seqs = total_seqs / total_seqs_all * 100,
    percent_total_parent_seqs = total_parent_seqs / total_parent_seqs_all * 100,
    percent_total_ids = total_ids / total_ids_all * 100,
    percent_total_supergroups = total_supergroups / total_supergroups_all * 100,
    
    # Categorize total_seqs with exclusive ranges
    og_size_category_seqs = case_when(
      total_seqs > 10000 ~ ">10000",
      total_seqs > 1000  ~ "1000-10000",
      total_seqs > 500   ~ "500-1,000",
      total_seqs > 200   ~ "200-500",
      total_seqs > 100   ~ "100-200",
      total_seqs >= 50   ~ "50-100",
      total_seqs < 50    ~ "10-50",
      TRUE ~ "<10"
    ),
    # Categorize total_parent_seqs with exclusive ranges
    og_size_category_parent_seqs = case_when(
      total_parent_seqs > 10000 ~ ">10000",
      total_parent_seqs > 1000  ~ "1000-10000",
      total_parent_seqs > 500   ~ "500-1,000",
      total_parent_seqs > 200   ~ "200-500",
      total_parent_seqs > 100   ~ "100-200",
      total_parent_seqs >= 50   ~ "50-100",
      total_parent_seqs < 50    ~ "10-50",
      TRUE ~ "<10"
    )
  ) %>%
  ungroup()


# jaccard similarity between all OGs

# Create lists of sequences associated with each algorithm and OG
jaccard_in_seq <- def_fun_sam_ogs %>%
  select(clean_parent_seq, Orthogroup, OG_source) %>%
  unique() %>%
  split(.$OG_source) %>%
  lapply(function(x) {
    split(x, x$Orthogroup) %>%
      lapply(function(group) {
        group %>%
          select(clean_parent_seq) %>%
          pull(clean_parent_seq)  # Pull the sequences as a vector
      })
  })


# Create an empty list to store the results
jaccard_mat <- list()
jaccard_metrics <- list()
jaccard_metrics_per_og <- list()

# Get all unique pairs of OG_sources
og_sources <- names(jaccard_in_seq)
og_sources <- expand.grid(og_sources, og_sources)

# Iterate over all pairs and calculate Jaccard similarity

for (i in 1:nrow(og_sources)) { 
  og_source1 <- og_sources[i, 1]
  og_source2 <- og_sources[i, 2]
  pair_name <- paste(og_source1, "vs", og_source2, sep = "_")
  
  jaccard_mat[[pair_name]] <- calculate_jaccard_similarity(jaccard_in_seq[[og_source1]], jaccard_in_seq[[og_source2]])

  jaccard_metrics[[pair_name]] <- calculate_jaccard_metrics(jaccard_mat[[pair_name]], og_source1, og_source2)
}

jaccard_metrics_all <- bind_rows(jaccard_metrics)

# Overall summary for tools
jaccard_metrics_all_long <- jaccard_metrics_all %>%
  gather(key = "Metric", value = "Value", -tool1, -tool2)

# Now for OG specific jaccard info
for (i in 1:nrow(og_sources)) { 
  og_source1 <- og_sources[i, 1]
  og_source2 <- og_sources[i, 2]
  pair_name <- paste(og_source2, "vs", og_source1, sep = "_")
  
  jaccard_metrics_per_og[[pair_name]] <- calculate_jaccard_metrics_per_og(jaccard_mat[[pair_name]], og_source2, og_source1)
}

jaccard_metrics_per_og_all <- bind_rows(jaccard_metrics_per_og) %>% 
  unique() %>%
  filter(OG_source != vs_OG_source)
rownames(jaccard_metrics_per_og_all) <- NULL

# Merge jaccard info with rest of data

def_fun_sam_ogs_jac <- def_fun_sam_ogs %>%
    left_join(jaccard_metrics_per_og_all, by = c("Orthogroup" = "Orthogroup", "OG_source" = "OG_source"))

main_seq_long = def_fun_sam_ogs_jac
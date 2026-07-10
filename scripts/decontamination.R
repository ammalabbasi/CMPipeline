#!/usr/bin/env Rscript

## ==============================================================================
## SECTION 1: Load Libraries and Define Helper Functions
## ==============================================================================


# Load necessary libraries
if (!require("pacman")) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse, phyloseq, decontam, argparse, tools)

    Loading required package: pacman
    


## --- Helper Function 1: Data Loading and Preprocessing ---
# This function encapsulates all initial data handling steps to ensure consistency.
load_and_prepare_data <- function(otu_path, meta_path) {
    # Log the action
  cat("INFO: Loading and preparing data...\n")
  
  # Read the raw OTU table and metadata from the provided paths
  otu_data_raw <- read.csv(otu_path, sep="\t", check.names = FALSE)
  metadata_raw <- read.csv(meta_path, sep="\t", check.names = FALSE)
  
  # --- OTU table formatting ---
  
  # Dynamically find the column to be used for row names (taxa names)
  name_col <- intersect(c("name", "species_name", "X", "species"), colnames(otu_data_raw))[1]
  if (is.na(name_col)) {
    stop("ERROR: A suitable column for taxon names was not found in the OTU table.")
  }
  rownames(otu_data_raw) <- otu_data_raw[[name_col]]
  
  # Isolate numeric columns and convert NAs to 0
  numeric_cols <- sapply(otu_data_raw, is.numeric)
  otu_data <- as.data.frame(otu_data_raw[, numeric_cols])
  otu_data[is.na(otu_data)] <- 0
  
  # --- Data Pruning and Quality Control ---

  # Prune empty samples (equivalent to prune_samples(sample_sums > 0))
  original_sample_count <- ncol(otu_data)
  otu_data <- otu_data[, (colSums(otu_data, na.rm=TRUE) != 0)]
  cat("INFO: Pruned", original_sample_count - ncol(otu_data), "empty samples.\n")

  # Prune very rare taxa (singletons), equivalent to prune_taxa(taxa_sums > 1)
  original_taxa_count <- nrow(otu_data)
  taxa_sums <- rowSums(otu_data)
  taxa_to_keep <- taxa_sums > 1
  otu_data <- otu_data[taxa_to_keep, ]
  cat("INFO: Pruned", original_taxa_count - nrow(otu_data), "rare taxa (singletons).\n")

  # --- Match OTU table and metadata ---
  
  # Find and filter for common samples
  common_samples <- intersect(colnames(otu_data), metadata_raw$sampleid)
  if (length(common_samples) == 0) {
    stop("ERROR: No common samples found. Please check sample IDs in both files.")
  }
  
  otu_data <- otu_data[, common_samples]
  metadata <- metadata_raw[metadata_raw$sampleid %in% common_samples, ]
  rownames(metadata) <- metadata$sampleid
  
  cat("INFO:", length(common_samples), "common samples prepared for final analysis.\n")
  
  # Return the processed data as a list
  return(list(otu = otu_data, meta = metadata, name_col = name_col))
}


## --- Helper Function 2: Manual Contaminant Identification ---
identify_contaminants_manual <- function(otu_df, meta_df) {
  cat("INFO: Identifying contaminants using Manual method...\n")
  
  microbiome_data_long <- otu_df %>%
    rownames_to_column(var = "name") %>%
    pivot_longer(cols = -name, names_to = "sample_id", values_to = "reads") %>%
    filter(reads > 0)
    
  microbiome_data_meta <- inner_join(microbiome_data_long,
                                     dplyr::select(meta_df, sampleid, shipment_batch, Type),
                                     by = c("sample_id" = "sampleid"))

  master_table <- data.frame(Contaminant = character(0), Batch_Count = character(0))

  for(batch_id in unique(microbiome_data_meta$shipment_batch)){
    microbiome_data_batch <- microbiome_data_meta %>% dplyr::filter(shipment_batch == batch_id)
    
    # Apply standardized validation check
    if(nrow(microbiome_data_batch) <= 2 || n_distinct(microbiome_data_batch$Type) < 2) {
        cat("WARN: [Manual] Skipping batch", batch_id, "- not enough samples or types for comparison.\n")
        next
    }
    
    microbiome_data_batch <- microbiome_data_batch %>% mutate(cohort = ifelse(tolower(Type) == "tumor", "tumor", "blood"))
    
    total_tumor_samples <- n_distinct((microbiome_data_batch %>% filter(cohort == "tumor"))$sample_id)
    total_blood_samples <- n_distinct((microbiome_data_batch %>% filter(cohort == "blood"))$sample_id)
    
    sample_detection_counts <- microbiome_data_batch %>%
      group_by(name, cohort) %>%
      summarise(n_detected = n_distinct(sample_id), .groups = 'drop') %>%
      pivot_wider(names_from = cohort, values_from = n_detected, values_fill = 0)

    if (!"tumor" %in% names(sample_detection_counts)) sample_detection_counts$tumor <- 0
    if (!"blood" %in% names(sample_detection_counts)) sample_detection_counts$blood <- 0
      
    species_stats <- sample_detection_counts %>%
      mutate(prevalence_tumor = tumor / total_tumor_samples, prevalence_blood = blood / total_blood_samples) %>%
      
      # Apply the minimum prevalence filter to match the original script's logic
      filter(tumor >= total_tumor_samples * 0.02 & blood >= total_blood_samples * 0.02) %>%
      
      rowwise() %>%
      mutate(p_value_fisher = fisher.test(matrix(c(tumor, blood, total_tumor_samples - tumor, total_blood_samples - blood), nrow = 2))$p.value) %>%
      ungroup() %>%
      mutate(q_value_fisher = p.adjust(p_value_fisher, method = "BH"),
             prevalence_FC = ifelse(prevalence_blood == 0, Inf, prevalence_tumor / prevalence_blood))

    contaminants_in_batch <- species_stats %>%
      filter(!(prevalence_FC > 1 & q_value_fisher < 0.05)) %>% # whitelist approach 
      # filter(prevalence_FC < 1 & q_value_fisher < 0.1) %>%  # blacklist approach 
      pull(name)
      
    if(length(contaminants_in_batch) > 0) {
        master_table <- rbind(master_table, data.frame(Contaminant = contaminants_in_batch, Batch_Count = batch_id))
    }
  }
  return(master_table)
}

## --- Helper Function 3: Decontam Package Identification ---
identify_contaminants_decontam <- function(otu_df, meta_df) {
  cat("INFO: Identifying contaminants using decontam package...\n")
  master_table <- data.frame(Contaminant = character(0), Batch_Count = character(0))

  for(batch_id in unique(meta_df$shipment_batch)) {
    meta_batch <- meta_df[meta_df$shipment_batch == batch_id, ]
    
    # Apply standardized validation check
    if(nrow(meta_batch) <= 2 || n_distinct(meta_batch$Type) < 2) {
        cat("WARN: [decontam] Skipping batch", batch_id, "- not enough samples or types for decontam.\n")
        next
    }
    
    otu_batch <- otu_df[, meta_batch$sampleid]
    
    OTU <- otu_table(as.matrix(otu_batch), taxa_are_rows = TRUE)
    meta_batch$is.neg <- meta_batch$Type == "Normal"
    META <- sample_data(meta_batch)
    
    if(sum(sample_data(META)$is.neg, na.rm=TRUE) == 0){
        cat("WARN: [decontam] Skipping batch", batch_id, "- no negative controls ('Normal' type) found.\n")
        next
    }
    
    physeq_batch <- phyloseq(OTU, META)
    contam_results <- isContaminant(physeq_batch, method="prevalence", neg="is.neg", threshold=0.1)
    
    contaminants_in_batch <- rownames(contam_results)[contam_results$contaminant == TRUE]
    
    if(length(contaminants_in_batch) > 0){
        master_table <- rbind(master_table, data.frame(Contaminant = contaminants_in_batch, Batch_Count = batch_id))
    }
  }
  return(master_table)
}

## --- Helper Function 4: Common Contaminant Removal ---
# This function takes an OTU table and a list of contaminants and returns a clean table.
remove_contaminants <- function(original_otu, master_contaminant_table, meta_df) {
  clean_otu <- original_otu
  
  if (nrow(master_contaminant_table) > 0) {
    cat("INFO: Removing", nrow(master_contaminant_table), "contaminant entries across all batches.\n")
    for (i in 1:nrow(master_contaminant_table)) {
      contaminant_name <- master_contaminant_table$Contaminant[i]
      batch_val <- master_contaminant_table$Batch_Count[i]
      
      # Find samples belonging to the specific batch
      samples_in_batch <- meta_df$sampleid[meta_df$shipment_batch == batch_val]
      
      # Set read counts to 0 for the contaminant in the specified samples
      if (contaminant_name %in% rownames(clean_otu)) {
          clean_otu[contaminant_name, samples_in_batch] <- 0
      }
    }
  } else {
    cat("INFO: No contaminants to remove.\n")
  }
  return(clean_otu)
}

## --- Helper Function 5: Save Decontaminated OTU Table ---
# This function provides a standardized way to write the final output files.
save_decontaminated_table <- function(clean_otu_table, original_name_col, file_prefix, method_suffix) {
  
  # Construct the full output file path
  output_file_path <- paste0(file_prefix, ".", method_suffix, "_decontaminated.csv")
  
  # Log the saving action
  cat("INFO: Writing", method_suffix, "decontaminated OTU table to:", output_file_path, "\n")
  
  # Convert the row names (taxa) back to the first column with its original name
  # and write the data frame to a CSV file.
  clean_otu_table %>%
    rownames_to_column(var = original_name_col) %>%
    write.csv(., output_file_path, row.names = FALSE)
}

## ==============================================================================
## SECTION 2: Main Workflow Execution
## ==============================================================================


if(T){
    force_interactive <- TRUE
}


# --- Argument Parsing ---
# Define and parse command-line arguments

if(interactive() || force_interactive) {
    args <- list(
        otu_table = "/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/02.Results/250210_CRC_CMPipeline/OTU_table_signature_metrics_after_consensus_taxa_Ammal/species/crc.bracken.raw.common.species.report.txt",
        metadata = "/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/99.Raw_data/CRC/00.Metadata_TMB_Signatures_SBS_ID_DBS_CN_SV_Drivers_2024AUG21_v2.tsv",
        prefix = "02.Results/test_run"
        )
} else{
    parser <- ArgumentParser(description="Perform manual and decontam-based decontamination on an OTU table.")
    parser$add_argument("--otu_table", type="character", required=TRUE, help="Path to the input OTU table file.")
    parser$add_argument("--metadata", type="character", required=TRUE, help="Path to the metadata file.")
    parser$add_argument("--prefix", type="character", required=TRUE, help="Prefix for the output files.")
    args <- parser$parse_args()
}


# --- Initial Data Preparation ---
# Use the helper function to load and preprocess data one time
prepared_data <- load_and_prepare_data(args$otu_table, args$metadata)
otu_original <- prepared_data$otu
meta_original <- prepared_data$meta
original_taxon_col_name <- prepared_data$name_col

# --- Run Decontamination Method 1: Manual Prevalence-Based ---
contaminants_manual <- identify_contaminants_manual(otu_original, meta_original)
otu_clean_manual <- remove_contaminants(otu_original, contaminants_manual, meta_original)
save_decontaminated_table(otu_clean_manual, original_taxon_col_name, args$prefix, "manual")

# --- Run Decontamination Method 2: Using 'decontam' Package ---
contaminants_decontam <- identify_contaminants_decontam(otu_original, meta_original)
otu_clean_decontam <- remove_contaminants(otu_original, contaminants_decontam, meta_original)
save_decontaminated_table(otu_clean_decontam, original_taxon_col_name, args$prefix, "decontam_pkg")

# --- End of Script ---
cat("INFO: Decontamination script finished successfully.\n")

#!/usr/bin/env Rscript
# ==============================================================================
#  BATCH-CORRECTION & NORMALISATION PIPELINE - Auto-Phased Tune_ConQuR
#  
#  Key improvements:
#  - Phase 1: Fast test (16 combinations, ~10-20 min)
#  - Phase 3: Comprehensive search (72 combinations, ~60-90 min)
#  - PERMANOVA-based stopping criteria (25% R² reduction)
#  - Automatic phase selection based on results






# Example usage (basic):
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_only

# Example usage (custom columns):
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --batch_column "batch_id" \
#   --covariates "age,gender,BMI" \
#   --type_column "sample_type" \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor" \
#   --tumor_only

# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE))
      install.packages("pacman", repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools", repos = "https://cloud.r-project.org")

  # Install ConQuR from GitHub if not available
  # Use ivartb/ConQuR_par fork which fixes the "batchid not found" error in foreach workers
  if (!requireNamespace("ConQuR", quietly = TRUE)) {
    message("INFO: Installing ConQuR (ivartb/ConQuR_par fork with batchid fix) from GitHub...")
    tryCatch({
      devtools::install_github("ivartb/ConQuR_par", quiet = TRUE, upgrade = "never")
    }, error = function(e) {
      message("INFO: Falling back to original ConQuR...")
      devtools::install_github("wdl2459/ConQuR", quiet = TRUE, upgrade = "never")
    })
  }

  pacman::p_load(tidyverse, vegan, ConQuR, argparse,
               zCompositions, devtools,
               doParallel, foreach, iterators) 

  ## DEICODE (for RCLR)
  if (!requireNamespace("DEICODE", quietly = TRUE)) {
    try(BiocManager::install("DEICODE", update = FALSE, ask = FALSE), silent = TRUE)
    if (!requireNamespace("DEICODE", quietly = TRUE))
      try(devtools::install_github("DEICODE-dev/DEICODE", quiet = TRUE), silent = TRUE)
  }
  if (requireNamespace("DEICODE", quietly = TRUE)) {
    library(DEICODE);  rclr_f <- DEICODE::rclr
  } else {
    message("WARN: DEICODE unavailable — using manual RCLR")
    rclr_f <- function(mat){
      t(apply(mat, 1, function(x){
        nz <- x>0; g <- exp(mean(log(x[nz]))); y <- rep(NA, length(x))
        y[nz] <- log(x[nz]/g); y
      }))
    }
  }
})

# ---------- Utility functions -------------------------------------------------
save_table <- function(df, path){
  if (!is.data.frame(df))             
      df <- as.data.frame(df)
  
  # Explicitly preserve column names
  col_names <- colnames(df)
  
  df_out <- df %>%
    rownames_to_column("sample_id")
  
  # Force column names to be written
  write_tsv(df_out, path, progress = FALSE, col_names = TRUE)
  
  cat("  → Saved", nrow(df_out), "rows x", ncol(df_out)-1, "taxa to", basename(path), "\n")
}

load_table <- function(path){
  # Auto-detect CSV vs TSV
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
  
  # Identify the first column (taxon names or sample IDs)
  first_col <- colnames(df)[1]
  
  # Check if this is OTU table format (taxa as rows) or already transposed
  if (first_col %in% c("species_name", "name", "taxa", "taxon", "genus", "species", "OTU", "ASV")) {
    # This is OTU table format (taxa as rows) - need to transpose
    cat("  → Detected OTU table format, transposing...\n")
    
    df_out <- df %>%
      column_to_rownames(first_col) %>%
      t() %>%
      as.data.frame()
    
  } else if (first_col == "sample_id") {
    # Already transposed format (samples as rows)
    df_out <- df %>%
      column_to_rownames("sample_id") %>%
      as.data.frame()
    
  } else {
    # Assume first column contains sample IDs
    cat("  → Assuming first column contains sample IDs\n")
    df_out <- df %>%
      column_to_rownames(first_col) %>%
      as.data.frame()
  }
  
  cat("  → Loaded", nrow(df_out), "samples x", ncol(df_out), "taxa from", basename(path), "\n")
  
  return(df_out)
}

# ---------- Load and prepare data ---------------------------------------------
load_and_prepare_data <- function(otu_path, meta_path) {
  cat("INFO: Loading and preparing data for batch correction...\n")
  
  otu_raw <- read_csv(otu_path, show_col_types = FALSE)

  # Taxon name column 자동 감지
  taxon_name_candidates <- c("species_name", "name", "X", "species", "taxon",
                             "taxa", "clade_name", "Taxon", "Species",
                             "organism", "OTU", "ASV", "#OTU ID", "OTUID")
  name_col <- taxon_name_candidates[taxon_name_candidates %in% colnames(otu_raw)][1]

  # 후보에서 찾지 못했으면 첫 번째 non-numeric 컬럼 사용
  if (is.na(name_col)) {
    cat("WARN: Standard taxon name column not found. Available OTU table columns:\n")
    cat("  ", paste(head(colnames(otu_raw), 10), collapse = ", "), "...\n")

    # 첫 번째 문자열 컬럼 찾기
    non_numeric_cols <- which(!sapply(otu_raw, is.numeric))
    if (length(non_numeric_cols) > 0) {
      name_col <- colnames(otu_raw)[non_numeric_cols[1]]
      cat("INFO: Using first non-numeric column '", name_col, "' as taxon name.\n", sep = "")
    } else {
      stop("ERROR: Cannot identify taxon name column in OTU table.\n",
           "Available columns: ", paste(colnames(otu_raw), collapse = ", "))
    }
  } else {
    cat("INFO: Using '", name_col, "' as taxon name column.\n", sep = "")
  }

  otu_data <- otu_raw %>%
    column_to_rownames(name_col) %>%
    t() %>%
    as.data.frame()
  
  original_sample_count <- nrow(otu_data)
  otu_data <- otu_data[rowSums(otu_data) > 0, ]
  cat("INFO: Removed", original_sample_count - nrow(otu_data), "empty samples.\n")

  # Sample ID normalization - remove common suffixes from OTU sample IDs
  cat("\nINFO: Normalizing sample IDs...\n")
  original_sample_ids <- rownames(otu_data)
  normalized_sample_ids <- gsub("\\.bracken\\.[GS]\\.krakenreport\\.txt$", "", original_sample_ids)
  normalized_sample_ids <- gsub("\\.bracken\\.[GS]\\.mpa\\.krakenreport\\.txt$", "", normalized_sample_ids)
  cat("INFO: Example OTU sample IDs (before):", paste(head(original_sample_ids, 3), collapse = ", "), "\n")
  cat("INFO: Example OTU sample IDs (after):", paste(head(normalized_sample_ids, 3), collapse = ", "), "\n")
  rownames(otu_data) <- normalized_sample_ids

  meta_raw <- read_tsv(meta_path, show_col_types = FALSE)

  # 빈 컬럼명 처리 (Empty column names)
  col_names <- colnames(meta_raw)
  empty_cols <- which(col_names == "" | is.na(col_names))
  if(length(empty_cols) > 0) {
    cat("WARN: Found", length(empty_cols), "empty column name(s). Removing...\n")
    meta_raw <- meta_raw[, -empty_cols, drop = FALSE]
  }

  # Sample ID column 자동 감지
  sample_id_candidates <- c("donor_id", "sample_id", "sampleid", "Sample.ID", "Sample_ID",
                            "SampleID", "Sample ID", "SAMPLE_ID", "Sample ID_x",
                            "subject_id", "SubjectID", "Subject_ID", "patient", "Patient")
  sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(meta_raw)][1]

  if(is.na(sample_id_col)) {
    stop("ERROR: Sample ID column not found in metadata. Tried: ",
         paste(sample_id_candidates, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(meta_raw), collapse = ", "))
  }

  cat("INFO: Using '", sample_id_col, "' as sample ID column.\n", sep = "")

  # sample_id로 컬럼명 표준화
  if(sample_id_col != "sample_id") {
    meta_raw <- meta_raw %>% rename(sample_id = !!sym(sample_id_col))
  }

  # Sample Type 처리: 여러 가능한 컬럼명 지원
  type_candidates <- c("Type", "type", "Sample Type_x", "Sample_Type", "sample_type",
                       "cohort", "Cohort", "group", "Group", "condition", "Condition")
  type_col <- type_candidates[type_candidates %in% colnames(meta_raw)][1]

  if(!is.na(type_col) && type_col != "Type") {
    cat("INFO: Using '", type_col, "' column for Type information\n", sep = "")
    meta_raw$Type <- meta_raw[[type_col]]
  } else if(is.na(type_col)) {
    cat("WARN: No Type column found in metadata. Creating default Type='Sample'\n")
    meta_raw$Type <- "Sample"
  }

  if (any(duplicated(meta_raw$sample_id))) {
    dup_ids <- unique(meta_raw$sample_id[duplicated(meta_raw$sample_id)])
    cat("WARN:", length(dup_ids), "duplicate sample_ids. Keeping first occurrence.\n")
    meta_raw <- meta_raw[!duplicated(meta_raw$sample_id), ]
  }

  common_samples <- intersect(rownames(otu_data), meta_raw$sample_id)

  if (length(common_samples) == 0) {
    stop("No common samples between OTU table and metadata")
  }

  otu_final <- otu_data[common_samples, , drop = FALSE]
  meta_final <- meta_raw %>%
    filter(sample_id %in% common_samples) %>%
    column_to_rownames("sample_id")
  meta_final <- meta_final[common_samples, , drop = FALSE]
  
  cat("INFO:", length(common_samples), "common samples prepared.\n")
  
  return(list(otu = otu_final, meta = meta_final))
}

# ---------- Select optimal reference batch ------------------------------------
select_optimal_reference_batch <- function(tab, meta, batch_col) {
  
  cat("\n=== SELECTING OPTIMAL REFERENCE BATCH ===\n")
  
  batches <- unique(meta[[batch_col]])
  batch_distances <- matrix(NA, length(batches), length(batches))
  rownames(batch_distances) <- batches
  colnames(batch_distances) <- batches
  
  batch_centroids <- list()
  for (batch in batches) {
    batch_samples <- rownames(meta[meta[[batch_col]] == batch, ])
    batch_data <- tab[batch_samples, ]
    batch_centroids[[batch]] <- colMeans(batch_data)
  }
  
  for (i in 1:length(batches)) {
    for (j in 1:length(batches)) {
      if (i != j) {
        dist_matrix <- vegdist(rbind(
          batch_centroids[[batches[i]]], 
          batch_centroids[[batches[j]]]
        ), method = "bray")
        batch_distances[i, j] <- as.numeric(dist_matrix)
      } else {
        batch_distances[i, j] <- 0
      }
    }
  }
  
  mean_distances <- rowMeans(batch_distances, na.rm = TRUE)
  
  result_df <- data.frame(
    Batch = names(mean_distances),
    Mean_Distance = mean_distances,
    N_Samples = sapply(batches, function(b) sum(meta[[batch_col]] == b))
  )
  result_df <- result_df[order(result_df$Mean_Distance), ]
  
  cat("\nBatch centrality ranking:\n")
  print(result_df, row.names = FALSE)
  
  top_batches <- result_df$Batch[1:min(3, nrow(result_df))]
  
  cat("\nTop", length(top_batches), "most central batches selected\n")
  
  return(as.character(top_batches))
}

# ---------- PERMANOVA ---------------------------------------------------------
run_permanova <- function(mat, meta, rhs_formula, dist_meth){
  smp  <- intersect(rownames(mat), rownames(meta))
  if (anyNA(mat[smp, ])) mat[smp, ][is.na(mat[smp, ])] <- 0
  mat <- mat[ , colSums(mat) > 0, drop = FALSE]
  dist <- vegdist(mat[smp, ], dist_meth)
  res  <- adonis2(dist ~ ., data = meta[smp, all.vars(rhs_formula)],
                  permutations = 999, by = "margin")
  as_tibble(res, rownames = "Term") %>% dplyr::select(Term, Df, R2, `Pr(>F)`)
}

# ---------- Evaluate batch effect (PERMANOVA-based) ---------------------------
evaluate_batch_effect <- function(data, meta, batch_col, covariates) {
  
  formula <- reformulate(c(batch_col, covariates))
  perm_result <- run_permanova(data, meta, formula, "bray")
  
  batch_r2 <- perm_result %>% 
    filter(Term == batch_col) %>% 
    pull(R2)
  
  if (length(batch_r2) == 0) {
    warning("Batch term not found in PERMANOVA results")
    return(list(
      batch_r2 = NA,
      batch_p = NA,
      covariate_r2 = NA,
      full_result = perm_result
    ))
  }
  
  batch_p <- perm_result %>% 
    filter(Term == batch_col) %>% 
    pull(`Pr(>F)`)
  
  # Covariate R² (for overcorrection check)
  covariate_r2 <- perm_result %>%
    filter(Term %in% covariates) %>%
    pull(R2) %>%
    sum()
  
  return(list(
    batch_r2 = batch_r2,
    batch_p = batch_p,
    covariate_r2 = covariate_r2,
    full_result = perm_result
  ))
}

# ---------- Run single phase --------------------------------------------------
run_phase <- function(phase_num, tab, meta, batch_col, covars, 
                      ref_batch_pool, final_covars) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("PHASE", phase_num, "EXECUTION\n")
  cat(strrep("=", 70), "\n")
  
  # Define parameter pools
  if (phase_num == 1) {
    # PHASE 1: 16 COMBINATIONS (~10-20 minutes)
    batch_pool <- ref_batch_pool[1:min(2, length(ref_batch_pool))]
    logistic_pool <- c(FALSE, TRUE)
    quantile_pool <- c("standard", "lasso")
    lambda_pool <- c(NA, "2p/logn")
    interplt_pool <- c(FALSE)
    
    cat("Strategy: Fast essential coverage\n")
    
  } else {
    # PHASE 3: 72 COMBINATIONS (~60-90 minutes)
    batch_pool <- ref_batch_pool
    logistic_pool <- c(FALSE, TRUE)
    quantile_pool <- c("standard", "lasso")
    lambda_pool <- c(NA, "2p/n", "2p/logn")
    interplt_pool <- c(FALSE, TRUE)
    
    cat("Strategy: Comprehensive grid search\n")
  }
  
  n_combinations <- length(batch_pool) * length(logistic_pool) * 
                    length(quantile_pool) * length(lambda_pool) * 
                    length(interplt_pool)
  
  cat("\nParameter space:\n")
  cat("  Reference batches:", length(batch_pool), "\n")
  cat("  Logistic lasso:   ", length(logistic_pool), "(FALSE, TRUE)\n")
  cat("  Quantile type:    ", length(quantile_pool), "(standard, lasso)\n")
  cat("  Lambda:           ", length(lambda_pool), "\n")
  cat("  Interpolation:    ", length(interplt_pool), "\n")
  cat("\nTotal combinations:", n_combinations, "\n")
  
  if (phase_num == 1) {
    cat("Estimated time:    10-20 minutes\n")
  } else {
    cat("Estimated time:    60-90 minutes\n")
  }
  cat(strrep("=", 70), "\n\n")
  
  start_time <- Sys.time()

  # Run Tune_ConQuR - Keep it simple like the working original
  options(warn = -1)
  registerDoSEQ()

  # CRITICAL: Assign batchid to global environment for foreach workers
  # ConQuR uses foreach internally and workers need access to batchid
  batchid <- factor(meta[[batch_col]])
  assign("batchid", batchid, envir = .GlobalEnv)

  tryCatch({
    result <- Tune_ConQuR(
      tax_tab = tab,
      batchid = batchid,
      covariates = final_covars,
      batch_ref_pool = batch_pool,
      logistic_lasso_pool = logistic_pool,
      quantile_type_pool = quantile_pool,
      simple_match_pool = FALSE,
      lambda_quantile_pool = lambda_pool,
      interplt_pool = interplt_pool,
      frequencyL = 0,
      frequencyU = 1
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)
    
    cat("\n✓ Phase", phase_num, "completed in", round(elapsed, 1), "minutes\n")
    
    return(list(
      data = as.data.frame(result$tax_final),
      method_matrix = result$method_final,
      elapsed = elapsed,
      phase = phase_num,
      success = TRUE
    ))
    
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)

    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! PHASE", phase_num, "FAILED after", round(elapsed, 1), "minutes !!!\n")
    cat(strrep("!", 70), "\n")
    cat("Error message:\n")
    cat("  ", e$message, "\n")
    cat(strrep("!", 70), "\n\n")

    return(list(
      data = NULL,
      method_matrix = NULL,
      elapsed = elapsed,
      phase = phase_num,
      success = FALSE,
      error = e$message
    ))
  })
}

# ---------- Vanilla ConQuR (default, non-tuned) --------------------------------
run_vanilla_conqur <- function(tab_raw, meta, batch_col, covars, corr_dir) {

  cat("\n", strrep("=", 70), "\n")
  cat("VANILLA ConQuR MODE (non-tuned)\n")
  cat(strrep("=", 70), "\n\n")

  # Data preparation (same as auto-phase)
  smp <- intersect(rownames(tab_raw), rownames(meta))
  tab <- tab_raw[smp, ]
  meta <- meta[smp, ]

  # Remove unknown batches
  unknown_count <- sum(meta[[batch_col]] == "Unknown" | is.na(meta[[batch_col]]))
  if (unknown_count > 0) {
    cat("WARN: Removing", unknown_count, "samples with Unknown/NA batch\n")
    valid_batch <- meta[[batch_col]] != "Unknown" & !is.na(meta[[batch_col]])
    meta <- meta[valid_batch, ]
    tab <- tab[rownames(meta), ]
  }

  # Remove samples with NA covariates
  complete_rows <- complete.cases(meta[, c(batch_col, covars)])
  n_incomplete <- sum(!complete_rows)
  if (n_incomplete > 0) {
    cat("WARN: Removing", n_incomplete, "samples with NA in covariates\n")
    meta <- meta[complete_rows, ]
    tab <- tab[rownames(meta), ]
  }

  # Batch validation
  batch_counts <- table(meta[[batch_col]])
  small_batches <- names(batch_counts)[batch_counts < 5]
  if (length(small_batches) > 0) {
    cat("WARN: Removing batches with <5 samples:", paste(small_batches, collapse=", "), "\n")
    keep_batches <- !(meta[[batch_col]] %in% small_batches)
    meta <- meta[keep_batches, ]
    tab <- tab[rownames(meta), ]
  }

  # Covariate validation
  keep_covars <- rep(TRUE, length(covars))
  names(keep_covars) <- covars
  batches <- unique(meta[[batch_col]])
  for (i in seq_along(covars)) {
    cov <- covars[i]
    if (!cov %in% colnames(meta)) { keep_covars[i] <- FALSE; next }
    for (batch in batches) {
      batch_data <- meta[meta[[batch_col]] == batch, cov]
      if (is.numeric(batch_data)) {
        v <- var(batch_data, na.rm = TRUE)
        if (is.na(v) || v == 0) { keep_covars[i] <- FALSE; break }
      } else {
        if (length(unique(batch_data[!is.na(batch_data)])) < 2) { keep_covars[i] <- FALSE; break }
      }
    }
  }
  valid_covars <- covars[keep_covars]
  if (length(valid_covars) == 0) {
    final_covars <- NULL
  } else {
    final_covars <- meta[, valid_covars, drop = FALSE]
  }

  meta <- droplevels(meta)

  # Baseline evaluation
  baseline <- evaluate_batch_effect(tab, meta, batch_col, valid_covars)
  cat(sprintf("Baseline Batch R2: %.4f (p = %.4f)\n", baseline$batch_r2, baseline$batch_p))

  # Select reference batch
  ref_batch_pool <- select_optimal_reference_batch(tab, meta, batch_col)
  batch_ref <- ref_batch_pool[1]
  cat("Using reference batch:", batch_ref, "\n\n")

  # Run vanilla ConQuR
  batchid <- factor(meta[[batch_col]])
  assign("batchid", batchid, envir = .GlobalEnv)

  start_time <- Sys.time()
  options(warn = -1)
  registerDoSEQ()

  tryCatch({
    result <- ConQuR(
      tax_tab = tab,
      batchid = batchid,
      covariates = final_covars,
      batch_ref = batch_ref,
      logistic_lasso = FALSE,
      quantile_type = "lasso",
      interplt = FALSE
    )

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)

    corrected <- as.data.frame(result)

    # Evaluate
    post_eval <- evaluate_batch_effect(corrected, meta, batch_col, valid_covars)
    delta_r2 <- baseline$batch_r2 - post_eval$batch_r2
    r2_reduction <- delta_r2 / baseline$batch_r2

    cat(sprintf("\nVanilla ConQuR completed in %.1f minutes\n", elapsed))
    cat(sprintf("Before: Batch R2 = %.4f\n", baseline$batch_r2))
    cat(sprintf("After:  Batch R2 = %.4f\n", post_eval$batch_r2))
    cat(sprintf("Reduction: %.1f%%\n", r2_reduction * 100))

    return(list(
      data = corrected,
      method_matrix = NULL,
      phase_used = "vanilla",
      baseline_r2 = baseline$batch_r2,
      final_r2 = post_eval$batch_r2,
      r2_reduction = r2_reduction,
      elapsed = elapsed,
      success = TRUE
    ))

  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    options(warn = 0)
    cat("\nVanilla ConQuR FAILED after", round(elapsed, 1), "minutes\n")
    cat("Error:", e$message, "\n")
    return(list(data = tab, phase_used = "vanilla", success = FALSE, message = e$message))
  })
}

# ---------- Auto-phase workflow -----------------------------------------------
run_auto_phase <- function(tab_raw, meta, batch_col, covars,
                            corr_dir, r2_threshold = 0.25) {

  cat("\n", strrep("=", 70), "\n")
  cat("AUTO-PHASE MODE: PERMANOVA-BASED STRATEGY\n")
  cat(strrep("=", 70), "\n")
  cat("Decision rule:\n")
  cat("  IF R² reduction ≥", r2_threshold * 100, "% → STOP at Phase 1\n", sep="")
  cat("  ELSE → Run Phase 3 and compare results\n")
  cat(strrep("=", 70), "\n\n")

  # Data preparation
  smp <- intersect(rownames(tab_raw), rownames(meta))
  tab <- tab_raw[smp, ]
  meta <- meta[smp, ]
  
  # Remove unknown batches
  unknown_count <- sum(meta[[batch_col]] == "Unknown" | is.na(meta[[batch_col]]))
  if (unknown_count > 0) {
    cat("WARN: Removing", unknown_count, "samples with Unknown/NA batch\n")
    valid_batch <- meta[[batch_col]] != "Unknown" & !is.na(meta[[batch_col]])
    meta <- meta[valid_batch, ]
    tab <- tab[rownames(meta), ]
  }
  
  # Remove samples with NA covariates
  complete_rows <- complete.cases(meta[, c(batch_col, covars)])
  n_incomplete <- sum(!complete_rows)
  
  if (n_incomplete > 0) {
    cat("WARN: Removing", n_incomplete, "samples with NA in covariates\n")
    meta <- meta[complete_rows, ]
    tab <- tab[rownames(meta), ]
  }
  
  # Check minimum samples
  if (nrow(meta) < 20) {
    cat("ERROR: Too few samples after filtering (", nrow(meta), ")\n")
    return(list(
      data = tab,
      phase_used = 0,
      success = FALSE,
      message = "Insufficient samples"
    ))
  }
  
  # Batch validation
  batch_counts <- table(meta[[batch_col]])
  cat("INFO: Batch distribution:\n")
  print(batch_counts)
  
  small_batches <- names(batch_counts)[batch_counts < 5]
  if (length(small_batches) > 0) {
    cat("WARN: Removing batches with <5 samples:", paste(small_batches, collapse=", "), "\n")
    keep_batches <- !(meta[[batch_col]] %in% small_batches)
    meta <- meta[keep_batches, ]
    tab <- tab[rownames(meta), ]
    batch_counts <- table(meta[[batch_col]])
  }
  
  # Covariate validation (batch-level)
  cat("\nINFO: Validating covariates at BATCH level...\n")
  
  batches <- unique(meta[[batch_col]])
  keep_covars <- rep(TRUE, length(covars))
  names(keep_covars) <- covars
  
  for (i in seq_along(covars)) {
    cov <- covars[i]
    
    if (!cov %in% colnames(meta)) {
      cat("WARN: Covariate", cov, "not found\n")
      keep_covars[i] <- FALSE
      next
    }
    
    problematic_batches <- c()
    
    for (batch in batches) {
      batch_data <- meta[meta[[batch_col]] == batch, cov]
      
      if (is.numeric(batch_data)) {
        batch_var <- var(batch_data, na.rm = TRUE)
        if (is.na(batch_var) || batch_var == 0) {
          problematic_batches <- c(problematic_batches, batch)
        }
      } else {
        batch_levels <- length(unique(batch_data[!is.na(batch_data)]))
        if (batch_levels < 2) {
          problematic_batches <- c(problematic_batches, batch)
        }
      }
    }
    
    if (length(problematic_batches) > 0) {
      cat("  Covariate", cov, "problematic in", length(problematic_batches), "batch(es)\n")
      cat("  Removing:", cov, "\n")
      keep_covars[i] <- FALSE
    }
  }
  
  valid_covars <- covars[keep_covars]
  
  if (length(valid_covars) == 0) {
    cat("\nWARN: No valid covariates after batch-level filtering\n")
    cat("  Continuing WITHOUT covariates\n\n")
    final_covars <- NULL
  } else {
    cat("INFO: Using", length(valid_covars), "covariate(s):", 
        paste(valid_covars, collapse=", "), "\n\n")
    final_covars <- meta[, valid_covars, drop = FALSE]
  }
  
  meta <- droplevels(meta)

  # Baseline evaluation
  cat("\n", strrep("-", 70), "\n")
  cat("BASELINE EVALUATION (Before correction)\n")
  cat(strrep("-", 70), "\n")

  baseline <- evaluate_batch_effect(tab, meta, batch_col, valid_covars)
  
  cat(sprintf("Batch R²:     %.4f (p = %.4f)\n", 
              baseline$batch_r2, baseline$batch_p))
  if (!is.na(baseline$covariate_r2) && baseline$covariate_r2 > 0) {
    cat(sprintf("Covariate R²: %.4f\n", baseline$covariate_r2))
  }
  cat(strrep("-", 70), "\n\n")
  
  # Select reference batches
  ref_batch_pool <- select_optimal_reference_batch(tab, meta, batch_col)
  
  # === PHASE 1 ===
  cat("\n>>> STARTING PHASE 1 <<<\n")
  phase1_result <- run_phase(1, tab, meta, batch_col, valid_covars, 
                              ref_batch_pool, final_covars)
  
  if (!phase1_result$success) {
    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! ConQuR PHASE 1 FAILED - RETURNING UNCORRECTED DATA !!!\n")
    cat(strrep("!", 70), "\n\n")
    return(list(
      data = tab,
      phase_used = 0,
      success = FALSE,
      message = phase1_result$error
    ))
  }
  
  # Evaluate Phase 1
  cat("\n>>> EVALUATING PHASE 1 <<<\n")
  phase1_eval <- evaluate_batch_effect(phase1_result$data, meta, 
                                        batch_col, valid_covars)
  
  delta_r2 <- baseline$batch_r2 - phase1_eval$batch_r2
  r2_reduction <- delta_r2 / baseline$batch_r2
  
  cat("\n", strrep("-", 70), "\n")
  cat("PHASE 1 EVALUATION\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Before:       Batch R² = %.4f (p = %.4f)\n", 
              baseline$batch_r2, baseline$batch_p))
  cat(sprintf("After:        Batch R² = %.4f (p = %.4f)\n", 
              phase1_eval$batch_r2, phase1_eval$batch_p))
  cat(sprintf("Change:       ΔR² = %.4f\n", delta_r2))
  cat(sprintf("Reduction:    %.1f%%\n", r2_reduction * 100))
  
  # Overcorrection check
  if (!is.na(baseline$covariate_r2) && baseline$covariate_r2 > 0) {
    covar_reduction <- (baseline$covariate_r2 - phase1_eval$covariate_r2) / 
                       baseline$covariate_r2
    if (covar_reduction > 0.2) {
      cat(sprintf("⚠ WARNING:    Covariate R² reduced by %.1f%% (possible overcorrection)\n", 
                  covar_reduction * 100))
    }
  }
  cat(strrep("-", 70), "\n\n")
  
  # Decision point
  if (phase1_eval$batch_r2 < baseline$batch_r2 && r2_reduction >= r2_threshold) {
    cat("✓ DECISION: STOP at Phase 1\n")
    cat(sprintf("  R² reduction (%.1f%%) meets threshold (≥%.1f%%)\n", 
                r2_reduction * 100, r2_threshold * 100))
    cat("  Phase 3 not needed.\n\n")
    
    # Save Phase 1 result
    final_result <- list(
      data = phase1_result$data,
      method_matrix = phase1_result$method_matrix,
      phase_used = 1,
      baseline_r2 = baseline$batch_r2,
      final_r2 = phase1_eval$batch_r2,
      r2_reduction = r2_reduction,
      elapsed = phase1_result$elapsed,
      success = TRUE
    )
    
  } else {
    # Need Phase 3
    if (phase1_eval$batch_r2 >= baseline$batch_r2) {
      cat("⚠ DECISION: RUN PHASE 3\n")
      cat("  Phase 1 did not improve batch effect (R² unchanged or worse)\n")
    } else {
      cat("→ DECISION: RUN PHASE 3\n")
      cat(sprintf("  Phase 1 R² reduction (%.1f%%) below threshold (%.1f%%)\n", 
                  r2_reduction * 100, r2_threshold * 100))
      cat("  Phase 3 may find better parameters.\n")
    }
    
    cat("\nWaiting 5 seconds before starting Phase 3...\n")
    Sys.sleep(5)
    
    # === PHASE 3 ===
    cat("\n>>> STARTING PHASE 3 <<<\n")
    phase3_result <- run_phase(3, tab, meta, batch_col, valid_covars, 
                                ref_batch_pool, final_covars)
    
    if (!phase3_result$success) {
      cat("\nPhase 3 failed. Using Phase 1 result.\n")
      final_result <- list(
        data = phase1_result$data,
        method_matrix = phase1_result$method_matrix,
        phase_used = 1,
        baseline_r2 = baseline$batch_r2,
        final_r2 = phase1_eval$batch_r2,
        r2_reduction = r2_reduction,
        elapsed = phase1_result$elapsed,
        success = TRUE,
        note = "Phase 3 failed, used Phase 1"
      )
    } else {
      # Evaluate Phase 3
      cat("\n>>> EVALUATING PHASE 3 <<<\n")
      phase3_eval <- evaluate_batch_effect(phase3_result$data, meta, 
                                            batch_col, valid_covars)
      
      phase3_delta <- baseline$batch_r2 - phase3_eval$batch_r2
      phase3_reduction <- phase3_delta / baseline$batch_r2
      
      cat("\n", strrep("-", 70), "\n")
      cat("PHASE 3 EVALUATION\n")
      cat(strrep("-", 70), "\n")
      cat(sprintf("Baseline:     Batch R² = %.4f\n", baseline$batch_r2))
      cat(sprintf("Phase 1:      Batch R² = %.4f (%.1f%% reduction)\n", 
                  phase1_eval$batch_r2, r2_reduction * 100))
      cat(sprintf("Phase 3:      Batch R² = %.4f (%.1f%% reduction)\n", 
                  phase3_eval$batch_r2, phase3_reduction * 100))
      cat(strrep("-", 70), "\n\n")
      
      # Compare and select
      if (phase3_eval$batch_r2 < phase1_eval$batch_r2) {
        improvement <- ((phase1_eval$batch_r2 - phase3_eval$batch_r2) / 
                       phase1_eval$batch_r2) * 100
        cat("✓ DECISION: Use Phase 3 result\n")
        cat(sprintf("  Phase 3 achieved %.1f%% lower Batch R² than Phase 1\n", 
                    improvement))
        
        final_result <- list(
          data = phase3_result$data,
          method_matrix = phase3_result$method_matrix,
          phase_used = 3,
          baseline_r2 = baseline$batch_r2,
          final_r2 = phase3_eval$batch_r2,
          r2_reduction = phase3_reduction,
          elapsed = phase1_result$elapsed + phase3_result$elapsed,
          success = TRUE,
          phase1_r2 = phase1_eval$batch_r2
        )
        
      } else {
        cat("✓ DECISION: Use Phase 1 result\n")
        cat("  Phase 3 did not improve upon Phase 1\n")
        
        final_result <- list(
          data = phase1_result$data,
          method_matrix = phase1_result$method_matrix,
          phase_used = 1,
          baseline_r2 = baseline$batch_r2,
          final_r2 = phase1_eval$batch_r2,
          r2_reduction = r2_reduction,
          elapsed = phase1_result$elapsed + phase3_result$elapsed,
          success = TRUE,
          note = "Phase 1 was better",
          phase3_r2 = phase3_eval$batch_r2
        )
      }
    }
  }
  
  cat("\n", strrep("=", 70), "\n")
  cat("AUTO-PHASE COMPLETED\n")
  cat(strrep("=", 70), "\n")
  cat("Phase used:       ", final_result$phase_used, "\n", sep="")
  cat(sprintf("Baseline R²:      %.4f\n", final_result$baseline_r2))
  cat(sprintf("Final R²:         %.4f\n", final_result$final_r2))
  cat(sprintf("Total reduction:  %.1f%%\n", final_result$r2_reduction * 100))
  cat(sprintf("Total time:       %.1f minutes\n", final_result$elapsed))
  cat(strrep("=", 70), "\n\n")

  return(final_result)
}

# ---------- PCoA plots --------------------------------------------------------
generate_conqur_plots <- function(tab, meta, batch_col, title, out_pref){
  smp <- intersect(rownames(tab), rownames(meta))
  fac <- meta[smp, batch_col]
  png(paste0(out_pref,"_bray.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac,
            main = paste(title,", Bray-Curtis"))
  dev.off()
  png(paste0(out_pref,"_aitch.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac, dissimilarity = "Aitch",
            main = paste(title,", Aitchison"))
  dev.off()
}

# ---------- Normalize ---------------------------------------------------------
normalize_all <- function(mat, stub){
  dir.create(dirname(stub), showWarnings = FALSE, recursive = TRUE)

  # Preserve original column and row names
  original_colnames <- colnames(mat)
  original_rownames <- rownames(mat)
  cat("\nINFO: Processing", nrow(mat), "samples x", length(original_colnames), "taxa\n")
  
  # CRITICAL: Ensure mat is a proper matrix/data.frame
  if (!is.data.frame(mat) && !is.matrix(mat)) {
    cat("  → Converting input to data.frame\n")
    mat <- as.data.frame(mat)
  }

  paths <- c(pseudo = paste0(stub,"_clr_pseudo.tsv"),
             czm    = paste0(stub,"_clr_czm.tsv"),
             rclr   = paste0(stub,"_rclr.tsv"))
  out <- list()

  ## CLR + pseudo
  cat("\n[1/3] CLR + pseudo-count\n")
  if (file.exists(paths["pseudo"])) {
    cat("  → Loading cached file\n")
    out$CLR_pseudo <- load_table(paths["pseudo"])
  } else {
    tryCatch({
      # Convert to matrix for numerical operations
      mat_num <- as.matrix(mat)
      positive_vals <- mat_num[mat_num > 0]
      
      if (length(positive_vals) == 0) {
        stop("No positive values in data")
      }
      
      pseudo_candidate <- quantile(positive_vals, 0.05) * 0.5
      pseudo <- max(pseudo_candidate, 0.1)
      
      cat("  → Calculated pseudocount:", round(pseudo, 6), "\n")
      
      # Add pseudocount and perform CLR
      mat_pseudo <- mat_num + pseudo
      result <- decostand(mat_pseudo, "clr")
      
      # Convert back to data.frame and restore names
      result <- as.data.frame(result)
      colnames(result) <- original_colnames
      rownames(result) <- original_rownames
      
      if (any(!is.finite(as.matrix(result)))) {
        stop("Non-finite values in CLR result")
      }
      
      out$CLR_pseudo <- result
      save_table(out$CLR_pseudo, paths["pseudo"])
      cat("  ✓ CLR+pseudo successful\n")
      
    }, error = function(e) {
      cat("  ✗ CLR+pseudo failed:", e$message, "\n")
      out$CLR_pseudo <<- NULL
    })
  }

  ## CLR + CZM
  cat("\n[2/3] CLR + CZM (zero imputation)\n")
  if (file.exists(paths["czm"])) {
    cat("  → Loading cached file\n")
    out$CLR_CZM <- load_table(paths["czm"])
  } else {
    tryCatch({
      # cmultRepl requires matrix input
      czm_mat <- cmultRepl(as.matrix(mat), 
                           method = "CZM", 
                           label = 0,
                           output = "p-counts",
                           z.warning = 1, 
                           z.delete = FALSE,
                           suppress.print = TRUE)
      
      # Perform CLR on imputed data
      result <- decostand(czm_mat, "clr")
      
      # Convert to data.frame and restore names
      result <- as.data.frame(result)
      colnames(result) <- original_colnames
      rownames(result) <- original_rownames
      
      out$CLR_CZM <- result   
      save_table(out$CLR_CZM, paths["czm"])
      cat("  ✓ CLR+CZM successful\n")
      
    }, error = function(e) {
      cat("  ✗ CLR+CZM failed:", e$message, "\n")
      out$CLR_CZM <<- NULL
    })
  }

  ## RCLR
  cat("\n[3/3] RCLR (robust CLR)\n")
  if (file.exists(paths["rclr"])) {
    cat("  → Loading cached file\n")
    out$RCLR <- load_table(paths["rclr"])
  } else {
    tryCatch({
      # rclr_f requires matrix input
      result <- rclr_f(as.matrix(mat))
      result[is.na(result)] <- 0
      
      # Convert to data.frame and restore names
      result <- as.data.frame(result)
      colnames(result) <- original_colnames
      rownames(result) <- original_rownames
      
      out$RCLR <- result   
      save_table(out$RCLR, paths["rclr"])
      cat("  ✓ RCLR successful\n")
      
    }, error = function(e) {
      cat("  ✗ RCLR failed:", e$message, "\n")
      out$RCLR <<- NULL
    })
  }
  
  # Remove NULL entries
  out <- out[!sapply(out, is.null)]
  
  cat("\n", strrep("-", 60), "\n")
  cat("Normalization summary:\n")
  cat("  Methods successful:", length(out), "/3\n")
  if (length(out) > 0) {
    cat("  Available methods:", paste(names(out), collapse=", "), "\n")
  }
  cat(strrep("-", 60), "\n")
  
  return(out)
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

parser <- ArgumentParser()
parser$add_argument("--otu",    required = TRUE)
parser$add_argument("--meta",   required = TRUE)
parser$add_argument("--prefix", required = TRUE)
parser$add_argument("--tumor_only", action = "store_true")
parser$add_argument("--phase", default = "auto",
                   help = "Phase mode: 'auto', '1', or '3'")
parser$add_argument("--r2_threshold", type = "double", default = 0.25,
                   help = "R2 reduction threshold for stopping at Phase 1 (default: 0.25)")
parser$add_argument("--batch_column", type = "character", default = "shipment_batch",
                   help = "Name of the column in metadata that specifies batch (default: shipment_batch)")
parser$add_argument("--covariates", type = "character", default = "tcga_patient_age",
                   help = "Comma-separated list of covariate column names (default: tcga_patient_age)")
parser$add_argument("--type_column", type = "character", default = "Type",
                   help = "Name of the column in metadata that specifies sample type for tumor filtering (default: Type)")
parser$add_argument("--tumor_value", type = "character", default = "Tumor",
                   help = "Comma-separated list of values to recognize as tumor samples (default: Tumor). Case-insensitive matching.")
parser$add_argument("--method", type = "character", default = "tune",
                   help = "Correction method: 'tune' (default, auto-phased Tune_ConQuR) or 'vanilla' (default ConQuR)")
args <- parser$parse_args()

dir.create(dirname(args$prefix), showWarnings = FALSE, recursive = TRUE)
corr_dir <- file.path(dirname(args$prefix), "corrected"); dir.create(corr_dir, showWarnings=FALSE)
norm_dir <- file.path(dirname(args$prefix), "normalized"); dir.create(norm_dir, showWarnings=FALSE)
plot_dir <- file.path(dirname(args$prefix), "pcoa_plots"); dir.create(plot_dir, showWarnings=FALSE)

# Parse covariates (comma-separated string to vector)
covariate_vector <- trimws(unlist(strsplit(args$covariates, ",")))
cat("INFO: Using batch column:", args$batch_column, "\n")
cat("INFO: Using covariates:", paste(covariate_vector, collapse = ", "), "\n")
if (args$tumor_only) {
  # Parse tumor values (comma-separated string to vector)
  tumor_value_vector <- trimws(unlist(strsplit(args$tumor_value, ",")))
  cat("INFO: Using type column for filtering:", args$type_column, "\n")
  cat("INFO: Tumor values to match:", paste(tumor_value_vector, collapse = ", "), "\n")
}

# Load data
data <- load_and_prepare_data(args$otu, args$meta)
otu_raw <- data$otu
meta <- data$meta

# Filter for tumor samples
if (args$tumor_only) {
  message("INFO: Filtering for 'Tumor' samples only.")

  if (!args$type_column %in% colnames(meta)) {
    stop("ERROR: Type column '", args$type_column, "' not found in metadata")
  }

  if (sum(is.na(meta[[args$type_column]])) > 0) {
    cat("WARN:", sum(is.na(meta[[args$type_column]])), "samples with NA in", args$type_column, "column\n")
    meta <- meta[!is.na(meta[[args$type_column]]), ]
    otu_raw <- otu_raw[rownames(meta), ]
  }

  meta[[args$type_column]] <- trimws(meta[[args$type_column]])

  # Case-insensitive matching with any of the tumor values
  tumor_samples <- rownames(meta[toupper(meta[[args$type_column]]) %in% toupper(tumor_value_vector), ])

  if (length(tumor_samples) == 0) {
    cat("ERROR: No tumor samples found.\n")
    cat("Available values in", args$type_column, "column:\n")
    print(table(meta[[args$type_column]]))
    cat("\nSearched for (case-insensitive):", paste(tumor_value_vector, collapse = ", "), "\n")
    stop("No tumor samples found")
  }

  meta <- meta[tumor_samples, , drop = FALSE]
  otu_raw <- otu_raw[tumor_samples, , drop = FALSE]

  message("INFO: ", nrow(meta), " tumor samples remaining.")
}

# Verify required columns
required_cols <- c(args$batch_column, covariate_vector)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols, collapse=", "))
}

predictors  <- c(args$batch_column, covariate_vector)
rhs_formula <- reformulate(predictors)

## Run ConQuR with auto-phase
tables <- list(raw = otu_raw)

cat("\n", strrep("=", 70), "\n")
cat("STARTING BATCH CORRECTION WITH AUTO-PHASE STRATEGY\n")
cat(strrep("=", 70), "\n")
cat("Mode:", args$phase, "\n")
cat("R² threshold:", args$r2_threshold * 100, "%\n")
cat(strrep("=", 70), "\n\n")

# Check cache
cache_file <- file.path(corr_dir, "ConQuR_tuned.tsv")
param_file <- file.path(corr_dir, "ConQuR_tuned_parameters.txt")

if (file.exists(cache_file)) {
  message("SKIP: ConQuR_tuned.tsv already exists (cached)")
  tables$ConQuR_tuned <- load_table(cache_file)
} else {
  # Dispatch based on --method argument
  if (args$method == "vanilla") {
    cat("INFO: Using VANILLA ConQuR method\n")
    result <- run_vanilla_conqur(
      otu_raw, meta, args$batch_column,
      predictors[-1],
      corr_dir = corr_dir
    )
  } else {
    # Default: auto-phased Tune_ConQuR
    result <- run_auto_phase(
      otu_raw, meta, args$batch_column,
      predictors[-1],
      corr_dir = corr_dir,
      r2_threshold = args$r2_threshold
    )
  }
  
  if (result$success) {
    tables$ConQuR_tuned <- result$data
    
    # Save results
    save_table(result$data, cache_file)
    
    # Save parameters
    sink(param_file)
    cat("Auto-Phase Tune_ConQuR Results\n")
    cat("==============================\n\n")
    cat("Phase used:       ", result$phase_used, "\n", sep="")
    cat("Total time:       ", round(result$elapsed, 1), " minutes\n", sep="")
    cat("\nQuality Metrics:\n")
    cat(sprintf("  Baseline Batch R²:     %.4f\n", result$baseline_r2))
    cat(sprintf("  Final Batch R²:        %.4f\n", result$final_r2))
    cat(sprintf("  R² reduction:          %.1f%%\n", result$r2_reduction * 100))
    if (!is.null(result$phase1_r2)) {
      cat(sprintf("  Phase 1 Batch R²:      %.4f\n", result$phase1_r2))
    }
    if (!is.null(result$phase3_r2)) {
      cat(sprintf("  Phase 3 Batch R²:      %.4f\n", result$phase3_r2))
    }
    cat("\nSelected parameters by taxa frequency:\n")
    cat("=====================================\n")
    print(result$method_matrix)
    cat("\n")
    if (!is.null(result$note)) {
      cat("Note:", result$note, "\n")
    }
    sink()
    
    cat("\n✓ Parameters saved to:", param_file, "\n")
    
  } else {
    # ==================== BATCH CORRECTION FAILED ====================
    cat("\n")
    cat(strrep("!", 70), "\n")
    cat("!!! ERROR: BATCH CORRECTION (ConQuR) FAILED !!!\n")
    cat(strrep("!", 70), "\n")
    cat("Reason:", ifelse(is.null(result$message), "Unknown error", result$message), "\n")
    cat(strrep("!", 70), "\n\n")

    # Save uncorrected data for pipeline continuity
    tables$ConQuR_tuned <- otu_raw
    save_table(otu_raw, cache_file)

    # Save failure info to parameters file
    sink(param_file)
    cat("Auto-Phase Tune_ConQuR Results\n")
    cat("==============================\n\n")
    cat("Status: FAILED\n")
    cat("Reason:", ifelse(is.null(result$message), "Unknown error", result$message), "\n")
    cat("\nUsing uncorrected data for downstream analysis.\n")
    sink()

    cat("WARNING: Saved UNCORRECTED data to corrected folder.\n")
    cat("         The 'ConQuR_tuned' results are identical to 'raw' data!\n\n")
  }
}

cat("\n", strrep("=", 70), "\n")
cat("BATCH CORRECTION COMPLETED\n")
cat(strrep("=", 70), "\n\n")

## PERMANOVA summary
sum_path <- file.path(dirname(args$prefix), "permanova_summary.tsv")
if (file.exists(sum_path)){
  message("SKIP: permanova_summary.tsv already exists")
  quit(status = 0)
}

## Main loop
permsum <- list()
for (tag in names(tables)){
  mat <- tables[[tag]]

  generate_conqur_plots(mat, meta, args$batch_column,
                        tag, file.path(plot_dir, tag))

  permsum[[length(permsum)+1]] <- run_permanova(
        mat, meta, rhs_formula, "bray") %>%
        mutate(Source = tag, Transform = "counts")

  norm_res <- normalize_all(mat,
        file.path(norm_dir, paste0(basename(args$prefix),"_",tag)))

  if (length(norm_res) > 0) {
    for (tr in names(norm_res)){
      if (!is.null(norm_res[[tr]])) {
        permsum[[length(permsum)+1]] <- run_permanova(
              norm_res[[tr]], meta, rhs_formula, "euclidean") %>%
              mutate(Source = tag, Transform = tr)
      }
    }
  } else {
    cat("WARN: All normalization failed for", tag, "\n")
  }
}

bind_rows(permsum) %>% save_table(sum_path)

cat("\n", strrep("=", 70), "\n")
cat("PIPELINE COMPLETED\n")
cat(strrep("=", 70), "\n")
cat("Corrected data:   ", cache_file, "\n", sep="")
cat("Parameters:       ", param_file, "\n", sep="")
cat("PERMANOVA results:", sum_path, "\n")
cat("PCoA plots:       ", plot_dir, "\n", sep="")
cat("\nPlease review PERMANOVA R² values to assess batch correction quality.\n")
cat(strrep("=", 70), "\n\n")

message("DONE")

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# Auto mode (basic, recommended)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_only

# Custom columns and covariates
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --batch_column "batch_id" \
#   --covariates "age,gender,BMI" \
#   --type_column "sample_type" \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor" \
#   --tumor_only

# Custom R² threshold (default 25%)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --r2_threshold 0.30  # 30% reduction required

# Manual Phase 1 only
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --phase 1

# Manual Phase 3 only
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --phase 3

# Multiple tumor types (case-insensitive)
# Rscript script.R \
#   --otu data.csv \
#   --meta meta.tsv \
#   --prefix output/results \
#   --tumor_value "Tumor,Primary Tumor,Recurrent Tumor,Metastatic" \
#   --tumor_only

# ==============================================================================

#!/usr/bin/env Rscript
# ==============================================================================
#  MaAsLin3 Differential Abundance Analysis
#
#  Runs MaAsLin3 on genus- or species-level corrected count data.
#  Uses TSS normalization + LOG transformation with linear models.
#  Falls back to MaAsLin2 if MaAsLin3 is not installed.
#
#  Usage:
#    Rscript scripts/run_maaslin.R \
#      --otu RESULTS/05_BATCH_CORRECTION/shipment_rerun/corrected/ConQuR_tuned.tsv \
#      --meta RESULTS/04_DECONTAMINATION/bracken_run/EAC_GEJ_METADATA.txt \
#      --prefix results/maaslin_genus \
#      --formula "Type + country + age_diag" \
#      --level genus
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
})

# ---------- Detect MaAsLin version -------------------------------------------
# MaAsLin3 package name is lowercase 'maaslin3'; MaAsLin2 is 'Maaslin2'
USE_MAASLIN3 <- FALSE
tryCatch({
  library(maaslin3)
  USE_MAASLIN3 <- TRUE
  cat("INFO: Using MaAsLin3 v", as.character(packageVersion("maaslin3")), "\n")
}, error = function(e) {
  tryCatch({
    library(Maaslin2)
    cat("INFO: MaAsLin3 not found, using MaAsLin2 v", as.character(packageVersion("Maaslin2")), "\n")
  }, error = function(e2) {
    stop("ERROR: Neither MaAsLin3 nor MaAsLin2 is installed")
  })
})

# ---------- Argument parsing --------------------------------------------------
parser <- ArgumentParser(description = "MaAsLin3/2 differential abundance analysis")
parser$add_argument("--otu", required = TRUE, help = "Corrected count table (TSV, samples x taxa)")
parser$add_argument("--meta", required = TRUE, help = "Metadata TSV")
parser$add_argument("--prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--formula", type = "character", default = "Type",
                    help = "Fixed effects formula (default: Type)")
parser$add_argument("--level", type = "character", default = "genus",
                    choices = c("genus", "species"),
                    help = "Taxonomic level for filtering (default: genus)")
parser$add_argument("--normalization", type = "character", default = "TSS",
                    help = "Normalization method (default: TSS)")
parser$add_argument("--transform", type = "character", default = "LOG",
                    help = "Transform method (default: LOG)")
parser$add_argument("--min_prevalence", type = "double", default = 0.10,
                    help = "Minimum prevalence filter (default: 0.10)")
parser$add_argument("--correction", type = "character", default = "BH",
                    help = "P-value correction method (default: BH)")
parser$add_argument("--max_significance", type = "double", default = 0.25,
                    help = "Max q-value for significance (default: 0.25)")
parser$add_argument("--cores", type = "integer", default = 1,
                    help = "Number of cores (default: 1)")
args <- parser$parse_args()

# ---------- Create output directory -------------------------------------------
outdir <- dirname(args$prefix)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# MaAsLin writes to a directory; use prefix basename as the output subdir
maaslin_outdir <- paste0(args$prefix, "_maaslin_output")
dir.create(maaslin_outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load data ---------------------------------------------------------
cat("INFO: Loading corrected count data...\n")

# Auto-detect CSV vs TSV
if (grepl("\\.csv$", args$otu, ignore.case = TRUE)) {
  otu_raw <- read_csv(args$otu, show_col_types = FALSE)
} else {
  otu_raw <- read_tsv(args$otu, show_col_types = FALSE)
}

# First column is sample_id
first_col <- colnames(otu_raw)[1]
otu_df <- otu_raw %>%
  column_to_rownames(first_col) %>%
  as.data.frame()

cat("INFO: OTU table:", nrow(otu_df), "samples x", ncol(otu_df), "taxa\n")

# Load metadata
cat("INFO: Loading metadata...\n")
meta_raw <- read_tsv(args$meta, show_col_types = FALSE)

# Auto-detect sample ID column
sample_id_candidates <- c("donor_id", "sample_id", "sampleid", "Sample.ID",
                           "Sample_ID", "SampleID", "patient", "Patient")
sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(meta_raw)][1]

if (is.na(sample_id_col)) {
  stop("ERROR: Cannot find sample ID column in metadata")
}

meta_df <- meta_raw %>%
  column_to_rownames(sample_id_col) %>%
  as.data.frame()

# Common samples
common <- intersect(rownames(otu_df), rownames(meta_df))
cat("INFO: Common samples:", length(common), "\n")

otu_df <- otu_df[common, ]
meta_df <- meta_df[common, , drop = FALSE]

# Auto-compute log_total_reads if formula needs it and column is missing
if ("Total_Reads" %in% colnames(meta_df) && !("log_total_reads" %in% colnames(meta_df))) {
  meta_df$log_total_reads <- log(meta_df$Total_Reads)
  cat("INFO: Computed log_total_reads from Total_Reads\n")
}

# Remove all-zero taxa
nonzero <- colSums(otu_df) > 0
otu_df <- otu_df[, nonzero]
cat("INFO: Non-zero taxa:", sum(nonzero), "\n")

# Ensure counts are numeric (round to integer for count data)
otu_df <- as.data.frame(lapply(otu_df, function(x) round(as.numeric(x))))
rownames(otu_df) <- common

# Parse formula variables
formula_vars <- all.vars(as.formula(paste("~", args$formula)))
cat("INFO: Formula variables:", paste(formula_vars, collapse = ", "), "\n")

missing_vars <- setdiff(formula_vars, colnames(meta_df))
if (length(missing_vars) > 0) {
  stop("ERROR: Missing metadata columns: ", paste(missing_vars, collapse = ", "))
}

# Remove samples with NA in formula variables
complete <- complete.cases(meta_df[, formula_vars, drop = FALSE])
if (sum(!complete) > 0) {
  cat("WARN: Removing", sum(!complete), "samples with NA in formula variables\n")
  otu_df <- otu_df[complete, ]
  meta_df <- meta_df[complete, , drop = FALSE]
}

# Set reference levels for categorical variables
# MaAsLin2 requires explicit reference for variables with >2 levels
for (var in formula_vars) {
  if (var %in% colnames(meta_df) && is.character(meta_df[[var]])) {
    meta_df[[var]] <- factor(meta_df[[var]])
    cat("INFO: Converted", var, "to factor with levels:", paste(levels(meta_df[[var]]), collapse=", "), "\n")
  }
}

# Set specific reference levels
if ("Type" %in% colnames(meta_df) && is.factor(meta_df$Type)) {
  if ("Control" %in% levels(meta_df$Type)) {
    meta_df$Type <- relevel(meta_df$Type, ref = "Control")
    cat("INFO: Set Type reference level to 'Control'\n")
  } else if ("Normal" %in% levels(meta_df$Type)) {
    meta_df$Type <- relevel(meta_df$Type, ref = "Normal")
    cat("INFO: Set Type reference level to 'Normal'\n")
  }
}

# For country, set most common as reference (or alphabetically first)
if ("country" %in% colnames(meta_df) && is.factor(meta_df$country)) {
  most_common <- names(sort(table(meta_df$country), decreasing = TRUE))[1]
  meta_df$country <- relevel(meta_df$country, ref = most_common)
  cat("INFO: Set country reference level to '", most_common, "'\n")
}

# For sex, set first level as reference
if ("sex" %in% colnames(meta_df) && is.factor(meta_df$sex)) {
  cat("INFO: sex levels:", paste(levels(meta_df$sex), collapse=", "), "\n")
}

# For tobacco/alcohol, set first level as reference
for (var in c("tobacco", "alcohol")) {
  if (var %in% colnames(meta_df) && is.factor(meta_df[[var]])) {
    cat("INFO:", var, "levels:", paste(levels(meta_df[[var]]), collapse=", "), "\n")
  }
}

cat("INFO: Final dataset:", nrow(otu_df), "samples x", ncol(otu_df), "taxa\n")

# ---------- Run MaAsLin -----------------------------------------------------
cat("\nINFO: Running MaAsLin...\n")
cat("INFO: Formula variables:", paste(formula_vars, collapse = ", "), "\n")
cat("INFO: Normalization:", args$normalization, "\n")
cat("INFO: Transform:", args$transform, "\n")
cat("INFO: Min prevalence:", args$min_prevalence, "\n")
cat("INFO: Correction:", args$correction, "\n")
cat("INFO: Max significance:", args$max_significance, "\n\n")

start_time <- Sys.time()

if (USE_MAASLIN3) {
  # MaAsLin3 API (lowercase function name)
  # Build reference list for Type variable
  ref_list <- NULL
  if ("Type" %in% formula_vars) {
    ref_level <- levels(meta_df$Type)[1]  # already set to Control/Normal above
    ref_list <- paste0("Type,", ref_level)
  }
  result <- maaslin3::maaslin3(
    input_data = otu_df,
    input_metadata = meta_df,
    output = maaslin_outdir,
    fixed_effects = formula_vars,
    reference = ref_list,
    normalization = args$normalization,
    transform = args$transform,
    min_prevalence = args$min_prevalence,
    correction = args$correction,
    max_significance = args$max_significance,
    cores = args$cores,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    verbosity = "INFO"
  )
} else {
  # MaAsLin2 API
  result <- Maaslin2(
    input_data = otu_df,
    input_metadata = meta_df,
    output = maaslin_outdir,
    fixed_effects = formula_vars,
    normalization = args$normalization,
    transform = args$transform,
    min_prevalence = args$min_prevalence,
    correction = args$correction,
    max_significance = args$max_significance,
    cores = args$cores,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
}

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\nINFO: MaAsLin completed in %.1f minutes\n", elapsed))

# ---------- Extract and format results ----------------------------------------
cat("\nINFO: Extracting results...\n")

# MaAsLin2 writes all_results.tsv and significant_results.tsv
all_results_path <- file.path(maaslin_outdir, "all_results.tsv")
sig_results_path <- file.path(maaslin_outdir, "significant_results.tsv")

if (file.exists(all_results_path)) {
  all_res <- read_tsv(all_results_path, show_col_types = FALSE)
} else {
  # MaAsLin3 may use different naming
  candidates <- list.files(maaslin_outdir, pattern = "results.*\\.tsv$", full.names = TRUE)
  if (length(candidates) > 0) {
    all_res <- read_tsv(candidates[1], show_col_types = FALSE)
  } else {
    cat("WARN: No results TSV found in output directory\n")
    all_res <- data.frame()
  }
}

if (nrow(all_res) > 0) {
  # Normalize column names: MaAsLin2 uses qval_individual/pval_individual;
  # older MaAsLin2 uses qval/pval
  if ("qval_individual" %in% colnames(all_res) && !("qval" %in% colnames(all_res))) {
    all_res$qval <- all_res$qval_individual
    all_res$pval <- all_res$pval_individual
  }

  # Filter to abundance model only (for compound testing in newer MaAsLin2)
  if ("model" %in% colnames(all_res)) {
    abundance_res <- all_res %>% filter(model == "abundance")
    cat("INFO: MaAsLin compound testing detected. Using abundance model results.\n")
    cat("INFO: Abundance results:", nrow(abundance_res), "rows; Prevalence results:",
        nrow(all_res %>% filter(model == "prevalence")), "rows\n")
  } else {
    abundance_res <- all_res
  }

  # Save full results with our naming convention
  full_path <- paste0(args$prefix, "_maaslin_full.tsv")
  write_tsv(abundance_res, full_path)
  cat("INFO: Full results saved to:", full_path, "\n")

  # Extract significant results per variable
  for (var in formula_vars) {
    # MaAsLin uses 'metadata' column for variable name and 'value' for level
    var_res <- abundance_res %>%
      filter(metadata == var)

    if (nrow(var_res) == 0) {
      # Try matching with level appended (e.g., TypeTumor)
      var_res <- abundance_res %>%
        filter(grepl(paste0("^", var), metadata))
    }

    # Get unique values for this variable
    var_values <- unique(var_res$value)

    for (val in var_values) {
      val_res <- var_res %>%
        filter(value == val)

      sig_val <- val_res %>%
        filter(qval < args$max_significance) %>%
        arrange(qval)

      level_name <- paste0(var, val)
      cat(sprintf("  %s: %d significant taxa (q < %.2f)\n",
                  level_name, nrow(sig_val), args$max_significance))

      if (nrow(sig_val) > 0) {
        sig_path <- paste0(args$prefix, "_maaslin_sig_", level_name, ".tsv")
        write_tsv(sig_val, sig_path)
        cat("    Saved to:", sig_path, "\n")

        # Print top 10
        cat("    Top 10:\n")
        top10 <- head(sig_val, 10)
        print_cols <- intersect(c("feature", "coef", "stderr", "pval", "qval"), colnames(top10))
        print(as.data.frame(top10[, print_cols]))
        cat("\n")
      }
    }
  }
}

# ---------- Summary -----------------------------------------------------------
summary_path <- paste0(args$prefix, "_maaslin_summary.txt")
sink(summary_path)
cat("MaAsLin Analysis Summary\n")
cat("========================\n\n")
cat("MaAsLin version:", ifelse(USE_MAASLIN3, "3", "2"), "\n")
cat("Formula variables:", paste(formula_vars, collapse = " + "), "\n")
cat("Level:", args$level, "\n")
cat("Samples:", nrow(otu_df), "\n")
cat("Taxa (after prevalence filter):", ncol(otu_df), "\n")
cat("Normalization:", args$normalization, "\n")
cat("Transform:", args$transform, "\n")
cat("Min prevalence:", args$min_prevalence, "\n")
cat("P-value correction:", args$correction, "\n")
cat("Max significance:", args$max_significance, "\n")
cat("Runtime:", round(elapsed, 1), "minutes\n\n")

if (nrow(abundance_res) > 0) {
  cat("Results per term (abundance model):\n")
  for (var in unique(abundance_res$metadata)) {
    var_data <- abundance_res %>% filter(metadata == var)
    for (val in unique(var_data$value)) {
      val_data <- var_data %>% filter(value == val)
      n_sig <- sum(val_data$qval < args$max_significance, na.rm = TRUE)
      cat(sprintf("  %s [%s]: %d significant taxa (q < %.2f)\n",
                  var, val, n_sig, args$max_significance))
    }
  }
}
sink()

cat("\nINFO: Summary saved to:", summary_path, "\n")
cat("INFO: MaAsLin output directory:", maaslin_outdir, "\n")
cat("DONE\n")

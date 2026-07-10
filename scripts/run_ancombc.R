#!/usr/bin/env Rscript
# ==============================================================================
#  ANCOM-BC2 Differential Abundance Analysis
#
#  Runs ANCOM-BC2 on genus- or species-level corrected count data.
#  Handles compositionality via bias correction.
#
#  Usage:
#    Rscript scripts/run_ancombc.R \
#      --otu RESULTS/05_BATCH_CORRECTION/corrected/ConQuR_tuned.tsv \
#      --meta metadata.tsv \
#      --prefix results/ancombc_genus \
#      --formula "Type" \
#      --level genus
# ==============================================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(TreeSummarizedExperiment)
  library(tidyverse)
  library(argparse)
})

# ---------- Argument parsing --------------------------------------------------
parser <- ArgumentParser(description = "ANCOM-BC2 differential abundance analysis")
parser$add_argument("--otu", required = TRUE, help = "Corrected count table (TSV, samples x taxa)")
parser$add_argument("--meta", required = TRUE, help = "Metadata TSV")
parser$add_argument("--prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--formula", type = "character", default = "Type",
                    help = "Fixed effects formula (default: Type)")
parser$add_argument("--level", type = "character", default = "genus",
                    choices = c("genus", "species"),
                    help = "Taxonomic level for filtering (default: genus)")
parser$add_argument("--p_adj_method", type = "character", default = "holm",
                    help = "P-value adjustment method (default: holm)")
parser$add_argument("--alpha", type = "double", default = 0.05,
                    help = "Significance threshold (default: 0.05)")
parser$add_argument("--prv_cut", type = "double", default = 0.10,
                    help = "Minimum prevalence fraction (default: 0.10)")
parser$add_argument("--lib_cut", type = "integer", default = 1000L,
                    help = "Minimum library size (default: 1000)")
args <- parser$parse_args()

# ---------- Create output directory -------------------------------------------
outdir <- dirname(args$prefix)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load data ---------------------------------------------------------
cat("INFO: Loading corrected count data...\n")

if (grepl("\\.csv$", args$otu, ignore.case = TRUE)) {
  otu_raw <- read_csv(args$otu, show_col_types = FALSE)
} else {
  otu_raw <- read_tsv(args$otu, show_col_types = FALSE)
}

first_col <- colnames(otu_raw)[1]
otu_mat <- otu_raw %>%
  column_to_rownames(first_col) %>%
  as.matrix()

cat("INFO: OTU table:", nrow(otu_mat), "samples x", ncol(otu_mat), "taxa\n")

# Load metadata
cat("INFO: Loading metadata...\n")
meta_raw <- read_tsv(args$meta, show_col_types = FALSE)

sample_id_candidates <- c("donor_id", "sample_id", "sampleid", "Sample.ID",
                           "Sample_ID", "SampleID", "patient", "Patient")
sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(meta_raw)][1]

if (is.na(sample_id_col)) {
  stop("ERROR: Cannot find sample ID column in metadata")
}

meta_df <- meta_raw %>%
  column_to_rownames(sample_id_col) %>%
  as.data.frame()

common <- intersect(rownames(otu_mat), rownames(meta_df))
cat("INFO: Common samples:", length(common), "\n")

otu_mat <- otu_mat[common, ]
meta_df <- meta_df[common, ]

# Remove all-zero taxa
nonzero <- colSums(otu_mat) > 0
otu_mat <- otu_mat[, nonzero]
cat("INFO: Non-zero taxa:", sum(nonzero), "\n")

otu_mat <- round(otu_mat)

# Remove samples with NA in formula variables
formula_vars <- all.vars(as.formula(paste("~", args$formula)))
cat("INFO: Formula variables:", paste(formula_vars, collapse = ", "), "\n")

missing_vars <- setdiff(formula_vars, colnames(meta_df))
if (length(missing_vars) > 0) {
  stop("ERROR: Missing metadata columns: ", paste(missing_vars, collapse = ", "))
}

complete <- complete.cases(meta_df[, formula_vars, drop = FALSE])
if (sum(!complete) > 0) {
  cat("WARN: Removing", sum(!complete), "samples with NA in formula variables\n")
  otu_mat <- otu_mat[complete, ]
  meta_df <- meta_df[complete, ]
}

cat("INFO: Final dataset:", nrow(otu_mat), "samples x", ncol(otu_mat), "taxa\n")

# ---------- Create phyloseq object --------------------------------------------
cat("\nINFO: Creating phyloseq object...\n")

otu_ps <- otu_table(t(otu_mat), taxa_are_rows = TRUE)
sam_ps <- sample_data(meta_df)
ps <- phyloseq(otu_ps, sam_ps)
cat("INFO: phyloseq object:", ntaxa(ps), "taxa,", nsamples(ps), "samples\n")

# ---------- Run ANCOM-BC2 -----------------------------------------------------
cat("\nINFO: Running ANCOM-BC2...\n")
cat("INFO: Formula: ~", args$formula, "\n")
cat("INFO: P-value adjustment:", args$p_adj_method, "\n")
cat("INFO: Alpha:", args$alpha, "\n")
cat("INFO: Prevalence cutoff:", args$prv_cut, "\n")
cat("INFO: Library size cutoff:", args$lib_cut, "\n\n")

start_time <- Sys.time()

result <- ancombc2(
  data = ps,
  fix_formula = args$formula,
  p_adj_method = args$p_adj_method,
  alpha = args$alpha,
  prv_cut = args$prv_cut,
  lib_cut = args$lib_cut,
  neg_lb = TRUE,
  pseudo_sens = TRUE,
  verbose = TRUE
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\nINFO: ANCOM-BC2 completed in %.1f minutes\n", elapsed))

# ---------- Extract results ---------------------------------------------------
cat("\nINFO: Extracting results...\n")

res_df <- result$res

full_path <- paste0(args$prefix, "_ancombc2_full.tsv")
write_tsv(res_df, full_path)
cat("INFO: Full results saved to:", full_path, "\n")

for (var in formula_vars) {
  matching_lfc  <- grep(paste0("^lfc_",  var), colnames(res_df), value = TRUE)
  matching_q    <- grep(paste0("^q_",    var), colnames(res_df), value = TRUE)
  matching_diff <- grep(paste0("^diff_", var), colnames(res_df), value = TRUE)

  for (k in seq_along(matching_lfc)) {
    lfc_k  <- matching_lfc[k]
    q_k    <- matching_q[k]
    diff_k <- matching_diff[k]

    if (!all(c(lfc_k, q_k, diff_k) %in% colnames(res_df))) next

    sig_df <- res_df %>%
      filter(!!sym(diff_k) == TRUE) %>%
      dplyr::select(taxon, !!sym(lfc_k), !!sym(q_k)) %>%
      arrange(!!sym(q_k))

    level_name <- gsub("^lfc_", "", lfc_k)
    cat(sprintf("  %s: %d significant taxa (q < %.2f)\n",
                level_name, nrow(sig_df), args$alpha))

    if (nrow(sig_df) > 0) {
      sig_path <- paste0(args$prefix, "_ancombc2_sig_", level_name, ".tsv")
      write_tsv(sig_df, sig_path)
      cat("    Saved to:", sig_path, "\n")
    }
  }
}

# ---------- Summary -----------------------------------------------------------
summary_path <- paste0(args$prefix, "_ancombc2_summary.txt")
sink(summary_path)
cat("ANCOM-BC2 Analysis Summary\n")
cat("==========================\n\n")
cat("Formula: ~", args$formula, "\n")
cat("Level:", args$level, "\n")
cat("Samples:", nrow(otu_mat), "\n")
cat("Taxa tested:", ncol(otu_mat), "\n")
cat("P-value adjustment:", args$p_adj_method, "\n")
cat("Alpha:", args$alpha, "\n")
cat("Prevalence cutoff:", args$prv_cut, "\n")
cat("Library size cutoff:", args$lib_cut, "\n")
cat("Runtime:", round(elapsed, 1), "minutes\n\n")

cat("Results per term:\n")
for (diff_col in grep("^diff_", colnames(res_df), value = TRUE)) {
  n_sig <- sum(res_df[[diff_col]], na.rm = TRUE)
  cat(sprintf("  %s: %d significant taxa\n", gsub("^diff_", "", diff_col), n_sig))
}
sink()

cat("\nINFO: Summary saved to:", summary_path, "\n")
cat("DONE\n")

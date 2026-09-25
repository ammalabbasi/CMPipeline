#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(vegan)})
a <- commandArgs(TRUE); repo <- a[1]; fixture <- a[2]; output <- a[3]
source(file.path(repo, "scripts/batch_contract.R"))
args <- list(taxon_column = "clade_name", sample_id_column = "patient")
input <- load_and_prepare_data(file.path(fixture, "counts.tsv"), file.path(fixture, "metadata.csv"))$otu
ids <- sprintf("S%02d", 1:24); raw <- as.matrix(input[ids, , drop = FALSE])
read_matrix <- function(file) {
  x <- as.data.frame(read_tsv(file, show_col_types = FALSE, name_repair = "minimal"))
  rn <- x[[1]]; x <- x[, -1, drop = FALSE]; rownames(x) <- rn
  matrix_ok(x); as.matrix(x)
}
corrected <- read_matrix(file.path(output, "corrected/ConQuR_tuned.tsv"))
assert_same_matrix_ids(corrected, raw)
status <- read_tsv(file.path(output, "corrected/BATCH_CORRECTION_STATUS.tsv"), show_col_types = FALSE)
stopifnot(status$status == "corrected", status$n_samples == 24L)
ledger <- read_tsv(file.path(output, "batch_correction_sample_status.tsv"), show_col_types = FALSE)
stopifnot(sum(ledger$status == "included") == 24L, sum(ledger$status == "excluded") == 4L)
for (source in c("raw", "ConQuR_tuned")) {
  counts <- if (source == "raw") raw else corrected
  for (method in c("clr_pseudo", "clr_czm", "rclr")) {
    file <- file.path(output, "normalized", paste0("batch_corr_run_", source, "_", method, ".tsv"))
    normalized <- read_matrix(file)
    assert_same_matrix_ids(normalized, raw)
    if (method == "clr_pseudo") {
      pseudo <- max(as.numeric(quantile(counts[counts > 0], .05)) * .5, .1)
      oracle <- log(counts + pseudo) - rowMeans(log(counts + pseudo))
      stopifnot(max(abs(normalized - oracle)) < 1e-10)
    }
  }
}
summary <- read_tsv(file.path(output, "permanova_summary.tsv"), show_col_types = FALSE)
stopifnot(setequal(unique(summary$Source), c("raw", "ConQuR_tuned")),
          setequal(unique(summary$Transform), c("counts", "CLR_pseudo", "CLR_CZM", "RCLR")))
cat("PASS: official batch entry; same 24 samples/taxa; four recorded exclusions; all six transforms; numerical CLR assertions\n")

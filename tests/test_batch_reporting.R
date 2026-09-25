#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
source(file.path(a[1], "scripts/batch_reporting.R"))
dir.create(a[2], recursive=TRUE, showWarnings=FALSE)
configure_batch_rng(17, 2)
tab <- matrix(c(1,2,3,4,5,6), 3, 2, dimnames=list(c("s1","s2","s3"), c("t1","t2")))
args <- list(prefix=file.path(a[2], "fixture"), seed=17, threads=2, method="tune", batch_column="batch")
mm <- matrix(c(NA, "TRUE", rep(NA,5)), ncol=1,
    dimnames=list(c("batch_ref","original","simple_match","logistic_lasso","quantile_type","lambda","interplt"),"0-1"))
write_batch_records(list(success=TRUE, baseline_r2=.4, final_r2=.4, method_matrix=mm), tab, c("batch","type"), args)
methods <- read.delim(file.path(a[2], "batch_correction_methods.tsv"))
quality <- read.delim(file.path(a[2], "batch_correction_quality.tsv"))
stopifnot(methods$selection == "original_selected", methods$taxa_in_bin == 2,
          quality$execution_status == "completed", quality$quality_status == "no_reduction_in_batch_R2")
cat("PASS: legitimate original selection is not relabelled a fitting failure\n")
write_batch_records(list(success=FALSE, message="injected failure"), tab, c("batch","type"), args)
methods <- read.delim(file.path(a[2], "batch_correction_methods.tsv"))
quality <- read.delim(file.path(a[2], "batch_correction_quality.tsv"))
stopifnot(methods$selection == "uncorrected_explicit_fallback",
          quality$execution_status == "uncorrected_explicit_fallback",
          quality$taxon_fit_diagnostics == "not_exposed_by_pinned_dependency")
cat("PASS: explicit uncorrected fallback is not labelled corrected\n")
cat("SUMMARY 2 reporting tests passed\n")

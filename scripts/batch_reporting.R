# Reproducible execution and diagnostic reporting, not a new correction algorithm.
configure_batch_rng <- function(seed, threads) {
  if (length(seed) != 1L || is.na(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be an integer between 0 and 2147483647")
  if (length(threads) != 1L || is.na(threads) || threads < 1 || threads != floor(threads))
    stop("threads must be a positive integer")
  if (.Platform$OS.type != "unix" && threads > 1L)
    stop("Reproducible parallel ConQuR currently requires the Unix multicore backend")
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(as.integer(seed))
  invisible(TRUE)
}

quality_summary <- function(baseline, final) {
  ok <- length(baseline) == 1L && length(final) == 1L && is.finite(baseline) && is.finite(final)
  ratio <- if (ok && baseline != 0) (baseline - final) / baseline else NA_real_
  label <- if (!ok) "not_evaluable_nonfinite_R2" else if (baseline <= 0)
    "not_interpretable_nonpositive_baseline" else if (final < 0)
    "not_interpretable_negative_final_R2" else if (final < baseline)
    "lower_batch_R2_descriptive_only" else "no_reduction_in_batch_R2"
  interpretable <- ok && baseline > 0 && final >= 0
  data.frame(baseline_R2 = if (length(baseline) == 1) baseline else NA_real_,
             final_R2 = if (length(final) == 1) final else NA_real_,
             diagnostic_relative_change = ratio,
             reported_reduction_percent = if (interpretable) ratio * 100 else NA_real_,
             quality_status = label, biological_signal_preservation = "not_assessed",
             taxon_fit_diagnostics = "not_exposed_by_pinned_dependency", stringsAsFactors = FALSE)
}

report_reduction_percent <- function(baseline, final) quality_summary(baseline, final)$reported_reduction_percent

write_batch_records <- function(result, tab, predictors, args) {
  out <- dirname(args$prefix)
  q <- quality_summary(result$baseline_r2, result$final_r2)
  q$execution_status <- if (isTRUE(result$success)) "completed" else "uncorrected_explicit_fallback"
  q$n_samples <- nrow(tab)
  q$n_taxa <- ncol(tab)
  q$method <- args$method
  write.table(q, file.path(out, "batch_correction_quality.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  settings <- data.frame(key = c("seed", "threads", "rng_kind", "normal_kind", "sample_kind", "covariates", "batch_column", "method", "comparison_scope"),
    value = c(args$seed, args$threads, RNGkind(), paste(predictors[-1], collapse = ","), args$batch_column, args$method,
              "same_inputs_code_environment_seed_threads; no cross-thread-count guarantee"))
  write.table(settings, file.path(out, "batch_correction_run_settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  mm <- result$method_matrix
  if (is.null(mm)) {
    methods <- data.frame(prevalence_bin = "all", selection = if (isTRUE(result$success)) "vanilla_requested" else "uncorrected_explicit_fallback",
                          taxa_in_bin = ncol(tab), fit_diagnostics = "not_exposed_by_pinned_dependency")
  } else {
    methods <- as.data.frame(t(mm), stringsAsFactors = FALSE)
    methods$prevalence_bin <- rownames(methods)
    prevalence <- colMeans(tab > 0)
    methods$taxa_in_bin <- vapply(strsplit(methods$prevalence_bin, "-", fixed = TRUE), function(bounds) {
      b <- as.numeric(bounds); sum(prevalence > b[1] & prevalence <= b[2])
    }, integer(1))
    methods$selection <- ifelse(methods$taxa_in_bin == 0, "empty_bin",
                          ifelse(is.na(methods$original), "unknown_selection",
                            ifelse(methods$original == "TRUE", "original_selected", "correction_selected")))
    methods$fit_diagnostics <- "not_exposed_by_pinned_dependency"
    rownames(methods) <- NULL
  }
  write.table(methods, file.path(out, "batch_correction_methods.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("QUALITY: ", q$quality_status, "; biological-signal preservation not assessed; per-taxon fit outcomes unavailable")
}

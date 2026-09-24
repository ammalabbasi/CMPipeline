#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(decontam)
  library(digest)
})

options(stringsAsFactors = FALSE, warn = 1)
script_started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)

stopf <- function(fmt, ...) {
  stop(sprintf(paste0("ERROR: ", fmt), ...), call. = FALSE)
}

collapse_values <- function(x, limit = 10L) {
  x <- unique(as.character(x))
  shown <- head(x, limit)
  suffix <- if (length(x) > limit) sprintf(" ... (%d total)", length(x)) else ""
  paste0(paste(shown, collapse = ", "), suffix)
}

parse_value_list <- function(x, label) {
  if (is.null(x) || length(x) != 1L || is.na(x) || trimws(x) == "") {
    stopf("--%s is required", label)
  }
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (any(values == "") || anyDuplicated(values)) {
    stopf("--%s must contain unique, nonempty comma-delimited values", label)
  }
  values
}

detect_delimiter <- function(path, label) {
  header <- readLines(path, n = 1L, warn = FALSE)
  if (length(header) != 1L || header == "") {
    stopf("%s is empty or has no header: %s", label, path)
  }
  tab_fields <- length(strsplit(header, "\t", fixed = TRUE)[[1]])
  comma_fields <- length(strsplit(header, ",", fixed = TRUE)[[1]])
  if (tab_fields < 2L && comma_fields < 2L) {
    stopf("Could not detect tab or comma delimiter for %s: %s", label, path)
  }
  if (tab_fields >= comma_fields) "\t" else ","
}

read_table_strict <- function(path, label) {
  if (!file.exists(path)) stopf("%s does not exist: %s", label, path)
  delimiter <- detect_delimiter(path, label)
  table <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = delimiter,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character",
      na.strings = character(),
      fill = FALSE,
      strip.white = FALSE
    ),
    error = function(e) stopf("Failed to parse %s: %s", label, conditionMessage(e))
  )
  if (nrow(table) == 0L) stopf("%s contains no data rows", label)
  if (any(names(table) == "")) stopf("%s contains a blank column name", label)
  if (anyDuplicated(names(table))) {
    duplicates <- unique(names(table)[duplicated(names(table))])
    stopf("%s contains duplicate sample/field column names: %s", label, collapse_values(duplicates))
  }
  table
}

validate_exact_strings <- function(x, label, ids = NULL) {
  missing <- is.na(x) | trimws(x) == "" | x %in% c("NA", "NaN")
  if (any(missing)) {
    affected <- if (is.null(ids)) which(missing) else ids[missing]
    stopf("%s contains missing values at: %s", label, collapse_values(affected))
  }
  padded <- x != trimws(x)
  if (any(padded)) {
    affected <- if (is.null(ids)) x[padded] else ids[padded]
    stopf("%s contains leading/trailing whitespace at: %s", label, collapse_values(affected))
  }
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

option_list <- list(
  make_option("--otu_table", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--prefix", type = "character"),
  make_option("--start_mode", type = "character"),
  make_option("--environment_spec", type = "character"),
  make_option("--nextflow_version", type = "character"),
  make_option("--taxon_column", type = "character"),
  make_option("--taxonomic_rank", type = "character"),
  make_option("--sample_id_column", type = "character"),
  make_option("--batch_column", type = "character"),
  make_option("--type_column", type = "character"),
  make_option("--control_role", type = "character"),
  make_option("--positive_values", type = "character"),
  make_option("--control_values", type = "character"),
  make_option("--unselected_type_policy", type = "character"),
  make_option("--threshold", type = "double"),
  make_option("--min_prevalence", type = "double"),
  make_option("--min_abundance", type = "integer"),
  make_option("--min_batches", type = "integer"),
  make_option("--min_total_per_batch", type = "integer"),
  make_option("--min_positive_per_batch", type = "integer"),
  make_option("--min_control_per_batch", type = "integer"),
  # keep | drop: see the zero-total block below. Not in required_args: defaults to keep.
  make_option("--zero_total_policy", type = "character", default = "keep"),
  # Optional TSV with columns `sample` and `bracken_genus_total` (from consensus_taxa):
  # the sample's library size before any filtering.
  make_option("--library_table", type = "character", default = NULL),
  # skip | error: a taxon present in a batch but in fewer than 2 of its samples cannot be tested
  # by decontam's prevalence method (it returns NA). skip = record it as unevaluable and not
  # flagged in that batch; error = stop, the ver4 behaviour.
  make_option("--unevaluable_policy", type = "character", default = "skip"),
  # Minimum upstream library (Bracken total) for an all-zero sample to count as evidence of
  # ABSENCE. Below it the sample could not have shown presence either way, so it is excluded
  # as insufficient_library and leaves the prevalence denominator (audit A06/A18; Ammal: 1000).
  make_option("--min_library", type = "double", default = 1000),
  # A zero-total sample whose library size is unknown (no --library_table, or a missing row)
  # stops the run unless this is set (audit A07).
  make_option("--allow_unknown_library", action = "store_true", default = FALSE)
)

args <- parse_args(OptionParser(option_list = option_list))

required_args <- c(
  "otu_table", "metadata", "prefix", "start_mode", "environment_spec",
  "nextflow_version", "taxon_column", "taxonomic_rank",
  "sample_id_column", "batch_column", "type_column", "control_role",
  "positive_values", "control_values", "unselected_type_policy", "threshold",
  "min_prevalence", "min_abundance", "min_batches", "min_total_per_batch",
  "min_positive_per_batch", "min_control_per_batch"
)
for (name in required_args) {
  value <- args[[name]]
  if (is.null(value) || length(value) != 1L || is.na(value) ||
      (is.character(value) && trimws(value) == "")) {
    stopf("--%s is required", name)
  }
}

if (!grepl("^[A-Za-z0-9._-]+$", args$prefix)) {
  stopf("--prefix may contain only letters, numbers, dot, underscore, and hyphen")
}
if (!(args$start_mode %in% c("beginning", "decontam"))) {
  stopf("Invalid --start_mode '%s'; expected beginning or decontam", args$start_mode)
}
if (!file.exists(args$environment_spec)) {
  stopf("--environment_spec does not exist: %s", args$environment_spec)
}
if (!(args$taxonomic_rank %in% c("genus", "species"))) {
  stopf("Invalid --taxonomic_rank '%s'; expected genus or species", args$taxonomic_rank)
}
if (!(args$control_role %in% c("experimental_blank", "surrogate_biological"))) {
  stopf("Invalid --control_role '%s'; expected experimental_blank or surrogate_biological", args$control_role)
}
if (!(args$unselected_type_policy %in% c("error", "exclude"))) {
  stopf("Invalid --unselected_type_policy '%s'; expected error or exclude", args$unselected_type_policy)
}
if (!is.finite(args$threshold) || args$threshold <= 0 || args$threshold > 1) {
  stopf("--threshold must be in (0, 1]")
}
if (!is.finite(args$min_prevalence) || args$min_prevalence <= 0 || args$min_prevalence > 1) {
  stopf("--min_prevalence must be in (0, 1]")
}
integer_args <- c("min_abundance", "min_batches", "min_total_per_batch",
                  "min_positive_per_batch", "min_control_per_batch")
for (name in integer_args) {
  if (args[[name]] < 1L) stopf("--%s must be at least 1", name)
}

positive_values <- parse_value_list(args$positive_values, "positive_values")
control_values <- parse_value_list(args$control_values, "control_values")
overlap_values <- intersect(positive_values, control_values)
if (length(overlap_values) > 0L) {
  stopf("Positive and control values overlap: %s", collapse_values(overlap_values))
}

otu_raw <- read_table_strict(args$otu_table, "OTU table")
metadata_raw <- read_table_strict(args$metadata, "metadata")

if (!(args$taxon_column %in% names(otu_raw))) {
  stopf("OTU table is missing taxon column '%s'", args$taxon_column)
}
sample_columns <- setdiff(names(otu_raw), args$taxon_column)
if (length(sample_columns) == 0L) stopf("OTU table contains no sample columns")
validate_exact_strings(sample_columns, "OTU sample column names")

taxon_ids <- otu_raw[[args$taxon_column]]
validate_exact_strings(taxon_ids, "Taxon IDs", seq_along(taxon_ids))
if (anyDuplicated(taxon_ids)) {
  duplicates <- unique(taxon_ids[duplicated(taxon_ids)])
  stopf("Duplicate taxon IDs: %s", collapse_values(duplicates))
}

count_columns <- lapply(sample_columns, function(sample_id) {
  raw_values <- otu_raw[[sample_id]]
  trimmed <- trimws(raw_values)
  numeric_values <- suppressWarnings(as.numeric(trimmed))
  invalid <- is.na(raw_values) | trimmed == "" | is.na(numeric_values) |
    !is.finite(numeric_values) | numeric_values < 0 |
    abs(numeric_values - round(numeric_values)) > .Machine$double.eps^0.5
  if (any(invalid)) {
    locations <- paste0(taxon_ids[invalid], "/", sample_id, "=", raw_values[invalid])
    error_label <- if (any(!is.na(numeric_values) & is.finite(numeric_values) & numeric_values < 0)) {
      "Negative count values"
    } else {
      "Invalid count values"
    }
    stopf("%s: %s", error_label, collapse_values(locations))
  }
  as.numeric(round(numeric_values))
})
count_matrix <- do.call(cbind, count_columns)
rownames(count_matrix) <- taxon_ids
colnames(count_matrix) <- sample_columns
storage.mode(count_matrix) <- "numeric"

required_metadata_columns <- c(args$sample_id_column, args$batch_column, args$type_column)
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata_raw))
if (length(missing_metadata_columns) > 0L) {
  stopf("Metadata is missing required columns: %s", collapse_values(missing_metadata_columns))
}
if ("decontam_role" %in% names(metadata_raw)) {
  stopf("Metadata contains reserved column 'decontam_role'")
}

metadata_ids <- metadata_raw[[args$sample_id_column]]
validate_exact_strings(metadata_ids, "Metadata sample IDs")
duplicate_metadata_ids <- unique(metadata_ids[duplicated(metadata_ids)])
if (length(duplicate_metadata_ids) > 0L) {
  conflict_records <- character()
  for (sample_id in duplicate_metadata_ids) {
    rows <- metadata_raw[metadata_ids == sample_id, , drop = FALSE]
    conflict_columns <- names(rows)[vapply(rows, function(x) length(unique(x)) > 1L, logical(1))]
    if (length(conflict_columns) > 0L) {
      conflict_records <- c(conflict_records,
                            sprintf("%s[%s]", sample_id, paste(conflict_columns, collapse = ",")))
    }
  }
  if (length(conflict_records) > 0L) {
    stopf("Conflicting duplicate metadata sample IDs: %s", collapse_values(conflict_records))
  }
  stopf("Duplicate metadata sample IDs: %s", collapse_values(duplicate_metadata_ids))
}

if (!(args$unevaluable_policy %in% c("skip", "error"))) {
  stopf("--unevaluable_policy must be skip or error (got '%s')", args$unevaluable_policy)
}
if (!(args$zero_total_policy %in% c("keep", "drop"))) {
  stopf("--zero_total_policy must be keep or drop (got '%s')", args$zero_total_policy)
}

# ZERO-TOTAL SAMPLES. A sample whose counts are all zero in THIS table (e.g. none of a short
# consensus whitelist is present in it) is evidence of absence for the prevalence method, not
# missing data: excluding zero-total controls makes the remaining controls look more
# contaminated than they are. So with policy `keep` such samples stay in, and only a truly EMPTY
# library (zero Bracken reads upstream, per --library_table) is excluded. Policy `drop`
# restores the earlier behaviour of excluding every zero-total sample.
# "Zero-total" is judged AFTER the --min_abundance filter (decision 24, audit A06): a near-empty
# library whose few reads all fall below min_abundance is all-zero in what the prevalence test
# sees, so it gets the same --min_library check as a sample that was all-zero on input.
sample_totals <- colSums(count_matrix * (count_matrix >= args$min_abundance))
zero_total_samples <- names(sample_totals)[sample_totals == 0]
library_totals <- NULL
if (!is.null(args$library_table) && nzchar(args$library_table)) {
  lib <- read_table_strict(args$library_table, "library table")
  if (all(c("sample", "bracken_genus_total") %in% names(lib))) {
    # bracken_library_sizes.tsv from consensus_taxa
    library_totals <- setNames(suppressWarnings(as.numeric(lib$bracken_genus_total)), lib$sample)
  } else {
    # An UNMASKED taxon table (first column = taxon, one column per sample): library = column sum.
    # Used when consensus is skipped, where no library-sizes file exists (audit A10/A07).
    lib_cols <- names(lib)[-1]
    lib_names <- sub("\\.bracken\\.[GS](\\.mpa)?\\.krakenreport\\.txt$", "", lib_cols)
    library_totals <- setNames(vapply(lib_cols, function(c) sum(suppressWarnings(as.numeric(lib[[c]])), na.rm = TRUE),
                                      numeric(1)), lib_names)
  }
}
if (!is.finite(args$min_library) || args$min_library < 0) stopf("--min_library must be >= 0")
insufficient_library_samples <- character()
if (args$zero_total_policy == "drop") {
  excluded_zero_samples <- zero_total_samples
} else {
  lib_size <- if (is.null(library_totals)) setNames(rep(NA_real_, length(zero_total_samples)), zero_total_samples)
              else library_totals[zero_total_samples]
  unknown <- zero_total_samples[is.na(lib_size)]
  if (length(unknown) > 0L && !args$allow_unknown_library) {
    stopf(paste0("%d zero-total sample(s) have no library size, so an empty or near-empty library cannot be ",
                 "told apart from real absence: %s. Pass --library_table (bracken_library_sizes.tsv from ",
                 "consensus; --decontam_library_sizes for --start_from decontam), or --allow_unknown_library"),
          length(unknown), collapse_values(unknown))
  }
  # Known and below --min_library (including 0): excluded. Unknown + --allow_unknown_library: kept.
  insufficient_library_samples <- zero_total_samples[!is.na(lib_size) & lib_size < args$min_library]
  excluded_zero_samples <- insufficient_library_samples
}
kept_zero_samples <- setdiff(zero_total_samples, excluded_zero_samples)
if (length(zero_total_samples) > 0L) {
  cat(sprintf("INFO: zero-total samples: %d (policy=%s): kept %d, excluded %d%s\n",
              length(zero_total_samples), args$zero_total_policy, length(kept_zero_samples),
              length(excluded_zero_samples),
              if (length(excluded_zero_samples)) paste0(" [", collapse_values(excluded_zero_samples), "]") else ""))
}
input_samples <- setdiff(names(sample_totals), excluded_zero_samples)
missing_metadata <- setdiff(input_samples, metadata_ids)
if (length(missing_metadata) > 0L) {
  stopf("OTU samples missing metadata: %s", collapse_values(missing_metadata))
}

matched_metadata <- metadata_raw[match(input_samples, metadata_ids), , drop = FALSE]
rownames(matched_metadata) <- input_samples
matched_types <- matched_metadata[[args$type_column]]
validate_exact_strings(matched_types, "Selected sample type", input_samples)

is_control <- matched_types %in% control_values
is_positive <- matched_types %in% positive_values
is_unselected <- !(is_control | is_positive)
if (any(is_unselected) && args$unselected_type_policy == "error") {
  unselected_values <- sort(unique(matched_types[is_unselected]))
  unselected_ids <- input_samples[is_unselected]
  stopf("Unselected sample types [%s] affect samples: %s",
        paste(unselected_values, collapse = ", "), collapse_values(unselected_ids))
}

selected_mask <- !is_unselected
selected_samples <- input_samples[selected_mask]
validated_metadata <- matched_metadata[selected_mask, , drop = FALSE]
validated_roles <- ifelse(is_control[selected_mask], "control", "positive")
validated_metadata$decontam_role <- validated_roles

if (length(selected_samples) == 0L) stopf("No samples remain after cohort selection")
if (!any(validated_roles == "control")) stopf("No selected controls remain after cohort selection")
if (!any(validated_roles == "positive")) stopf("No selected positives remain after cohort selection")

selected_batches <- validated_metadata[[args$batch_column]]
validate_exact_strings(selected_batches, "Selected sample batch", selected_samples)

batch_levels <- unique(selected_batches)
batch_eligibility <- do.call(rbind, lapply(batch_levels, function(batch_id) {
  in_batch <- selected_batches == batch_id
  n_positive <- sum(validated_roles[in_batch] == "positive")
  n_control <- sum(validated_roles[in_batch] == "control")
  n_total <- sum(in_batch)
  failures <- character()
  if (n_total < args$min_total_per_batch) {
    failures <- c(failures, sprintf("total_%d_below_minimum_%d", n_total, args$min_total_per_batch))
  }
  if (n_positive < args$min_positive_per_batch) {
    failures <- c(failures, sprintf("positive_%d_below_minimum_%d", n_positive, args$min_positive_per_batch))
  }
  if (n_control < args$min_control_per_batch) {
    failures <- c(failures, sprintf("control_%d_below_minimum_%d", n_control, args$min_control_per_batch))
  }
  data.frame(
    batch = batch_id,
    n_positive = n_positive,
    n_control = n_control,
    n_total = n_total,
    low_power_total_le_4 = n_total <= 4L,
    eligible = length(failures) == 0L,
    reason = if (length(failures) == 0L) "eligible" else paste(failures, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

eligible_batches <- batch_eligibility$batch[batch_eligibility$eligible]
if (length(eligible_batches) < args$min_batches) {
  write_tsv(batch_eligibility, sprintf("%s.batch_eligibility.tsv", args$prefix))
  details <- paste0(
    "batch ", batch_eligibility$batch, "(positive=", batch_eligibility$n_positive,
    ",control=", batch_eligibility$n_control, ",total=", batch_eligibility$n_total,
    ",eligible=", batch_eligibility$eligible, ",reason=", batch_eligibility$reason, ")"
  )
  stopf("Eligible batches %d below required %d: %s [before batching: %d zero-total sample(s) excluded, %d kept, policy=%s]",
        length(eligible_batches), args$min_batches, paste(details, collapse = "; "),
        length(excluded_zero_samples), length(kept_zero_samples), args$zero_total_policy)
}

selected_counts <- count_matrix[, selected_samples, drop = FALSE]
initial_taxon_sums <- rowSums(selected_counts)
initial_taxon_presence <- rowSums(selected_counts > 0)
filtered_counts <- selected_counts
filtered_counts[filtered_counts < args$min_abundance] <- 0
filtered_taxon_sums <- rowSums(filtered_counts)
filtered_taxon_presence <- rowSums(filtered_counts > 0)
filtered_prevalence <- filtered_taxon_presence / ncol(filtered_counts)
passes_prevalence <- filtered_taxon_sums > 0 & filtered_prevalence >= args$min_prevalence
evaluated_taxa <- taxon_ids[passes_prevalence]
if (length(evaluated_taxa) == 0L) {
  stopf("Zero taxa remain for contaminant evaluation after declared filters")
}

# decontam::isContaminant() removes every sample with zero total counts before testing ("Removed
# N samples with zero total counts"), because it first normalises counts to frequencies. For the
# prevalence method that throws away real evidence: an all-zero control shows those taxa ABSENT.
# With --zero_total_policy keep, a batch that contains an all-zero sample is therefore scored with
# decontam's own per-taxon prevalence test, decontam:::isContaminantPrevalence(), the function
# isContaminant() calls internally: same test, same p-values, zero samples retained. Other
# batches use isContaminant() unchanged. The self-check below proves the two agree whenever no
# sample is zero, and fails the run if a decontam release ever changes the internal function.
prevalence_test <- get0("isContaminantPrevalence", envir = asNamespace("decontam"), inherits = FALSE)
per_taxon_prevalence <- function(counts, is_neg) {
  p <- apply(counts, 1, function(x) prevalence_test(x, is_neg, method = "auto"))
  data.frame(p = p, p.prev = p, contaminant = !is.na(p) & p < args$threshold,
             row.names = rownames(counts), check.names = FALSE)
}
if (args$zero_total_policy == "keep") {
  if (is.null(prevalence_test)) {
    stopf("decontam %s has no internal isContaminantPrevalence(); --zero_total_policy keep cannot run. Use drop, or pin decontam 1.22.0",
          as.character(packageVersion("decontam")))
  }
  check_counts <- rbind(t1 = c(10, 20, 0, 30, 4, 0, 7, 0), t2 = c(0, 5, 8, 0, 9, 6, 0, 3),
                        t3 = c(12, 0, 0, 9, 0, 0, 5, 0), t4 = c(1, 1, 1, 1, 2, 2, 2, 2))
  colnames(check_counts) <- paste0("c", 1:8)
  check_neg <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  check_ps <- phyloseq::phyloseq(phyloseq::otu_table(check_counts, taxa_are_rows = TRUE),
                                 phyloseq::sample_data(data.frame(is.neg = check_neg,
                                                                  row.names = colnames(check_counts))))
  wrapper <- suppressWarnings(decontam::isContaminant(check_ps, method = "prevalence", neg = "is.neg",
                                                      threshold = args$threshold, detailed = TRUE))
  direct <- suppressWarnings(per_taxon_prevalence(check_counts, check_neg))
  if (!isTRUE(all.equal(as.numeric(wrapper$p.prev), as.numeric(direct$p.prev))) ||
      !identical(as.logical(wrapper$contaminant), as.logical(direct$contaminant))) {
    stopf("self-check failed: decontam::isContaminant and its per-taxon prevalence test disagree (decontam %s)",
          as.character(packageVersion("decontam")))
  }
}

contaminants_by_batch <- do.call(rbind, lapply(eligible_batches, function(batch_id) {
  batch_samples <- selected_samples[selected_batches == batch_id]
  batch_roles <- validated_roles[selected_batches == batch_id]
  batch_counts <- filtered_counts[evaluated_taxa, batch_samples, drop = FALSE]
  taxon_present <- rowSums(batch_counts) > 0
  otu <- phyloseq::otu_table(batch_counts, taxa_are_rows = TRUE)
  sample_frame <- data.frame(
    is.neg = batch_roles == "control",
    row.names = batch_samples,
    check.names = FALSE
  )
  physeq <- phyloseq::phyloseq(otu, phyloseq::sample_data(sample_frame))
  batch_has_zero_sample <- any(colSums(batch_counts) == 0)
  result <- tryCatch(
    if (args$zero_total_policy == "keep" && batch_has_zero_sample) {
      cat(sprintf("INFO: batch '%s' has %d all-zero sample(s); scored with decontam's per-taxon prevalence test so they count as absent\n",
                  batch_id, sum(colSums(batch_counts) == 0)))
      per_taxon_prevalence(batch_counts, batch_roles == "control")
    } else {
      decontam::isContaminant(
        physeq,
        method = "prevalence",
        neg = "is.neg",
        threshold = args$threshold,
        normalize = TRUE,
        detailed = TRUE
      )
    },
    error = function(e) stopf("decontam failed for batch '%s': %s", batch_id, conditionMessage(e))
  )
  result_taxa <- rownames(result)
  if (!setequal(result_taxa, evaluated_taxa)) {
    stopf("decontam returned a different taxon set for batch '%s'", batch_id)
  }
  result <- result[match(evaluated_taxa, result_taxa), , drop = FALSE]
  raw_p <- if ("p" %in% names(result)) as.numeric(result$p) else rep(NA_real_, length(evaluated_taxa))
  raw_p_prev <- if ("p.prev" %in% names(result)) as.numeric(result$p.prev) else rep(NA_real_, length(evaluated_taxa))
  raw_contaminant <- as.logical(result$contaminant)
  evaluation_status <- ifelse(
    !taxon_present,
    "absent_in_batch",
    ifelse(is.na(raw_contaminant) | is.na(raw_p), "unevaluable", "evaluated")
  )
  raw_contaminant[!taxon_present] <- FALSE
  data.frame(
    taxon_id = evaluated_taxa,
    batch = batch_id,
    n_positive = sum(batch_roles == "positive"),
    n_control = sum(batch_roles == "control"),
    p = raw_p,
    p_prev = raw_p_prev,
    evaluation_status = evaluation_status,
    contaminant = raw_contaminant,
    stringsAsFactors = FALSE
  )
}))
unevaluable_rows <- contaminants_by_batch$evaluation_status == "unevaluable"
if (any(unevaluable_rows) && args$unevaluable_policy == "skip") {
  # Untestable in that batch (present in < 2 of its samples) -> not flagged there. The final
  # call still needs >= --min_batches batches that DID flag it, so nothing is called on this.
  contaminants_by_batch$contaminant[unevaluable_rows] <- FALSE
  cat(sprintf("INFO: %d taxon/batch pair(s) unevaluable (present in < 2 samples of the batch); recorded, not flagged: %s\n",
              sum(unevaluable_rows),
              collapse_values(paste0(contaminants_by_batch$taxon_id[unevaluable_rows], "/",
                                     contaminants_by_batch$batch[unevaluable_rows]))))
}
if (any(unevaluable_rows) && args$unevaluable_policy == "error") {
  write_tsv(batch_eligibility, sprintf("%s.batch_eligibility.tsv", args$prefix))
  write_tsv(contaminants_by_batch, sprintf("%s.contaminants_by_batch.tsv", args$prefix))
  affected <- paste0(
    contaminants_by_batch$taxon_id[unevaluable_rows], "/",
    contaminants_by_batch$batch[unevaluable_rows]
  )
  stopf("Unevaluable decontam results for present taxa: %s", collapse_values(affected))
}

detected_batch_count <- vapply(evaluated_taxa, function(taxon_id) {
  sum(contaminants_by_batch$taxon_id == taxon_id & contaminants_by_batch$contaminant)
}, integer(1))
detected_batch_list <- vapply(evaluated_taxa, function(taxon_id) {
  paste(contaminants_by_batch$batch[
    contaminants_by_batch$taxon_id == taxon_id & contaminants_by_batch$contaminant
  ], collapse = ",")
}, character(1))
final_contaminant_flag <- detected_batch_count >= args$min_batches
final_contaminants <- evaluated_taxa[final_contaminant_flag]

# How many eligible batches could actually TEST each taxon (present in >= 2 samples there). With
# unevaluable_policy skip, a taxon testable in fewer than min_batches batches can never be called,
# however contaminant-like it is; `callable` makes that visible (audit A19).
n_evaluable_batches <- vapply(evaluated_taxa, function(taxon_id) {
  sum(contaminants_by_batch$taxon_id == taxon_id & contaminants_by_batch$evaluation_status == "evaluated")
}, integer(1))
contaminants_final <- data.frame(
  taxon_id = evaluated_taxa,
  detected_batch_count = detected_batch_count,
  detected_batches = detected_batch_list,
  min_batches = args$min_batches,
  n_evaluable_batches = n_evaluable_batches,
  callable = n_evaluable_batches >= args$min_batches,
  final_contaminant = final_contaminant_flag,
  stringsAsFactors = FALSE
)

final_taxa <- evaluated_taxa[!final_contaminant_flag]
final_counts <- filtered_counts[final_taxa, , drop = FALSE]
post_removal_zero_taxa <- rownames(final_counts)[rowSums(final_counts) == 0]
if (length(post_removal_zero_taxa) > 0L) {
  final_counts <- final_counts[setdiff(rownames(final_counts), post_removal_zero_taxa), , drop = FALSE]
  final_taxa <- rownames(final_counts)
}

filter_report <- data.frame(
  taxon_id = taxon_ids,
  initial_sum = initial_taxon_sums,
  initial_presence = initial_taxon_presence,
  after_abundance_sum = filtered_taxon_sums,
  after_abundance_presence = filtered_taxon_presence,
  after_abundance_prevalence = filtered_prevalence,
  passes_prevalence = passes_prevalence,
  detected_batch_count = 0L,
  detected_batches = "",
  final_action = "remove",
  terminal_reason = "",
  stringsAsFactors = FALSE
)
eval_index <- match(evaluated_taxa, filter_report$taxon_id)
filter_report$detected_batch_count[eval_index] <- detected_batch_count
filter_report$detected_batches[eval_index] <- detected_batch_list
for (i in seq_len(nrow(filter_report))) {
  taxon_id <- filter_report$taxon_id[i]
  if (initial_taxon_sums[i] == 0) {
    filter_report$terminal_reason[i] <- "all_zero_input"
  } else if (!passes_prevalence[i]) {
    filter_report$terminal_reason[i] <- "below_min_prevalence_after_abundance_filter"
  } else if (taxon_id %in% final_contaminants) {
    filter_report$terminal_reason[i] <- sprintf("contaminant_in_at_least_%d_batches", args$min_batches)
  } else if (taxon_id %in% post_removal_zero_taxa) {
    filter_report$terminal_reason[i] <- "all_zero_after_declared_filters"
  } else {
    filter_report$final_action[i] <- "retain"
    if (filter_report$detected_batch_count[i] > 0L) {
      filter_report$terminal_reason[i] <- sprintf(
        "contaminant_recurrence_below_%d_batches", args$min_batches
      )
    } else {
      filter_report$terminal_reason[i] <- "retained_noncontaminant"
    }
  }
}

all_validation_ids <- c(sample_columns, metadata_ids[!(metadata_ids %in% sample_columns)])
sample_validation <- data.frame(
  sample_id = all_validation_ids,
  in_otu = all_validation_ids %in% sample_columns,
  in_metadata = all_validation_ids %in% metadata_ids,
  initial_total_count = NA_real_,
  sample_type = NA_character_,
  role = NA_character_,
  final_total_count = NA_real_,
  final_status = NA_character_,
  reason = NA_character_,
  stringsAsFactors = FALSE
)
otu_match <- match(sample_validation$sample_id, sample_columns)
meta_match <- match(sample_validation$sample_id, metadata_ids)
sample_validation$initial_total_count[!is.na(otu_match)] <- sample_totals[otu_match[!is.na(otu_match)]]
sample_validation$sample_type[!is.na(meta_match)] <- metadata_raw[[args$type_column]][meta_match[!is.na(meta_match)]]
for (i in seq_len(nrow(sample_validation))) {
  sample_id <- sample_validation$sample_id[i]
  if (sample_id %in% excluded_zero_samples) {
    sample_validation$final_status[i] <- "excluded"
    sample_validation$reason[i] <- if (args$zero_total_policy == "drop") "zero_total_input" else "insufficient_library"
  } else if (!(sample_id %in% sample_columns)) {
    sample_validation$final_status[i] <- "excluded"
    sample_validation$reason[i] <- "metadata_only"
  } else if (!(sample_id %in% selected_samples)) {
    sample_validation$final_status[i] <- "excluded"
    sample_validation$reason[i] <- "unselected_sample_type"
  } else {
    selected_index <- match(sample_id, selected_samples)
    sample_validation$role[i] <- validated_roles[selected_index]
    sample_validation$final_total_count[i] <- colSums(final_counts)[selected_index]
    sample_validation$final_status[i] <- "included"
    sample_validation$reason[i] <- paste0("included_", validated_roles[selected_index])
  }
}

# Which included samples were all-zero in the input table (kept as absence evidence).
sample_validation$zero_total_input <- sample_validation$sample_id %in% zero_total_samples

output_files <- c(
  decontaminated = sprintf("%s.decontam_pkg_decontaminated.csv", args$prefix),
  validated_metadata = sprintf("%s.validated_metadata.tsv", args$prefix),
  sample_validation = sprintf("%s.sample_validation.tsv", args$prefix),
  batch_eligibility = sprintf("%s.batch_eligibility.tsv", args$prefix),
  contaminants_by_batch = sprintf("%s.contaminants_by_batch.tsv", args$prefix),
  contaminants_final = sprintf("%s.contaminants_final.tsv", args$prefix),
  filter_report = sprintf("%s.filter_report.tsv", args$prefix),
  run_manifest = sprintf("%s.run_manifest.tsv", args$prefix),
  summary = sprintf("%s.decontamination_summary.txt", args$prefix)
)
plot_directory <- "decontamination_plots"
plot_file <- file.path(plot_directory, sprintf("%s.decontamination_diagnostics.pdf", args$prefix))
existing_outputs <- c(output_files[file.exists(output_files)],
                      if (dir.exists(plot_directory)) plot_directory else character())
if (length(existing_outputs) > 0L) {
  stopf("Refusing to overwrite existing outputs: %s", collapse_values(existing_outputs))
}

decontaminated_table <- data.frame(
  setNames(list(final_taxa), args$taxon_column),
  as.data.frame(final_counts, check.names = FALSE),
  check.names = FALSE
)
write.csv(decontaminated_table, output_files[["decontaminated"]], row.names = FALSE, quote = TRUE)
write_tsv(validated_metadata, output_files[["validated_metadata"]])
write_tsv(sample_validation, output_files[["sample_validation"]])
write_tsv(batch_eligibility, output_files[["batch_eligibility"]])
write_tsv(contaminants_by_batch, output_files[["contaminants_by_batch"]])
write_tsv(contaminants_final, output_files[["contaminants_final"]])
write_tsv(filter_report, output_files[["filter_report"]])

dir.create(plot_directory, recursive = FALSE)
grDevices::pdf(plot_file, width = 8, height = 6, onefile = TRUE)
tryCatch({
  batch_counts <- rbind(batch_eligibility$n_positive, batch_eligibility$n_control)
  colnames(batch_counts) <- batch_eligibility$batch
  barplot(batch_counts, beside = FALSE, col = c("#B54A35", "#2E6F9E"),
          xlab = "Batch", ylab = "Selected samples",
          main = "Decontamination batch eligibility")
  legend("topright", legend = c("Positive", "Control"),
         fill = c("#B54A35", "#2E6F9E"), bty = "n")
  detected_per_batch <- tapply(contaminants_by_batch$contaminant,
                               contaminants_by_batch$batch, sum)
  barplot(detected_per_batch, col = "#C8902F", xlab = "Eligible batch",
          ylab = "Contaminant calls", main = "Per-batch contaminant calls")
}, finally = {
  grDevices::dev.off()
})

summary_lines <- c(
  "DECONTAMINATION SUMMARY",
  "=======================",
  sprintf("Control role: %s", args$control_role),
  sprintf("Taxonomic rank: %s", args$taxonomic_rank),
  sprintf("Input taxa: %d", nrow(count_matrix)),
  sprintf("Input samples: %d", ncol(count_matrix)),
  sprintf("Zero-total samples: %d (policy %s): kept %d as absence evidence, excluded %d",
          length(zero_total_samples), args$zero_total_policy, length(kept_zero_samples),
          length(excluded_zero_samples)),
  sprintf("Unevaluable taxon/batch pairs (policy %s): %d", args$unevaluable_policy,
          sum(contaminants_by_batch$evaluation_status == "unevaluable")),
  sprintf("Metadata-only samples excluded: %d", sum(sample_validation$reason == "metadata_only")),
  sprintf("Unselected-type samples excluded: %d", sum(sample_validation$reason == "unselected_sample_type")),
  sprintf("Selected samples: %d", length(selected_samples)),
  sprintf("Selected positives: %d", sum(validated_roles == "positive")),
  sprintf("Selected controls: %d", sum(validated_roles == "control")),
  sprintf("Eligible batches: %d", length(eligible_batches)),
  sprintf("Taxa after abundance/prevalence filters: %d", length(evaluated_taxa)),
  sprintf("Final contaminants: %d", length(final_contaminants)),
  sprintf("Final retained taxa: %d", nrow(final_counts))
)
writeLines(summary_lines, output_files[["summary"]])

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else NA_character_
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
conda_history <- if (nzchar(conda_prefix)) file.path(conda_prefix, "conda-meta", "history") else ""
script_completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest <- data.frame(
  key = c(
    "status", "prefix", "start_mode", "taxonomic_rank", "taxon_column", "sample_id_column",
    "batch_column", "type_column", "control_role", "positive_values", "control_values",
    "unselected_type_policy", "zero_total_policy", "unevaluable_policy", "min_library",
    "allow_unknown_library", "library_table_sha256",
    "threshold", "min_prevalence", "min_abundance",
    "min_batches", "min_total_per_batch", "min_positive_per_batch",
    "min_control_per_batch", "otu_sha256", "metadata_sha256", "script_sha256",
    "environment_spec_sha256", "environment_identifier", "conda_history_sha256",
    "nextflow_version", "slurm_job_id", "r_version", "decontam_version",
    "phyloseq_version", "optparse_version", "digest_version", "started_utc", "completed_utc"
  ),
  value = c(
    "completed", args$prefix, args$start_mode, args$taxonomic_rank, args$taxon_column,
    args$sample_id_column, args$batch_column, args$type_column, args$control_role,
    paste(positive_values, collapse = ","), paste(control_values, collapse = ","),
    args$unselected_type_policy, args$zero_total_policy, args$unevaluable_policy, args$min_library,
    args$allow_unknown_library,
    if (!is.null(args$library_table) && nzchar(args$library_table)) sha256_file(args$library_table) else "not_provided",
    args$threshold, args$min_prevalence,
    args$min_abundance, args$min_batches, args$min_total_per_batch,
    args$min_positive_per_batch, args$min_control_per_batch,
    sha256_file(args$otu_table), sha256_file(args$metadata),
    if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_,
    sha256_file(args$environment_spec),
    if (nzchar(conda_prefix)) basename(conda_prefix) else "not_set",
    if (nzchar(conda_history) && file.exists(conda_history)) sha256_file(conda_history) else "not_available",
    args$nextflow_version, Sys.getenv("SLURM_JOB_ID", unset = "not_set"),
    R.version.string, as.character(packageVersion("decontam")),
    as.character(packageVersion("phyloseq")), as.character(packageVersion("optparse")),
    as.character(packageVersion("digest")), script_started_utc, script_completed_utc
  ),
  stringsAsFactors = FALSE
)
hashed_outputs <- output_files[names(output_files) != "run_manifest"]
manifest <- rbind(
  manifest,
  data.frame(
    key = paste0("output_sha256.", names(hashed_outputs)),
    value = vapply(hashed_outputs, sha256_file, character(1)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    key = "output_sha256.diagnostics_pdf",
    value = sha256_file(plot_file),
    stringsAsFactors = FALSE
  )
)
write_tsv(manifest, output_files[["run_manifest"]])

cat(sprintf(
  "Decontamination completed: %d selected samples, %d eligible batches, %d final contaminants, %d retained taxa\n",
  length(selected_samples), length(eligible_batches), length(final_contaminants), nrow(final_counts)
))

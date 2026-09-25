# Exact input/cohort/output contracts for the batch module.
read_batch_table <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  if (length(header) != 1L || !nzchar(header)) stop("Empty table: ", path)
  delim <- if (grepl("\t", header, fixed = TRUE)) "\t" else ","
  x <- readr::read_delim(path, delim = delim, quote = '"',
                        col_types = readr::cols(.default = readr::col_character()),
                        name_repair = "minimal", na = character(), trim_ws = FALSE,
                        skip_empty_rows = FALSE, show_col_types = FALSE, progress = FALSE)
  if (nrow(readr::problems(x))) stop("Malformed table records: ", path)
  ids_ok(names(x), "column names")
  as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
}

ids_ok <- function(x, what) {
  if (!length(x) || anyNA(x) || any(!nzchar(trimws(x))) ||
      any(x != trimws(x)) || anyDuplicated(x))
    stop("Invalid or duplicate ", what)
  invisible(TRUE)
}

matrix_ok <- function(x, counts = FALSE, integer_counts = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("Matrix must be nonempty, numeric and finite")
  ids_ok(rownames(x), "sample IDs")
  ids_ok(colnames(x), "taxon IDs")
  if (counts && any(x < 0)) stop("Counts must be nonnegative")
  if (integer_counts && (any(x != floor(x)) || any(x > 2^53)))
    stop("Input counts must be exactly representable nonnegative integers")
  invisible(TRUE)
}

load_and_prepare_data <- function(otu_path, meta_path) {
  otu <- read_batch_table(otu_path)
  tax_col <- args$taxon_column
  id_col <- args$sample_id_column
  if (is.null(tax_col) || !nzchar(tax_col) || !tax_col %in% names(otu))
    stop("Declared taxon column not found: ", tax_col)
  ids_ok(otu[[tax_col]], "taxon IDs")
  samples <- setdiff(names(otu), tax_col)
  ids_ok(samples, "sample IDs")
  vals <- as.matrix(otu[, samples, drop = FALSE])
  numbers <- suppressWarnings(as.numeric(vals))
  mat <- matrix(numbers, nrow = nrow(otu), ncol = length(samples),
                dimnames = list(otu[[tax_col]], samples))
  mat <- t(mat)
  matrix_ok(mat, counts = TRUE, integer_counts = TRUE)
  meta <- read_batch_table(meta_path)
  if (is.null(id_col) || !nzchar(id_col) || !id_col %in% names(meta))
    stop("Declared sample ID column not found: ", id_col)
  ids_ok(meta[[id_col]], "metadata sample IDs")
  missing <- setdiff(samples, meta[[id_col]])
  if (length(missing)) stop("Metadata missing for count-table samples: ", paste(missing, collapse = ", "))
  extra <- setdiff(meta[[id_col]], samples)
  if (length(extra)) message("Unused metadata rows: ", paste(extra, collapse = ", "))
  meta <- meta[match(samples, meta[[id_col]]), , drop = FALSE]
  rownames(meta) <- meta[[id_col]]
  list(otu = as.data.frame(mat, check.names = FALSE), meta = meta, unused_metadata = extra)
}

prepare_batch_cohort <- function(data, args, covars) {
  tab <- data$otu; meta <- data$meta
  if (anyDuplicated(covars) || args$batch_column %in% covars || any(!nzchar(covars)))
    stop("Covariates must be unique, nonempty and distinct from the batch column")
  needed <- c(args$batch_column, covars, if (args$tumor_only) args$type_column)
  missing <- setdiff(needed, names(meta))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  ledger <- data.frame(sample = rownames(tab), status = "included", reason = "eligible",
                       stringsAsFactors = FALSE)
  exclude <- function(mask, reason) {
    ix <- which(mask & ledger$status == "included")
    ledger$status[ix] <<- "excluded"; ledger$reason[ix] <<- reason
  }
  missing_value <- function(x) is.na(x) | trimws(x) %in% c("", "NA")
  exclude(rowSums(tab) == 0, "all-zero after decontamination")
  if (args$tumor_only) {
    types <- meta[[args$type_column]]
    exclude(missing_value(types), "missing sample type")
    wanted <- toupper(trimws(strsplit(args$tumor_value, ",", fixed = TRUE)[[1]]))
    exclude(!toupper(trimws(types)) %in% wanted, "outside requested tumor subset")
  }
  batch <- trimws(meta[[args$batch_column]])
  exclude(missing_value(batch) | batch == "Unknown", "missing or Unknown batch")
  meta[[args$batch_column]] <- factor(batch)
  for (v in covars) {
    values <- trimws(meta[[v]])
    exclude(missing_value(values), paste0("missing covariate: ", v))
    values[missing_value(values)] <- NA_character_
    values <- type.convert(values, as.is = TRUE)
    if (is.numeric(values) && any(!is.finite(values[!is.na(values)])))
      stop("Non-finite covariate: ", v)
    meta[[v]] <- if (is.character(values) || is.logical(values)) factor(values) else values
  }
  sizes <- table(droplevels(meta[ledger$status == "included", args$batch_column]))
  small <- names(sizes)[sizes < 5L]
  exclude(as.character(meta[[args$batch_column]]) %in% small, "batch has fewer than 5 eligible samples")
  if (length(data$unused_metadata)) {
    ledger <- rbind(ledger, data.frame(sample = data$unused_metadata,
                         status = "unused", reason = "metadata row absent from count table"))
  }
  write.table(ledger, file.path(dirname(args$prefix), "batch_correction_sample_status.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  keep <- ledger$sample[ledger$status == "included"]
  if (!length(keep)) stop("No eligible samples after batch-cohort selection")
  tab <- tab[keep, , drop = FALSE]; meta <- droplevels(meta[keep, , drop = FALSE])
  if (nlevels(meta[[args$batch_column]]) < 2L) stop("At least two eligible batches are required")
  # Do not change the requested biological adjustment model to make a fit succeed.
  for (v in covars) {
    for (b in levels(meta[[args$batch_column]])) {
      x <- meta[meta[[args$batch_column]] == b, v]
      if (length(unique(x)) < 2L)
        stop("Requested covariate does not vary within batch: ", v, " / ", b,
             "; change the model explicitly rather than silently dropping the covariate")
    }
  }
  zero_taxa <- colSums(tab) == 0
  tab <- tab[, !zero_taxa, drop = FALSE]
  matrix_ok(tab, counts = TRUE, integer_counts = TRUE)
  write.table(data.frame(taxon = colnames(data$otu),
                         status = ifelse(colnames(data$otu) %in% colnames(tab), "included", "all-zero in analysis cohort")),
              file.path(dirname(args$prefix), "batch_correction_taxon_status.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  assign("dropped_empty_samples", ledger$sample[ledger$reason == "all-zero after decontamination"], envir = globalenv())
  list(otu = tab, meta = meta)
}

assert_same_matrix_ids <- function(x, reference) {
  if (!identical(rownames(x), rownames(reference)) || !identical(colnames(x), colnames(reference)))
    stop("Output sample/taxon IDs or order differ from the analysis input")
  invisible(TRUE)
}

batch_distance <- function(mat, meta, dist_meth, transformed = FALSE) {
  matrix_ok(mat, counts = !transformed)
  if (!identical(rownames(mat), rownames(meta))) stop("Count/metadata sample order differs")
  if (!transformed) {
    if (any(rowSums(mat) <= 0)) stop("Zero-total analysis sample cannot be evaluated")
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
  }
  if (!ncol(mat)) stop("No features available for distance calculation")
  vegan::vegdist(mat, method = dist_meth)
}

reject_existing_batch_outputs <- function(root) {
  dirs <- file.path(root, c("corrected", "normalized", "pcoa_plots"))
  existing <- unlist(lapply(dirs, list.files, all.files = TRUE, no.. = TRUE, recursive = TRUE))
  if (length(existing) || any(file.exists(file.path(root, c("permanova_summary.tsv",
                                  "batch_correction_sample_status.tsv", "batch_correction_taxon_status.tsv", "batch_correction_quality.tsv",
                                  "batch_correction_methods.tsv", "batch_correction_run_settings.tsv")))))
    stop("Existing batch outputs: use a fresh output directory; filename-only reuse is not allowed")
}

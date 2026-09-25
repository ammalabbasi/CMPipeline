#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble); library(vegan); library(zCompositions)})
argv <- commandArgs(TRUE); root <- argv[1]; out <- argv[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts/batch_contract.R"))
source(file.path(root, "scripts/batch_reporting.R"))
code <- parse(file.path(root, "scripts/2500703_batch_correction_normalization.r"))
for (expr in code) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr)
}
passed <- 0L
test <- function(name, expr) {
  force(expr); passed <<- passed + 1L; cat("PASS", name, "\n")
}
fails <- function(expr, pattern) {
  message <- tryCatch({force(expr); NA_character_}, error = function(e) conditionMessage(e))
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE))
    stop("Expected error '", pattern, "', got: ", message)
  invisible(TRUE)
}
args <- list(taxon_column = "clade_name", sample_id_column = "patient", batch_column = "batch",
             type_column = "type", tumor_only = FALSE, tumor_value = "Tumor", prefix = file.path(out, "case", "fixture"))
dir.create(dirname(args$prefix), recursive = TRUE)
taxa <- c("k__Bacteria|g__A|s__A_one", "k__Bacteria|g__B|s__B_one")
otu <- data.frame(clade_name = taxa, S1 = c(2, 4), S2 = c(6, 8), check.names = FALSE)
meta <- data.frame(patient = c("S1", "S2"), type = c("Tumor", "Normal"), batch = c("b1", "b2"))
o <- file.path(out, "otu.csv"); m <- file.path(out, "meta.tsv")
write.csv(otu, o, row.names = FALSE); write_tsv(meta, m)
loaded <- load_and_prepare_data(o, m)
test("CSV counts + TSV metadata", stopifnot(identical(rownames(loaded$otu), c("S1", "S2"))))
o2 <- file.path(out, "otu.txt"); m2 <- file.path(out, "meta.csv")
write_tsv(otu, o2); write.csv(meta, m2, row.names = FALSE)
test("TSV counts + quoted CSV metadata", stopifnot(identical(load_and_prepare_data(o2, m2)$otu, loaded$otu)))
write_tsv(meta[1, ], m)
test("missing metadata fatal", fails(load_and_prepare_data(o, m), "Metadata missing"))
z <- otu; z$S2 <- 0; write.csv(z, o, row.names = FALSE)
test("missing metadata for zero sample fatal", fails(load_and_prepare_data(o, m), "Metadata missing"))
write_tsv(rbind(meta, meta[1, ]), m)
test("duplicate metadata fatal", fails(load_and_prepare_data(o, m), "duplicate metadata sample IDs"))
write_tsv(meta, m)
for (bad in c("-1", "1.2", "Inf", "NA", "oops", "")) {
  z <- otu; z$S1 <- as.character(z$S1); z$S1[1] <- bad; write.csv(z, o, row.names = FALSE)
  test(paste("invalid count", bad), fails(load_and_prepare_data(o, m), if (bad == "-1") "nonnegative" else if (bad == "1.2") "integers" else "finite"))
}
writeLines(c('clade_name,S1,S1', 'a,1,2'), o)
test("duplicate count header fatal", fails(load_and_prepare_data(o, m), "column names"))
z <- otu; z$clade_name <- c("a", "a"); write.csv(z, o, row.names = FALSE)
test("duplicate taxon fatal", fails(load_and_prepare_data(o, m), "taxon IDs"))
write.csv(otu, o, row.names = FALSE); args$taxon_column <- "absent"
test("declared taxon cannot fall back", fails(load_and_prepare_data(o, m), "Declared taxon column"))
args$taxon_column <- "clade_name"
z <- otu; names(z)[2] <- "S1.bracken.G.krakenreport.txt"; write.csv(z, o, row.names = FALSE)
test("canonical ID suffix not stripped", fails(load_and_prepare_data(o, m), "Metadata missing"))
write.csv(otu[1, ], o, row.names = FALSE)
test("single taxon preserves dimensions", stopifnot(identical(dim(load_and_prepare_data(o, m)$otu), c(2L, 1L))))
write.csv(otu, o, row.names = FALSE); extra <- rbind(meta, data.frame(patient = "UNUSED", type = "Normal", batch = "b1")); write_tsv(extra, m)
test("extra metadata explicitly unused", stopifnot(identical(load_and_prepare_data(o, m)$unused_metadata, "UNUSED")))

# Count-independent CLR distance: the first coordinate has a negative column sum.
x <- matrix(c(-2, 1, 1, -1, -1, 2, -3, 2, 1), nrow = 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("t1", "t2", "t3")))
md <- data.frame(batch = c("b1", "b1", "b2"), row.names = rownames(x))
test("CLR distance retains nonpositive-sum coordinates", stopifnot(isTRUE(all.equal(as.numeric(batch_distance(x, md, "euclidean", TRUE)), as.numeric(dist(x)), tolerance = 1e-12))))
test("sample order mismatch rejected", fails(batch_distance(x, md[3:1, , drop = FALSE], "euclidean", TRUE), "sample order"))

# Cohort selection is shared and every exclusion is recorded.
ids <- sprintf("T%02d", 1:24)
counts <- as.data.frame(matrix(1:72, 24, 3, dimnames = list(ids, c("a", "b", "c"))))
metadata <- data.frame(patient = ids, batch = rep(c("b1", "b2"), each = 12),
                       type = rep(c("Tumor", "Normal"), 12), row.names = ids)
counts[1, ] <- 0; metadata$type[2] <- NA
d <- list(otu = counts, meta = metadata, unused_metadata = "Unused")
p <- prepare_batch_cohort(d, args, "type")
test("one ordered analysis cohort", stopifnot(nrow(p$otu) == 22L, identical(rownames(p$otu), rownames(p$meta))))
ledger <- read_tsv(file.path(dirname(args$prefix), "batch_correction_sample_status.tsv"), show_col_types = FALSE)
test("exclusions and unused rows recorded", stopifnot(sum(ledger$status == "excluded") == 2L, sum(ledger$status == "unused") == 1L))
constant <- d; constant$meta$type[13:24] <- "Tumor"
test("requested covariate not silently removed", fails(prepare_batch_cohort(constant, args, "type"), "does not vary"))
test("existing output destination rejected", fails(reject_existing_batch_outputs(dirname(args$prefix)), "Existing batch outputs"))

mat <- matrix(c(1, 3, 0, 5, 2, 1, 3, 0, 2, 4, 1, 2), nrow = 4, byrow = TRUE,
              dimnames = list(paste0("N", 1:4), paste0("f", 1:3)))
write_status <- function(status, reason) cat("STATUS", status, reason, "\n")
norm <- normalize_all(mat, file.path(out, "norm", "sample"))
test("all required transforms returned", stopifnot(identical(names(norm), c("CLR_pseudo", "CLR_CZM", "RCLR"))))
pseudo <- max(as.numeric(quantile(mat[mat > 0], .05)) * .5, .1)
reference <- log(mat + pseudo) - rowMeans(log(mat + pseudo))
test("CLR numeric oracle", stopifnot(max(abs(as.matrix(norm$CLR_pseudo) - reference)) < 1e-12))
test("normalization filenames never silently reused", fails(normalize_all(mat, file.path(out, "norm", "sample")), "Existing normalization"))
positive_mat <- mat + 1
positive_norm <- normalize_all(positive_mat, file.path(out, "zero_free", "sample"))
positive_clr <- log(positive_mat) - rowMeans(log(positive_mat))
test("zero-free input returns all required transforms", stopifnot(identical(names(positive_norm), c("CLR_pseudo", "CLR_CZM", "RCLR"))))
test("zero-free CZM equals direct CLR", stopifnot(max(abs(as.matrix(positive_norm$CLR_CZM) - positive_clr)) < 1e-12))
test("zero-free RCLR equals direct CLR", stopifnot(max(abs(as.matrix(positive_norm$RCLR) - positive_clr)) < 1e-12))
cmultRepl <- function(...) stop("injected CZM failure")
test("one method failure is fatal", fails(normalize_all(mat, file.path(out, "fault", "sample")), "Required normalization CLR_CZM failed"))
test("partial method outputs not written", stopifnot(!file.exists(file.path(out, "fault", "sample_clr_pseudo.tsv"))))
test("quality preserves negative final R2 but hides percentage", {
  q <- quality_summary(.9, -.01)
  stopifnot(q$final_R2 == -.01, is.na(q$reported_reduction_percent), grepl("negative_final", q$quality_status))
})
test("ordinary quality percentage", stopifnot(abs(report_reduction_percent(.4, .3) - 25) < 1e-10))
test("zero baseline is not interpretable", stopifnot(is.na(report_reduction_percent(0, 0))))
test("nonfinite metric not evaluable", stopifnot(quality_summary(NA_real_, .1)$quality_status == "not_evaluable_nonfinite_R2"))
test("worse fit not called improved", stopifnot(quality_summary(.2, .3)$quality_status == "no_reduction_in_batch_R2"))
test("seed validation", fails(configure_batch_rng(-1, 2), "seed must"))
test("thread validation", fails(configure_batch_rng(1, 0), "threads must"))
test("master RNG reproducible", {
  configure_batch_rng(17, 2); first <- runif(10)
  configure_batch_rng(17, 2); stopifnot(identical(first, runif(10)), RNGkind()[1] == "L'Ecuyer-CMRG")
})
test("unobserved fit outcomes explicitly labelled", stopifnot(quality_summary(.4, .3)$taxon_fit_diagnostics == "not_exposed_by_pinned_dependency"))
cat("SUMMARY", passed, "tests passed\n")

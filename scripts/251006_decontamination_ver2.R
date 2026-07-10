#!/usr/bin/env Rscript

# ==============================================================================
# SECTION 1: Load Libraries and Define Helper Functions
# ==============================================================================

# Load necessary libraries

# Load necessary libraries
if (!require("pacman")) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse, phyloseq, decontam, argparse, tools, ggrepel)

theme_set(theme_minimal(base_size = 12) + 
          theme(plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)))

# --- Helper Function 1: Data Loading and Preprocessing ---

load_and_prepare_data <- function(otu_path, meta_path, output_prefix = NULL, batch_column = "shipment_batch") {
  cat("INFO: Loading and preparing data...\n")
  
  otu_data_raw <- read.csv(otu_path, sep="\t", check.names = FALSE)
  metadata_raw <- read.csv(meta_path, sep="\t", check.names = FALSE)

  # 빈 컬럼명 처리 (Empty column names)
  col_names <- colnames(metadata_raw)
  empty_cols <- which(col_names == "" | is.na(col_names))
  if(length(empty_cols) > 0) {
    cat("WARN: Found", length(empty_cols), "empty column name(s). Removing...\n")
    metadata_raw <- metadata_raw[, -empty_cols, drop = FALSE]
  }

  # Sample ID column 자동 감지
  sample_id_candidates <- c("donor_id", "sample_id", "sampleid", "Sample.ID", "Sample_ID",
                            "SampleID", "Sample ID", "SAMPLE_ID", "Sample ID_x",
                            "subject_id", "SubjectID", "Subject_ID", "patient", "Patient")
  sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(metadata_raw)][1]

  if(is.na(sample_id_col)) {
    stop("ERROR: Sample ID column not found in metadata. Tried: ",
         paste(sample_id_candidates, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(metadata_raw), collapse = ", "))
  }

  cat("INFO: Using '", sample_id_col, "' as sample ID column.\n", sep = "")

  # sample_id로 컬럼명 표준화
  if(sample_id_col != "sample_id") {
    metadata_raw <- metadata_raw %>% rename(sample_id = !!sym(sample_id_col))
  }

  # Sample Type 처리: type_column 파라미터 사용 (없으면 자동 감지)
  type_column <- getOption("decontam_type_column", "cohort")

  if(type_column %in% colnames(metadata_raw)) {
    cat("INFO: Using '", type_column, "' column for Type information\n", sep = "")
    metadata_raw$Type <- metadata_raw[[type_column]]
  } else {
    # Fallback: 자동 감지
    type_candidates <- c("Type", "type", "Sample Type_x", "Sample_Type", "sample_type",
                         "cohort", "Cohort", "group", "Group", "condition", "Condition")
    type_col <- type_candidates[type_candidates %in% colnames(metadata_raw)][1]

    if(!is.na(type_col)) {
      cat("INFO: Specified type_column '", type_column, "' not found. Using '", type_col, "' instead.\n", sep = "")
      metadata_raw$Type <- metadata_raw[[type_col]]
    } else {
      cat("WARN: No Type column found in metadata. Creating default Type='Sample'\n")
      metadata_raw$Type <- "Sample"
    }
  }

  # OTU table formatting - Taxon name column 자동 감지
  taxon_name_candidates <- c("name", "species_name", "X", "species", "taxon",
                             "taxa", "clade_name", "Taxon", "Species",
                             "organism", "OTU", "ASV", "#OTU ID", "OTUID")
  name_col <- taxon_name_candidates[taxon_name_candidates %in% colnames(otu_data_raw)][1]

  # 후보에서 찾지 못했으면 첫 번째 non-numeric 컬럼 사용
  if (is.na(name_col)) {
    cat("WARN: Standard taxon name column not found. Available OTU table columns:\n")
    cat("  ", paste(head(colnames(otu_data_raw), 10), collapse = ", "), "...\n")

    # 첫 번째 문자열 컬럼 찾기
    non_numeric_cols <- which(!sapply(otu_data_raw, is.numeric))
    if (length(non_numeric_cols) > 0) {
      name_col <- colnames(otu_data_raw)[non_numeric_cols[1]]
      cat("INFO: Using first non-numeric column '", name_col, "' as taxon name.\n", sep = "")
    } else {
      stop("ERROR: Cannot identify taxon name column in OTU table.\n",
           "Available columns: ", paste(colnames(otu_data_raw), collapse = ", "))
    }
  } else {
    cat("INFO: Using '", name_col, "' as taxon name column.\n", sep = "")
  }

  rownames(otu_data_raw) <- otu_data_raw[[name_col]]
  
  numeric_cols <- sapply(otu_data_raw, is.numeric)
  otu_data <- as.data.frame(otu_data_raw[, numeric_cols])
  otu_data[is.na(otu_data)] <- 0
  
  # Sample ID normalization
  cat("\nINFO: Normalizing sample IDs...\n")
  original_sample_ids <- colnames(otu_data)

  # Remove common suffixes from OTU sample IDs:
  # 1. .bracken.G.krakenreport.txt or .bracken.S.krakenreport.txt
  # 2. Other potential suffixes
  normalized_sample_ids <- gsub("\\.bracken\\.[GS]\\.krakenreport\\.txt$", "", original_sample_ids)
  normalized_sample_ids <- gsub("\\.bracken\\.[GS]\\.mpa\\.krakenreport\\.txt$", "", normalized_sample_ids)
  
  cat("INFO: Example OTU sample IDs (before):", paste(head(original_sample_ids, 3), collapse = ", "), "\n")
  
  # 중복 처리
  if (any(duplicated(normalized_sample_ids))) {
    cat("WARN: Duplicate sample IDs detected after normalization.\n")
    keep_first <- !duplicated(normalized_sample_ids)
    otu_data <- otu_data[, keep_first]
    normalized_sample_ids <- normalized_sample_ids[keep_first]
    cat("INFO: Kept first occurrence of duplicates. Remaining samples:", ncol(otu_data), "\n")
  }
  
  colnames(otu_data) <- normalized_sample_ids
  cat("INFO: Example OTU sample IDs (after):", paste(head(colnames(otu_data), 3), collapse = ", "), "\n")
  cat("INFO: Example metadata sample IDs:", paste(head(metadata_raw$sample_id, 3), collapse = ", "), "\n")
  
  # Metadata 확장
  metadata_expanded <- create_or_expand_metadata(otu_data, metadata_raw, batch_column)
  
  # QC
  original_sample_count <- ncol(otu_data)
  otu_data <- otu_data[, (colSums(otu_data, na.rm=TRUE) != 0)]
  cat("INFO: Pruned", original_sample_count - ncol(otu_data), "empty samples.\n")

  original_taxa_count <- nrow(otu_data)
  taxa_sums <- rowSums(otu_data)
  otu_data <- otu_data[taxa_sums > 1, ]
  cat("INFO: Pruned", original_taxa_count - nrow(otu_data), "rare taxa (singletons).\n")

  # Intersection
  common_samples <- intersect(colnames(otu_data), metadata_expanded$sample_id)
  
  if (length(common_samples) == 0) {
    stop("ERROR: No common samples found.")
  }
  
  cat("INFO: Found", length(common_samples), "common samples out of", 
      ncol(otu_data), "OTU samples and", nrow(metadata_expanded), "metadata samples.\n")
  
  otu_data <- otu_data[, common_samples]
  metadata <- metadata_expanded[metadata_expanded$sample_id %in% common_samples, ]
  rownames(metadata) <- metadata$sample_id
  metadata <- metadata[colnames(otu_data), ]
  
  cat("INFO:", length(common_samples), "common samples prepared for final analysis.\n")
  
  return(list(otu = otu_data, meta = metadata, name_col = name_col))
}

# --- Helper Function: Extract Sample Type from TCGA ID ---
extract_sample_type_from_id <- function(sample_ids) {
  # FLEXIBLE: 두 가지 형식 지원
  # 1. 표준 TCGA 형식: TCGA-XX-XXXX-01 (숫자 코드)
  # 2. Suffix 형식: TCGA-XX-XXXX_TUMOR (텍스트 suffix)

  # 먼저 텍스트 suffix 확인 (_TUMOR, _NORMAL, etc.)
  text_suffix <- toupper(str_extract(sample_ids, "_(TUMOR|NORMAL|BLOOD|CONTROL)$"))
  text_suffix <- gsub("_", "", text_suffix)

  # 숫자 코드 확인 (표준 TCGA)
  type_codes <- str_extract(sample_ids, "-[0-9]{2}$")
  type_codes <- gsub("-", "", type_codes)

  # Type 결정 (text suffix 우선, 없으면 숫자 코드 사용)
  sample_types <- case_when(
    # Text suffix 기반 (우선순위)
    !is.na(text_suffix) & text_suffix == "TUMOR" ~ "Primary Tumor",
    !is.na(text_suffix) & text_suffix == "NORMAL" ~ "Solid Tissue Normal",
    !is.na(text_suffix) & text_suffix == "BLOOD" ~ "Blood Derived Normal",
    !is.na(text_suffix) & text_suffix == "CONTROL" ~ "Blood Derived Normal",

    # 숫자 코드 기반 (fallback)
    !is.na(type_codes) & type_codes == "01" ~ "Primary Tumor",
    !is.na(type_codes) & type_codes == "02" ~ "Recurrent Tumor",
    !is.na(type_codes) & type_codes == "03" ~ "Primary Blood Derived Cancer - Peripheral Blood",
    !is.na(type_codes) & type_codes == "10" ~ "Blood Derived Normal",
    !is.na(type_codes) & type_codes == "11" ~ "Solid Tissue Normal",
    !is.na(type_codes) & type_codes == "06" ~ "Metastatic",

    TRUE ~ "Unknown"
  )

  return(sample_types)
}

create_or_expand_metadata <- function(otu_df, meta_df, batch_column = "shipment_batch") {
  cat("\n========================================\n")
  cat("METADATA EXPANSION\n")
  cat("========================================\n")

  otu_samples <- colnames(otu_df)
  meta_samples <- meta_df$sample_id

  cat("OTU samples: ", length(otu_samples), "\n", sep = "")
  cat("Metadata samples: ", length(meta_samples), "\n", sep = "")

  # OTU에만 있는 샘플들
  otu_only <- setdiff(otu_samples, meta_samples)
  cat("Samples in OTU but not in metadata: ", length(otu_only), "\n", sep = "")
  
  if(length(otu_only) > 0) {
    cat("\nExtracting Type from OTU-only sample IDs...\n")
    
    # Type 추출
    otu_only_types <- extract_sample_type_from_id(otu_only)
    
    cat("Type distribution of new samples:\n")
    print(table(otu_only_types))
    
    # 환자 ID 추출 (FLEXIBLE: - 또는 _ 구분자 지원)
    extract_patient_id <- function(sample_ids) {
      # TCGA-XX-XXXX 형식 추출 (뒤에 -## 또는 _TUMOR 등이 올 수 있음)
      patient_ids <- sub("(TCGA-[^-]+-[^-]+)[-_].*", "\\1", sample_ids)

      # DEBUG: 추출이 제대로 되었는지 확인
      unchanged <- sum(patient_ids == sample_ids)
      if(unchanged > 0) {
        cat("DEBUG: extract_patient_id - ", unchanged, " IDs unchanged (regex not matched)\n", sep="")
        cat("DEBUG: Example unchanged ID: ", sample_ids[patient_ids == sample_ids][1], "\n", sep="")
      }

      return(patient_ids)
    }

    new_patient_ids <- extract_patient_id(otu_only)

    cat("DEBUG: Unique patients in new samples:", length(unique(new_patient_ids)), "\n")
    cat("DEBUG: First 5 patient IDs:", paste(head(unique(new_patient_ids), 5), collapse = ", "), "\n")
    
    # 기존 metadata에서 환자 ID 추출
    meta_df$patient_id <- extract_patient_id(meta_df$sample_id)
    
    cat("DEBUG: Unique patients in original metadata:", length(unique(meta_df$patient_id)), "\n")
    cat("DEBUG: First 5 original patient IDs:", paste(head(unique(meta_df$patient_id), 5), collapse = ", "), "\n")
    
    # 환자별 batch 매핑
    patient_batch_map <- meta_df %>%
      group_by(patient_id) %>%
      summarise(
        !!batch_column := first(na.omit(.data[[batch_column]])),
        n_samples = n(),
        .groups = 'drop'
      )

    cat("DEBUG: Patient-batch mapping created:", nrow(patient_batch_map), "unique patients\n")
    cat("DEBUG: Batch distribution in original metadata:\n")
    print(table(meta_df[[batch_column]]))
    
    # 새 샘플들과 매핑
    new_rows <- data.frame(
      sample_id = otu_only,
      patient_id = new_patient_ids,
      Type = otu_only_types,
      stringsAsFactors = FALSE
    )
    
    # Batch 매핑
    new_rows <- new_rows %>%
      left_join(patient_batch_map %>% select(patient_id, !!batch_column), by = "patient_id")

    # 매핑 결과 확인
    cat("DEBUG: Mapping results:\n")
    cat("  - Mapped successfully:", sum(!is.na(new_rows[[batch_column]])), "\n")
    cat("  - Failed to map:", sum(is.na(new_rows[[batch_column]])), "\n")

    unmapped <- sum(is.na(new_rows[[batch_column]]))
    if(unmapped > 0) {
      cat("WARN:", unmapped, "samples have no matching patient in original metadata\n")
      cat("     These will be assigned to 'Unknown' batch\n")
      unmapped_patients <- unique(new_rows$patient_id[is.na(new_rows[[batch_column]])])
      cat("     Unmapped patients (first 5):", paste(head(unmapped_patients, 5), collapse = ", "), "\n")
      new_rows[[batch_column]][is.na(new_rows[[batch_column]])] <- "Unknown"
    }

    cat("\nBatch distribution of NEW samples:\n")
    print(table(new_rows[[batch_column]]))
    
    # patient_id 제거
    new_rows$patient_id <- NULL
    meta_df$patient_id <- NULL

    # 컬럼 동기화 (RCC 스타일 + Type 보존)
    # meta_df에 있는 추가 컬럼들을 new_rows에 NA로 추가
    other_cols <- setdiff(colnames(meta_df), colnames(new_rows))
    for(col in other_cols) {
      new_rows[[col]] <- NA
    }

    # 컬럼 순서 맞추기 (이제 meta_df의 모든 컬럼이 new_rows에도 있음)
    new_rows <- new_rows[, colnames(meta_df)]

    # 중복 체크 및 제거
    duplicated_samples <- new_rows$sample_id[new_rows$sample_id %in% meta_df$sample_id]
    if(length(duplicated_samples) > 0) {
      cat("WARN: Found", length(duplicated_samples), "sample_id(s) that already exist in metadata.\n")
      cat("      These will NOT be added (keeping original metadata entries).\n")
      cat("      Duplicated IDs (first 10):", paste(head(duplicated_samples, 10), collapse = ", "), "\n")
      new_rows <- new_rows[!new_rows$sample_id %in% meta_df$sample_id, ]
    }

    # 병합
    meta_expanded <- rbind(meta_df, new_rows)

    # 중복 확인 후 rownames 설정
    if(any(duplicated(meta_expanded$sample_id))) {
      cat("ERROR: Duplicate sample_id detected after merge!\n")
      dups <- meta_expanded$sample_id[duplicated(meta_expanded$sample_id)]
      cat("Duplicates:", paste(head(dups, 10), collapse = ", "), "\n")
      # 중복 제거 (첫 번째 것만 유지)
      meta_expanded <- meta_expanded[!duplicated(meta_expanded$sample_id), ]
    }
    rownames(meta_expanded) <- meta_expanded$sample_id

    cat("\nExpanded metadata summary:\n")
    cat("  Original: ", nrow(meta_df), " samples\n", sep = "")
    cat("  Added:    ", nrow(new_rows), " samples\n", sep = "")
    cat("  Total:    ", nrow(meta_expanded), " samples\n", sep = "")

    cat("\nBatch distribution AFTER expansion:\n")
    print(table(meta_expanded[[batch_column]]))
    
  } else {
    cat("\nNo expansion needed.\n")
    meta_expanded <- meta_df
  }
  
  # ===== Type 표준화 (파라미터 기반) =====
  cat("\n--- Type Standardization ---\n")
  cat("Before standardization:\n")
  print(table(meta_expanded$Type, useNA = "ifany"))

  # 파라미터에서 Tumor/Control 값 목록 가져오기
  tumor_values <- getOption("decontam_tumor_values", c("Cancer", "Precancer", "Tumor", "Primary Tumor", "Recurrent Tumor", "Metastatic"))
  control_values <- getOption("decontam_control_values", c("Normal", "Healthy", "Control", "Blood", "Solid Tissue Normal"))

  cat("INFO: Tumor values:", paste(tumor_values, collapse = ", "), "\n")
  cat("INFO: Control values:", paste(control_values, collapse = ", "), "\n")

  # Tumor 매핑
  for(val in tumor_values) {
    meta_expanded$Type[meta_expanded$Type == val] <- "Tumor"
  }

  # Control 매핑
  for(val in control_values) {
    meta_expanded$Type[meta_expanded$Type == val] <- "Control"
  }

  # Blood 패턴 매칭 (추가 안전장치)
  meta_expanded$Type[grepl("Blood", meta_expanded$Type, ignore.case = TRUE)] <- "Control"

  # Unknown 처리
  remaining_types <- unique(meta_expanded$Type[!meta_expanded$Type %in% c("Tumor", "Control")])
  if(length(remaining_types) > 0) {
    cat("WARN: Found unrecognized types:", paste(remaining_types, collapse = ", "), "\n")
    cat("      These will remain as-is. Consider adding them to --tumor_values or --control_values\n")
  }

  cat("\nAfter standardization:\n")
  print(table(meta_expanded$Type, useNA = "ifany"))
  
  cat("========================================\n\n")
  
  return(meta_expanded)
}



# --- Helper Function 2: Batch-specific Scatter Plot ---

create_and_save_plot <- function(plot_data, blacklist, title, plot_path) {
  if (!"is_contaminant" %in% colnames(plot_data)) {
    cat("WARN: 'is_contaminant' column not found for", title, ". Skipping plot.\n")
    return(invisible(NULL)) 
  }

  blacklist_flexible <- gsub(" ", "[ _]", blacklist)
  blacklist_pattern <- paste(blacklist_flexible, collapse = "|")

  # Species 이름 짧게 처리 (마지막 부분만 추출)
  # 예: "k__Bacteria|...|s__Streptococcus_parasanguinis" -> "Streptococcus parasanguinis"
  plot_data <- plot_data %>%
    mutate(
      short_name = case_when(
        grepl("\\|s__", name) ~ gsub(".*\\|s__", "", name),  # s__ 이후 추출
        grepl("\\|g__", name) ~ gsub(".*\\|g__", "", name),  # genus level
        TRUE ~ name  # 기타
      ),
      short_name = gsub("_", " ", short_name),  # underscore를 space로
      label = ifelse(grepl(blacklist_pattern, name, ignore.case = TRUE), short_name, "")
    )

  n_contam <- sum(plot_data$is_contaminant)
  
  p <- ggplot(plot_data, aes(x = prevalence_control, y = prevalence_sample)) +
    geom_point(data = plot_data %>% filter(!is_contaminant), 
               color = "gray50", alpha = 0.5, size = 2) +
    geom_point(data = plot_data %>% filter(is_contaminant), 
               color = "red", alpha = 0.7, size = 3) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                    box.padding = 0.5, segment.color = 'grey50') +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", size = 0.8) +
    annotate("text", x = 0.7, y = 0.9, label = "Tumor-enriched", 
             color = "blue", size = 3.5, fontface = "italic") +
    annotate("text", x = 0.9, y = 0.7, label = "Blood-enriched\n(Contaminants)", 
             color = "red", size = 3.5, fontface = "italic") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = title,
      subtitle = paste(n_contam, "contaminants identified (red points)"),
      x = "Prevalence in Blood (Control)",
      y = "Prevalence in Tumor (Sample)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "red", size = 10),
      panel.grid.minor = element_blank()
    )
  
  ggsave(plot_path, plot = p, width = 9, height = 8, dpi = 300, bg = "white")
  cat("INFO: Saved plot to", plot_path, "\n")
  
  return(p)
}

# --- Helper Function 3: Decontam Package Identification ---

identify_contaminants_decontam <- function(otu_df, meta_df, blacklist_species, plot_dir_path, threshold, batch_column = "shipment_batch") {
  cat("INFO: Identifying contaminants using decontam package...\n")

  # ===== 배치별 샘플 분포 확인 (수정) =====
  cat("\n========================================\n")
  cat("BATCH DISTRIBUTION ANALYSIS\n")
  cat("========================================\n")

  # Solid Tissue Normal 제외하고 계산
  batch_summary <- meta_df %>%
    filter(Type %in% c("Tumor", "Control")) %>%  # ✅ 이것만 포함
    group_by(!!sym(batch_column), Type) %>%
    summarise(n_samples = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Type, values_from = n_samples, values_fill = 0) %>%
    mutate(
      Tumor = if("Tumor" %in% names(.)) Tumor else 0,
      Control = if("Control" %in% names(.)) Control else 0,
      total = Tumor + Control,
      usable = (Control > 0 & Tumor > 0 & total > 2)
    )

  # Solid Tissue Normal 정보 추가 (참고용)
  solid_normal_count <- meta_df %>%
    filter(Type == "Solid Tissue Normal") %>%
    group_by(!!sym(batch_column)) %>%
    summarise(n_solid = n(), .groups = 'drop')

  batch_summary <- batch_summary %>%
    left_join(solid_normal_count, by = batch_column) %>%
    mutate(n_solid = ifelse(is.na(n_solid), 0, n_solid)) %>%
    select(!!sym(batch_column), Tumor, Control, n_solid, total, usable)
  
  cat("Batch-wise sample distribution:\n")
  cat("(Solid Tissue Normal shown separately - not used in decontam)\n\n")
  print(batch_summary)
  
  suitable_batches <- batch_summary %>% filter(usable == TRUE)
  cat("\n", nrow(suitable_batches), " batches suitable for decontam (out of ", nrow(batch_summary), ")\n", sep = "")
  
  if(nrow(suitable_batches) == 0) {
    cat("\nWARNING: NO batches meet criteria!\n")
    cat("All batches will be skipped. Consider using overall cohort approach.\n")
  }
  cat("========================================\n\n")
  
  master_table <- data.frame(Contaminant = character(0), Batch_Count = character(0))

  for(batch_id in unique(meta_df[[batch_column]])) {
    meta_batch <- meta_df[meta_df[[batch_column]] == batch_id, ]
    
    if(nrow(meta_batch) < 3) {
      cat("SKIP: Batch '", batch_id, "' - only ", nrow(meta_batch), " samples\n", sep = "")
      next
    }
    
    # is.neg 설정: Control로 표준화된 Type 사용
    meta_batch$is.neg <- (meta_batch$Type == "Control")
    
    n_neg <- sum(meta_batch$is.neg)
    n_pos <- sum(meta_batch$Type == "Tumor")
    n_solid_normal <- sum(meta_batch$Type == "Solid Tissue Normal")
    
    if(n_neg == 0) {
      cat("SKIP: Batch '", batch_id, "' - no control samples\n", sep = "")
      next
    }
    
    if(n_pos == 0) {
      cat("SKIP: Batch '", batch_id, "' - no tumor samples\n", sep = "")
      next
    }
    
    cat("PROCESS: Batch '", batch_id, "' - ", 
    n_pos, " tumor + ", n_neg, " blood control",
    if(n_solid_normal > 0) paste0(" (+ ", n_solid_normal, " solid normal, not used in decontam)") else "",
    "\n", sep = "")
    
    otu_batch <- otu_df[, meta_batch$sample_id]
    physeq_batch <- phyloseq(otu_table(as.matrix(otu_batch), taxa_are_rows = TRUE), sample_data(meta_batch))
    contam_results <- isContaminant(physeq_batch, method="prevalence", neg="is.neg", threshold=threshold, detailed = TRUE)
    
    # 플롯 생성
    ps_pa <- transform_sample_counts(physeq_batch, function(abund) 1*(abund>0))
    ps_pa_neg <- prune_samples(sample_data(ps_pa)$is.neg, ps_pa)
    ps_pa_pos <- prune_samples(!sample_data(ps_pa)$is.neg, ps_pa)

    plot_data_decontam <- data.frame(
      name = taxa_names(ps_pa),
      prevalence_control = (taxa_sums(ps_pa_neg) / nsamples(ps_pa_neg)),
      prevalence_sample = (taxa_sums(ps_pa_pos) / nsamples(ps_pa_pos)),
      is_contaminant = contam_results$contaminant
    )

    plot_title <- paste("Decontam Package - Batch:", batch_id)
    safe_batch_name <- gsub("[^A-Za-z0-9_-]", "_", batch_id)
    plot_path <- file.path(plot_dir_path, paste0("decontam_pkg_", safe_batch_name, ".png"))
    create_and_save_plot(plot_data_decontam, blacklist_species, plot_title, plot_path)
    
    contaminants_in_batch <- rownames(contam_results)[contam_results$contaminant == TRUE]
    
    if(length(contaminants_in_batch) > 0){
      cat("  → Found ", length(contaminants_in_batch), " contaminants\n", sep = "")
      master_table <- rbind(master_table, data.frame(Contaminant = contaminants_in_batch, Batch_Count = batch_id))
    } else {
      cat("  → No contaminants identified\n")
    }
  }
  
  return(master_table)
}

# --- Helper Function 4: Filter Rare OTUs by Prevalence ---

filter_rare_otu <- function(otu_table, prevalence_threshold = 0.05, min_abundance = 10, plots_dir = NULL) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("RARE TAXA FILTERING (2-STEP PROCESS)\n")
  cat(strrep("=", 70), "\n")
  cat("  Step 1: Sample-wise abundance threshold:", min_abundance, "reads\n")
  cat("  Step 2: Cohort-wide prevalence threshold:", prevalence_threshold * 100, "%\n")
  cat(strrep("=", 70), "\n\n")

  otu_table_numeric <- otu_table %>% mutate(across(everything(), as.numeric))
  original_taxa <- nrow(otu_table_numeric)
  total_samples <- ncol(otu_table_numeric)
  
  # ===== STEP 1: Sample-wise abundance filtering =====
  cat("STEP 1: Applying sample-wise abundance filter...\n")
  cat("────────────────────────────────────────────────\n")
  
  # 각 샘플에서 min_abundance 미만을 0으로 처리
  otu_abundance_filtered <- otu_table_numeric
  otu_abundance_filtered[otu_abundance_filtered < min_abundance] <- 0
  
  # 통계
  n_values_zeroed <- sum(otu_table_numeric > 0 & otu_table_numeric < min_abundance)
  total_nonzero_values <- sum(otu_table_numeric > 0)
  pct_zeroed <- 100 * n_values_zeroed / total_nonzero_values
  
  cat(sprintf("  - Original non-zero values:       %8d\n", total_nonzero_values))
  cat(sprintf("  - Values < %d reads (zeroed):     %8d (%.2f%% of non-zero)\n", 
              min_abundance, n_values_zeroed, pct_zeroed))
  
  # Taxa가 완전히 사라졌는지 확인
  taxa_present_before <- rowSums(otu_table_numeric > 0) > 0
  taxa_present_after <- rowSums(otu_abundance_filtered > 0) > 0
  taxa_lost_to_abundance <- taxa_present_before & !taxa_present_after
  n_lost_abundance <- sum(taxa_lost_to_abundance)
  
  if(n_lost_abundance > 0) {
    cat(sprintf("  - Taxa lost completely:           %8d\n", n_lost_abundance))
    lost_taxa_names <- rownames(otu_table_numeric)[taxa_lost_to_abundance]
    cat("    Examples:", paste(head(lost_taxa_names, 3), collapse=", "), "\n")
  } else {
    cat("  - Taxa lost completely:           %8d\n", 0)
  }
  
  # ===== STEP 2: Prevalence filtering =====
  cat("\nSTEP 2: Applying prevalence filter...\n")
  cat("────────────────────────────────────────────────\n")
  
  # Step 1 후 데이터에서 prevalence 계산
  presence_counts_before <- rowSums(otu_table_numeric > 0)
  presence_counts_after <- rowSums(otu_abundance_filtered > 0)
  prevalence_pct_before <- 100 * presence_counts_before / total_samples
  prevalence_pct_after <- 100 * presence_counts_after / total_samples
  
  min_sample_count <- total_samples * prevalence_threshold
  passes_prevalence <- presence_counts_after >= min_sample_count
  
  # 통계 출력
  cat(sprintf("  - Total samples in cohort:        %8d\n", total_samples))
  cat(sprintf("  - Minimum samples required:       %8.0f (%.1f%%)\n",
              min_sample_count, prevalence_threshold * 100))
  cat(sprintf("  - Taxa passing prevalence:        %8d\n", sum(passes_prevalence)))
  cat(sprintf("  - Taxa failing prevalence:        %8d\n", sum(!passes_prevalence)))
  
  # ===== Final filtering =====
  filtered_otu_table <- otu_abundance_filtered[passes_prevalence, ]
  
  # ===== Summary =====
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("FILTERING SUMMARY\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("Original taxa:                    %8d (100.0%%)\n", original_taxa))
  cat(sprintf("Lost to abundance filter:         %8d (%5.1f%%)\n", 
              n_lost_abundance, 100 * n_lost_abundance / original_taxa))
  cat(sprintf("Lost to prevalence filter:        %8d (%5.1f%%)\n", 
              sum(!passes_prevalence) - n_lost_abundance,
              100 * (sum(!passes_prevalence) - n_lost_abundance) / original_taxa))
  cat(sprintf("Final retained taxa:              %8d (%5.1f%%)\n", 
              nrow(filtered_otu_table),
              100 * nrow(filtered_otu_table) / original_taxa))
  cat(strrep("=", 70), "\n\n")
  
  # ===== 시각화 =====
  if(!is.null(plots_dir)) {
    cat("Generating visualization...\n")
    
    comparison_data <- data.frame(
      taxa = rownames(otu_table_numeric),
      prevalence_before = prevalence_pct_before,
      prevalence_after = prevalence_pct_after,
      mean_abundance_before = rowMeans(otu_table_numeric),
      mean_abundance_after = rowMeans(otu_abundance_filtered),
      status = case_when(
        !passes_prevalence ~ "Filtered out",
        prevalence_pct_after < prevalence_pct_before ~ "Prevalence decreased",
        TRUE ~ "Unchanged"
      )
    )
    
    # 색상 및 크기 설정
    comparison_data$status <- factor(comparison_data$status, 
                                     levels = c("Unchanged", "Prevalence decreased", "Filtered out"))
    
    # Plot 1: Prevalence Before vs After (main plot)
    p1 <- ggplot(comparison_data, aes(x = prevalence_before, y = prevalence_after, color = status, size = status, alpha = status)) +
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed", size = 0.8) +
      geom_point() +
      scale_color_manual(values = c("Unchanged" = "gray50", 
                                     "Prevalence decreased" = "orange", 
                                     "Filtered out" = "red")) +
      scale_size_manual(values = c("Unchanged" = 1.5, 
                                    "Prevalence decreased" = 2, 
                                    "Filtered out" = 2.5)) +
      scale_alpha_manual(values = c("Unchanged" = 0.3, 
                                     "Prevalence decreased" = 0.6, 
                                     "Filtered out" = 0.8)) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      annotate("text", x = 80, y = 15, 
               label = sprintf("Prevalence decreased\n(n=%d)", sum(comparison_data$status == "Prevalence decreased")),
               color = "orange", size = 3.5, fontface = "bold") +
      annotate("text", x = 15, y = 2, 
               label = sprintf("Filtered out\n(n=%d)", sum(comparison_data$status == "Filtered out")),
               color = "red", size = 3.5, fontface = "bold") +
      labs(
        title = "Effect of Sample-wise Abundance Filtering on Taxa Prevalence",
        subtitle = sprintf("Abundance threshold: %d reads | Prevalence threshold: %.1f%%", 
                          min_abundance, prevalence_threshold * 100),
        x = "Prevalence Before Filtering (%)",
        y = "Prevalence After Abundance Filtering (%)",
        color = "Taxa Status",
        size = "Taxa Status",
        alpha = "Taxa Status"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    plot_path_1 <- file.path(plots_dir, "rare_filtering_prevalence_impact.png")
    ggsave(plot_path_1, plot = p1, width = 10, height = 8, dpi = 300, bg = "white")
    cat("  → Saved:", plot_path_1, "\n")
    
    # Plot 2: Distribution of prevalence changes
    comparison_data$prevalence_change <- comparison_data$prevalence_after - comparison_data$prevalence_before
    
    p2 <- ggplot(comparison_data, aes(x = prevalence_change, fill = status)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 0.8) +
      scale_fill_manual(values = c("Unchanged" = "gray50", 
                                    "Prevalence decreased" = "orange", 
                                    "Filtered out" = "red")) +
      labs(
        title = "Distribution of Prevalence Changes After Abundance Filtering",
        subtitle = sprintf("%d taxa analyzed", nrow(comparison_data)),
        x = "Change in Prevalence (% points)",
        y = "Number of Taxa",
        fill = "Taxa Status"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    plot_path_2 <- file.path(plots_dir, "rare_filtering_prevalence_distribution.png")
    ggsave(plot_path_2, plot = p2, width = 10, height = 6, dpi = 300, bg = "white")
    cat("  → Saved:", plot_path_2, "\n\n")
    
  } else {
    cat("INFO: No plots_dir provided, skipping visualization.\n\n")
  }

  return(filtered_otu_table)
}
                                     
# --- Helper Function 5: Filter Contaminants by Batch Frequency ---

filter_contaminants_by_frequency <- function(master_contaminant_table, min_batches = 2) {
  cat("\nINFO: Filtering contaminants that appear in at least", min_batches, "batches...\n")
  
  if(nrow(master_contaminant_table) == 0) {
    cat("INFO: No contaminants to filter.\n")
    return(data.frame(Contaminant = character(0), n_batches = integer(0), batch_list = character(0)))
  }
  
  # 배치별 집계
  contaminant_summary <- master_contaminant_table %>%
    group_by(Contaminant) %>%
    summarise(
      n_batches = n_distinct(Batch_Count),
      batch_list = paste(unique(Batch_Count), collapse = ", "),
      .groups = 'drop'
    )
  
  # 임계값 적용
  final_contaminants <- contaminant_summary %>%
    filter(n_batches >= min_batches)
  
  cat("INFO: Initial contaminants:", nrow(contaminant_summary), "\n")
  cat("INFO: After filtering (>=", min_batches, "batches):", nrow(final_contaminants), "\n")
  
  # 제거된 contaminants 로깅
  removed <- contaminant_summary %>% filter(n_batches < min_batches)
  if(nrow(removed) > 0) {
    cat("INFO: Removed", nrow(removed), "contaminants found in <", min_batches, "batches\n")
  }
  
  return(final_contaminants)
}

# --- Helper Function 6: Remove Contaminants Globally ---

remove_contaminants_global <- function(original_otu, final_contaminants_df) {
  clean_otu <- original_otu
  
  if (nrow(final_contaminants_df) > 0) {
    cat("\nINFO: Removing", nrow(final_contaminants_df), "contaminants globally from all samples...\n")

    to_remove <- final_contaminants_df$Contaminant[final_contaminants_df$Contaminant %in% rownames(clean_otu)]
    clean_otu[to_remove, ] <- 0
    for (name in to_remove) {
      n_batches <- final_contaminants_df$n_batches[final_contaminants_df$Contaminant == name]
      cat("  - Removed:", name, "(found in", n_batches, "batches)\n")
    }
  } else {
    cat("INFO: No contaminants to remove.\n")
  }
  
  return(clean_otu)
}

# --- Helper Function 7: Save Decontaminated Table ---

save_decontaminated_table <- function(clean_otu_table, original_name_col, file_prefix, method_suffix) {
  output_file_path <- paste0(file_prefix, ".", method_suffix, "_decontaminated.csv")
  cat("INFO: Writing", method_suffix, "decontaminated OTU table to:", output_file_path, "\n")
  clean_otu_table %>%
    rownames_to_column(var = original_name_col) %>%
    write.csv(., output_file_path, row.names = FALSE)
}

# --- Helper Function 8: Overall Cohort Scatter Plot ---

plot_overall_contaminants_scatter <- function(otu_clean, meta_df, 
                                               master_contaminant_table_raw,
                                               final_contaminants_df,
                                               blacklist, plot_path) {
  cat("\nINFO: Generating overall cohort scatter plot...\n")
  
  # 1. OTU를 presence/absence로 변환
  otu_pa <- otu_clean > 0
  
  # 2. 샘플을 Blood/Tumor로 분리
  meta_df <- meta_df %>%
    mutate(cohort = case_when(
        tolower(Type) == "tumor" ~ "tumor",
        tolower(Type) == "control" ~ "control",  # Blood Derived Normal
        TRUE ~ "other"  # Solid Tissue Normal 등
  ))
  
  control_samples <- meta_df %>% filter(cohort == "control") %>% pull(sample_id)
  tumor_samples <- meta_df %>% filter(cohort == "tumor") %>% pull(sample_id)  

  if(length(control_samples) == 0 || length(tumor_samples) == 0) {  # ✅
    cat("WARN: Missing control or tumor samples. Skipping overall scatter plot.\n")
    return(invisible(NULL))
  }  

  # 3. 전체 코호트에서 prevalence 계산
  prevalence_control <- rowSums(otu_pa[, control_samples, drop = FALSE]) / length(control_samples)  # ✅
  prevalence_tumor <- rowSums(otu_pa[, tumor_samples, drop = FALSE]) / length(tumor_samples)  

  # 4. 플롯 데이터 생성
  plot_data <- data.frame(
    name = names(prevalence_control),  # ✅
    prevalence_control = prevalence_control,  # ✅
    prevalence_sample = prevalence_tumor
  )
  
  # 5. Contaminant 상태 추가
  all_batch_contaminants <- unique(master_contaminant_table_raw$Contaminant)
  final_contaminants <- final_contaminants_df$Contaminant
  
  plot_data <- plot_data %>%
    mutate(
      status = case_when(
        name %in% final_contaminants ~ "final",
        name %in% all_batch_contaminants ~ "batch_only",
        TRUE ~ "clean"
      )
    )
  
  # 6. Blacklist 매칭 및 Species 이름 짧게 처리
  blacklist_flexible <- gsub(" ", "[ _]", blacklist)
  blacklist_pattern <- paste(blacklist_flexible, collapse = "|")

  # Species 이름 짧게 처리
  plot_data <- plot_data %>%
    mutate(
      short_name = case_when(
        grepl("\\|s__", name) ~ gsub(".*\\|s__", "", name),  # s__ 이후 추출
        grepl("\\|g__", name) ~ gsub(".*\\|g__", "", name),  # genus level
        TRUE ~ name  # 기타
      ),
      short_name = gsub("_", " ", short_name),  # underscore를 space로
      label = case_when(
        status == "final" ~ short_name,
        grepl(blacklist_pattern, name, ignore.case = TRUE) ~ short_name,
        TRUE ~ ""
      )
    )
  
  # 7. Scatter plot 생성
  n_final <- sum(plot_data$status == "final")
  n_batch_only <- sum(plot_data$status == "batch_only")
  
  p <- ggplot(plot_data, aes(x = prevalence_control, y = prevalence_sample)) +
    geom_point(data = plot_data %>% filter(status == "clean"), 
               color = "gray70", alpha = 0.3, size = 1.5) +
    geom_point(data = plot_data %>% filter(status == "batch_only"), 
               color = "orange", alpha = 0.6, size = 2.5) +
    geom_point(data = plot_data %>% filter(status == "final"), 
               color = "red", alpha = 0.8, size = 4) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 30,
                    box.padding = 0.5, segment.color = 'grey50',
                    min.segment.length = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", size = 0.8) +
    annotate("text", x = 0.75, y = 0.95, label = "Tumor-enriched", 
             color = "blue", size = 4, fontface = "italic") +
    annotate("text", x = 0.95, y = 0.75, label = "Blood-enriched\n(Contaminants)", 
             color = "red", size = 4, fontface = "italic") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "Overall Cohort: Contaminant Identification",
      subtitle = paste0(
        "Red (n=", n_final, "): Final contaminants to be removed (≥2 batches) | ",
        "Orange (n=", n_batch_only, "): Identified in only 1 batch"
      ),
      x = "Prevalence in Blood (All Batches Combined)",
      y = "Prevalence in Tumor (All Batches Combined)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),  
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  ggsave(plot_path, plot = p, width = 10, height = 9, dpi = 300)
  cat("INFO: Saved overall scatter plot to", plot_path, "\n")
  
  return(invisible(NULL))
}

# --- Helper Function 9: Contaminant Proportions Comparison (NEW) ---

plot_contaminant_proportions_comparison <- function(original_otu,
                                                     master_contaminant_table_raw,
                                                     final_contaminants_df,
                                                     rare_taxa_list,
                                                     meta_df,
                                                     plots_dir,
                                                     batch_column = "shipment_batch") {
  cat("\nINFO: Generating contaminant proportion comparison plots...\n")
  
  total_reads <- colSums(original_otu)
  
  # ===== Rare taxa reads 계산 =====
  rare_reads <- original_otu %>%
    rownames_to_column(var = "name") %>%
    filter(name %in% rare_taxa_list) %>%
    pivot_longer(cols = -name, names_to = "sample_id", values_to = "reads") %>%
    group_by(sample_id) %>%
    summarise(rare_reads = sum(reads, na.rm = TRUE), .groups = 'drop')
  
  # ===== Decontam contaminant reads 계산 =====
  contaminant_list <- final_contaminants_df$Contaminant
  decontam_reads <- original_otu %>%
    rownames_to_column(var = "name") %>%
    filter(name %in% contaminant_list) %>%
    pivot_longer(cols = -name, names_to = "sample_id", values_to = "reads") %>%
    group_by(sample_id) %>%
    summarise(decontam_reads = sum(reads, na.rm = TRUE), .groups = 'drop')
  
  # ===== 통합 데이터 =====
  summary_combined <- data.frame(
    sample_id = names(total_reads),
    total_reads = total_reads
  ) %>%
    inner_join(dplyr::select(meta_df, sample_id, !!sym(batch_column), Type), by = "sample_id") %>%
    left_join(rare_reads, by = "sample_id") %>%
    left_join(decontam_reads, by = "sample_id") %>%
    mutate(
      rare_reads = ifelse(is.na(rare_reads), 0, rare_reads),
      decontam_reads = ifelse(is.na(decontam_reads), 0, decontam_reads),
      total_contaminant_reads = rare_reads + decontam_reads,
      
      rare_proportion = ifelse(total_reads > 0, rare_reads / total_reads, 0),
      decontam_proportion = ifelse(total_reads > 0, decontam_reads / total_reads, 0),
      total_proportion = ifelse(total_reads > 0, total_contaminant_reads / total_reads, 0),
      
      Type = case_when(
        Type == "Tumor" ~ "Tumor",
        Type == "Control" ~ "Blood",
        Type == "Normal" ~ "Normal Tissue",
        Type == "Solid Tissue Normal" ~ "Normal Tissue",
        TRUE ~ Type
      ),
      Type = factor(Type, levels = c("Blood", "Tumor", "Normal Tissue"))
    )
  
# ===== Plot 1: Stacked Bar - Mean proportions =====
summary_means <- summary_combined %>%
  group_by(Type) %>%
  summarise(
    n = n(),
    mean_rare = mean(rare_proportion),
    mean_decontam = mean(decontam_proportion),
    se_rare = sd(rare_proportion) / sqrt(n()),
    se_decontam = sd(decontam_proportion) / sqrt(n()),
    .groups = 'drop'
  )

summary_means_long <- summary_means %>%
  select(Type, mean_rare, mean_decontam) %>%
  pivot_longer(cols = c(mean_rare, mean_decontam), 
               names_to = "contaminant_type", 
               values_to = "proportion") %>%
  mutate(contaminant_type = factor(contaminant_type,
                                    levels = c("mean_rare", "mean_decontam"),
                                    labels = c("Rare Taxa", "Decontam")))

p1 <- ggplot(summary_means_long, aes(x = Type, y = proportion, fill = contaminant_type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.6) +
  geom_text(data = summary_means, 
            aes(x = Type, y = mean_rare + mean_decontam, label = sprintf("%.1f%%", (mean_rare + mean_decontam) * 100)),
            vjust = -0.5, size = 5, fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("Rare Taxa" = "#FF7F00", "Decontam" = "#984EA3")) +
  labs(
    title = "Mean Contaminant Proportion by Sample Type",
    subtitle = sprintf("Rare taxa: %d (low prevalence) | Decontam: %d (blood-enriched)", 
                      length(rare_taxa_list), nrow(final_contaminants_df)),
    x = "Sample Type",
    y = "Mean Contaminant Fraction",
    fill = "Source"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

ggsave(file.path(plots_dir, "proportion_stacked_by_type.png"), 
       plot = p1, width = 10, height = 7, dpi = 300)

# ===== Plot 2: Stacked Bar by Batch - 각 Type별 facet =====
summary_batch_means <- summary_combined %>%
  group_by(Type, !!sym(batch_column)) %>%
  summarise(
    n = n(),
    mean_rare = mean(rare_proportion),
    mean_decontam = mean(decontam_proportion),
    .groups = 'drop'
  )

summary_batch_long <- summary_batch_means %>%
  pivot_longer(cols = c(mean_rare, mean_decontam),
               names_to = "contaminant_type",
               values_to = "proportion") %>%
  mutate(contaminant_type = factor(contaminant_type,
                                    levels = c("mean_rare", "mean_decontam"),
                                    labels = c("Rare Taxa", "Decontam")))

p2 <- ggplot(summary_batch_long, aes(x = !!sym(batch_column), y = proportion, fill = contaminant_type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 0.7) +
  facet_wrap(~Type, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Rare Taxa" = "#FF7F00", "Decontam" = "#984EA3")) +
  labs(
    title = "Mean Contaminant Proportion by Batch and Type",
    x = "Batch ID",
    y = "Mean Contaminant Fraction",
    fill = "Source"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(file.path(plots_dir, "proportion_stacked_by_batch.png"), 
       plot = p2, width = 14, height = 10, dpi = 300)

# ===== Plot 3: Individual samples as stacked bar (선택적) =====
# 샘플이 너무 많으면 skip
if(nrow(summary_combined) <= 100) {
  summary_sample_long <- summary_combined %>%
    arrange(Type, desc(total_proportion)) %>%
    mutate(sample_order = row_number()) %>%
    select(sample_order, sample_id, Type, rare_proportion, decontam_proportion) %>%
    pivot_longer(cols = c(rare_proportion, decontam_proportion),
                 names_to = "contaminant_type",
                 values_to = "proportion") %>%
    mutate(contaminant_type = factor(contaminant_type,
                                      levels = c("rare_proportion", "decontam_proportion"),
                                      labels = c("Rare Taxa", "Decontam")))
  
  p3 <- ggplot(summary_sample_long, aes(x = reorder(sample_id, sample_order), y = proportion, fill = contaminant_type)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    facet_wrap(~Type, scales = "free_x", ncol = 1) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Rare Taxa" = "#FF7F00", "Decontam" = "#984EA3")) +
    labs(
      title = "Contaminant Proportion per Sample",
      x = "Sample (ordered by total contaminant)",
      y = "Contaminant Fraction",
      fill = "Source"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ggsave(file.path(plots_dir, "proportion_stacked_per_sample.png"), 
         plot = p3, width = 14, height = 10, dpi = 300)
}
  
  # ===== Summary statistics =====
  cat("\n========================================\n")
  cat("CONTAMINANT PROPORTION SUMMARY\n")
  cat("========================================\n")
  
  summary_stats <- summary_combined %>%
    group_by(Type) %>%
    summarise(
      n = n(),
      mean_rare = mean(rare_proportion) * 100,
      mean_decontam = mean(decontam_proportion) * 100,
      mean_total = mean(total_proportion) * 100,
      median_total = median(total_proportion) * 100,
      .groups = 'drop'
    )
  
  cat(sprintf("%-15s %8s %12s %12s %12s %12s\n", 
              "Type", "N", "Rare(%)", "Decontam(%)", "Total(%)", "Median(%)"))
  cat(strrep("-", 75), "\n")
  for(i in 1:nrow(summary_stats)) {
    cat(sprintf("%-15s %8d %12.2f %12.2f %12.2f %12.2f\n",
                summary_stats$Type[i],
                summary_stats$n[i],
                summary_stats$mean_rare[i],
                summary_stats$mean_decontam[i],
                summary_stats$mean_total[i],
                summary_stats$median_total[i]))
  }
  cat("========================================\n\n")
  
  cat("INFO: All contaminant proportion plots completed.\n")
  return(invisible(NULL))
}

# --- Helper Function: Generate Summary Report ---
generate_summary_report <- function(otu_original, otu_pre_filtered, otu_clean, otu_filtered,
                                     meta_df, master_contaminant_table_raw,
                                     final_contaminants_df, output_prefix,
                                     min_prevalence, min_abundance, threshold, min_batches,
                                     batch_column = "shipment_batch") {  # ✅ 파라미터 추가
  
  report_lines <- c()
  report_lines <- c(report_lines, "=============================================================")
  report_lines <- c(report_lines, "DECONTAMINATION PIPELINE SUMMARY REPORT")
  report_lines <- c(report_lines, "=============================================================")
  report_lines <- c(report_lines, paste("Generated:", Sys.time()))
  report_lines <- c(report_lines, "")
  
  # ===== Section 0: Parameters =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, "PIPELINE PARAMETERS")
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, sprintf("Prevalence threshold:              %.1f%%", min_prevalence * 100))
  report_lines <- c(report_lines, sprintf("Abundance threshold:               %d", min_abundance))
  report_lines <- c(report_lines, sprintf("Decontam threshold:                %.2f", threshold))
  report_lines <- c(report_lines, sprintf("Minimum batches:                   %d", min_batches))
  report_lines <- c(report_lines, "")
  
  # ===== Section 1: Taxa Filtering =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, "1. TAXA FILTERING STATISTICS")
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, sprintf("Initial taxa count:                %5d", nrow(otu_original)))
  report_lines <- c(report_lines, sprintf("After hard-coded removal:          %5d (-%d)", 
                                           nrow(otu_pre_filtered), 
                                           nrow(otu_original) - nrow(otu_pre_filtered)))
  report_lines <- c(report_lines, sprintf("After rare taxa filtering:         %5d (-%d)", 
                                           nrow(otu_clean), 
                                           nrow(otu_pre_filtered) - nrow(otu_clean)))
  report_lines <- c(report_lines, sprintf("  (prevalence < %.0f%% OR abundance < %d)", 
                                           min_prevalence * 100, min_abundance))  
  report_lines <- c(report_lines, sprintf("After decontamination:             %5d (-%d)", 
                                           nrow(otu_filtered[rowSums(otu_filtered) > 0, ]), 
                                           nrow(otu_clean) - nrow(otu_filtered[rowSums(otu_filtered) > 0, ])))
  report_lines <- c(report_lines, sprintf("Final taxa count:                  %5d", 
                                           nrow(otu_filtered[rowSums(otu_filtered) > 0, ])))
  report_lines <- c(report_lines, sprintf("Total taxa removed:                %5d (%.1f%%)", 
                                           nrow(otu_original) - nrow(otu_filtered[rowSums(otu_filtered) > 0, ]),
                                           100 * (nrow(otu_original) - nrow(otu_filtered[rowSums(otu_filtered) > 0, ])) / nrow(otu_original)))
  report_lines <- c(report_lines, "")
  
  # ===== Section 2: Sample Distribution =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, "2. SAMPLE DISTRIBUTION")
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, sprintf("Total samples:                     %5d", ncol(otu_original)))
  
  type_counts <- table(meta_df$Type)
  report_lines <- c(report_lines, sprintf("  - Tumor:                         %5d (%.1f%%)", 
                                           type_counts["Tumor"], 
                                           100 * type_counts["Tumor"] / ncol(otu_original)))
  report_lines <- c(report_lines, sprintf("  - Control (Blood):               %5d (%.1f%%)", 
                                           type_counts["Control"], 
                                           100 * type_counts["Control"] / ncol(otu_original)))
  if("Solid Tissue Normal" %in% names(type_counts)) {
    report_lines <- c(report_lines, sprintf("  - Solid Tissue Normal:           %5d (%.1f%%)", 
                                             type_counts["Solid Tissue Normal"], 
                                             100 * type_counts["Solid Tissue Normal"] / ncol(otu_original)))
  }
  report_lines <- c(report_lines, "")
  
  # ===== Section 3: Batch-wise Distribution =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, "3. BATCH-WISE SAMPLE DISTRIBUTION")
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  
  batch_summary <- meta_df %>%
    filter(Type %in% c("Tumor", "Control")) %>%
    group_by(!!sym(batch_column), Type) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Type, values_from = n, values_fill = 0) %>%
    mutate(
      Tumor = if("Tumor" %in% names(.)) Tumor else 0,
      Control = if("Control" %in% names(.)) Control else 0,
      Total = Tumor + Control,
      Usable = (Control > 0 & Tumor > 0 & Total > 2)
    ) %>%
    arrange(desc(Usable), desc(Total))
  
  report_lines <- c(report_lines, sprintf("Total batches: %d", nrow(batch_summary)))
  report_lines <- c(report_lines, sprintf("Usable batches (for decontam): %d", sum(batch_summary$Usable)))
  report_lines <- c(report_lines, "")
  
  report_lines <- c(report_lines, sprintf("%-50s %6s %8s %6s %7s", 
                                          "Batch Name", "Tumor", "Control", "Total", "Usable"))
  report_lines <- c(report_lines, strrep("-", 80))
  
  for(i in 1:nrow(batch_summary)) {
    row <- batch_summary[i, ]
    batch_name <- substr(row[[batch_column]], 1, 50)
    report_lines <- c(report_lines, sprintf("%-50s %6d %8d %6d %7s",
                                            batch_name,
                                            row$Tumor,
                                            row$Control,
                                            row$Total,
                                            ifelse(row$Usable, "Yes", "No")))
  }
  report_lines <- c(report_lines, "")
  
  # ===== Section 4: Contaminant Identification =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, "4. CONTAMINANT IDENTIFICATION")
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  
  if(nrow(master_contaminant_table_raw) > 0) {
    contam_per_batch <- master_contaminant_table_raw %>%
      group_by(Batch_Count) %>%
      summarise(n_contaminants = n_distinct(Contaminant), .groups = 'drop') %>%
      arrange(desc(n_contaminants))
    
    report_lines <- c(report_lines, sprintf("Total unique contaminants identified (any batch): %d", 
                                            length(unique(master_contaminant_table_raw$Contaminant))))
    report_lines <- c(report_lines, sprintf("Batches with contaminants found: %d", nrow(contam_per_batch)))
    report_lines <- c(report_lines, "")
    
    report_lines <- c(report_lines, "Contaminants per batch:")
    report_lines <- c(report_lines, sprintf("%-50s %15s", "Batch Name", "N Contaminants"))
    report_lines <- c(report_lines, strrep("-", 70))
    
    for(i in 1:nrow(contam_per_batch)) {
      batch_name <- substr(contam_per_batch$Batch_Count[i], 1, 50)
      report_lines <- c(report_lines, sprintf("%-50s %15d", 
                                              batch_name,
                                              contam_per_batch$n_contaminants[i]))
    }
  } else {
    report_lines <- c(report_lines, "No contaminants identified in any batch.")
  }
  report_lines <- c(report_lines, "")
  
  # ===== Section 5: Final Filtering =====
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  report_lines <- c(report_lines, sprintf("5. FINAL CONTAMINANT FILTERING (≥%d BATCHES)", min_batches))  # ✅ 동적
  report_lines <- c(report_lines, "-------------------------------------------------------------")
  
  if(nrow(final_contaminants_df) > 0) {
    report_lines <- c(report_lines, sprintf("Total contaminants after filtering: %d", 
                                            nrow(final_contaminants_df)))
    report_lines <- c(report_lines, sprintf("Contaminants removed (found in <%d batches): %d",  # ✅ 동적
                                            min_batches,
                                            length(unique(master_contaminant_table_raw$Contaminant)) - nrow(final_contaminants_df)))
    report_lines <- c(report_lines, "")
    
    report_lines <- c(report_lines, "Final contaminant list:")
    report_lines <- c(report_lines, sprintf("%-60s %8s %s", 
                                            "Contaminant Name", "N Batches", "Batch List"))
    report_lines <- c(report_lines, strrep("-", 100))
    
    final_sorted <- final_contaminants_df %>% arrange(desc(n_batches))
    for(i in 1:min(50, nrow(final_sorted))) {  # 최대 50개만 표시
      row <- final_sorted[i, ]
      contam_name <- substr(row$Contaminant, 1, 60)
      batch_list <- substr(row$batch_list, 1, 30)
      report_lines <- c(report_lines, sprintf("%-60s %8d %s", 
                                              contam_name,
                                              row$n_batches,
                                              batch_list))
    }
    
    if(nrow(final_sorted) > 50) {
      report_lines <- c(report_lines, sprintf("... and %d more contaminants", nrow(final_sorted) - 50))
    }
  } else {
    report_lines <- c(report_lines, "No contaminants passed the filtering criteria.")
  }
  report_lines <- c(report_lines, "")
  
  # ===== Section 6: Summary =====
  report_lines <- c(report_lines, "=============================================================")
  report_lines <- c(report_lines, "SUMMARY")
  report_lines <- c(report_lines, "=============================================================")
  report_lines <- c(report_lines, sprintf("Initial taxa:          %5d", nrow(otu_original)))
  report_lines <- c(report_lines, sprintf("Final taxa:            %5d (%.1f%% retained)", 
                                          nrow(otu_filtered[rowSums(otu_filtered) > 0, ]),
                                          100 * nrow(otu_filtered[rowSums(otu_filtered) > 0, ]) / nrow(otu_original)))
  report_lines <- c(report_lines, sprintf("Total samples:         %5d", ncol(otu_original)))
  report_lines <- c(report_lines, sprintf("Contaminants removed:  %5d", nrow(final_contaminants_df)))
  report_lines <- c(report_lines, "=============================================================")
  report_lines <- c(report_lines, "")
  
  # Print to console
  cat("\n")
  for(line in report_lines) {
    cat(line, "\n")
  }
  
  # Save to file
  report_file <- paste0(output_prefix, ".decontamination_summary.txt")
  writeLines(report_lines, report_file)
  cat("\nINFO: Summary report saved to:", report_file, "\n")
  
  return(invisible(report_lines))
}
                                     
# ==============================================================================
# SECTION 2: Main Workflow Execution
# ==============================================================================

# --- Define Blacklist ---
# Manually define a list of species to always highlight on plots.
# This could also be read from an external file.
blacklist_species <- c("Rhodococcus fascians", "Streptococcus sanguinis", "Kocuria palustris", "Pseudomonas putida", "Psychrobacter pulmonis", "Geobacillus vulcani", "Bosea lupini", "Ensifer adhaerens", "Novosphingobium pentaromativorans", "Methylobacterium aminovorans", "Brevundimonas aurantiaca", "Caulobacter vibrioides", "Psychrobacter maritimus", "Sphingomonas faeni", "Sphingomonas koreensis", "Acinetobacter schindleri", "Ensifer meliloti", "Algoriphagus aquaeductus", "Shinella zoogloeoides", "Sphingopyxis macrogoltabida", "Pseudomonas fragi", "Ochrobactrum tritici", "Phenylobacterium haematophilum", "Comamonas aquatica", "Acinetobacter baumannii", "Acinetobacter calcoaceticus", "Sphingopyxis alaskensis", "Bradyrhizobium liaoningense", "Flavobacterium lindanitolerans", "Ralstonia insidiosa", "Pseudomonas stutzeri", "Acinetobacter junii", "Microbacterium oxydans", "Microbacterium chocolatum", "Brevibacterium epidermidis", "Ochrobactrum intermedium", "Lactococcus lactis", "Staphylococcus hominis", "Shinella granuli", "Streptococcus salivarius", "Sphingomonas mucosissima", "Methylobacterium dankookense", "Haemophilus parainfluenzae", "Sphingobacterium spiritivorum", "Sphingomonas leidyi", "Caulobacter leidyi", "Acinetobacter parvus", "Methylobacterium jeotgali", "Methylobacterium adhaesivum", "Corynebacterium kroppenstedtii", "Pseudomonas brenneri", "Acinetobacter towneri", "Pelomonas aquatica", "Methylobacterium oryzae", "Acinetobacter lwoffii", "Methylobacterium oxalidis", "Streptococcus mitis", "Corynebacterium tuberculostearicum", "Kocuria rhizophila", "Acidovorax defluvii", "Pseudomonas veronii", "Delftia acidovorans", "Afipia broomeae", "Pedomicrobium australicum", "Rhizobium radiobacter", "Bosea vestrisii", "Bradyrhizobium daqingense", "Staphylococcus epidermidis", "Sphingomonas yabuuchiae", "Delftia tsuruhatensis", "Brevundimonas vesicularis", "Bradyrhizobium denitrificans", "Escherichia flexneri", "Shigella flexneri", "Halomonas axialensis", "Methylobacterium radiotolerans", "Enhydrobacter aerosaccus", "Ralstonia pickettii", "Halomonas aquamarina", "Micrococcus luteus", "Paracoccus aminovorans", "Halomonas sulfidaeris", "Bradyrhizobium japonicum", "Pseudomonas fluorescens", "Halomonas meridiana", "Shewanella algae", "Brevundimonas diminuta", "Acinetobacter johnsonii", "Bradyrhizobium elkanii", "Sphingomonas echinoides", "Acinetobacter guillouiae", "Stenotrophomonas maltophilia", "Variovorax paradoxus", "Propionibacterium acnes", "Pelomonas puraquae",
"Rhodococcus_fascians", "Streptococcus_sanguinis", "Kocuria_palustris", "Pseudomonas_putida", "Psychrobacter_pulmonis", "Geobacillus_vulcani", "Bosea_lupini", "Ensifer_adhaerens", "Novosphingobium_pentaromativorans", "Methylobacterium_aminovorans", "Brevundimonas_aurantiaca", "Caulobacter_vibrioides", "Psychrobacter_maritimus", "Sphingomonas_faeni", "Sphingomonas_koreensis", "Acinetobacter_schindleri", "Ensifer_meliloti", "Algoriphagus_aquaeductus", "Shinella_zoogloeoides", "Sphingopyxis_macrogoltabida", "Pseudomonas_fragi", "Ochrobactrum_tritici", "Phenylobacterium_haematophilum", "Comamonas_aquatica", "Acinetobacter_baumannii", "Acinetobacter_calcoaceticus", "Sphingopyxis_alaskensis", "Bradyrhizobium_liaoningense", "Flavobacterium_lindanitolerans", "Ralstonia_insidiosa", "Pseudomonas_stutzeri", "Acinetobacter_junii", "Microbacterium_oxydans", "Microbacterium_chocolatum", "Brevibacterium_epidermidis", "Ochrobactrum_intermedium", "Lactococcus_lactis", "Staphylococcus_hominis", "Shinella_granuli", "Streptococcus_salivarius", "Sphingomonas_mucosissima", "Methylobacterium_dankookense", "Haemophilus_parainfluenzae", "Sphingobacterium_spiritivorum", "Sphingomonas_leidyi", "Caulobacter_leidyi", "Acinetobacter_parvus", "Methylobacterium_jeotgali", "Methylobacterium_adhaesivum", "Corynebacterium_kroppenstedtii", "Pseudomonas_brenneri", "Acinetobacter_towneri", "Pelomonas_aquatica", "Methylobacterium_oryzae", "Acinetobacter_lwoffii", "Methylobacterium_oxalidis", "Streptococcus_mitis", "Corynebacterium_tuberculostearicum", "Kocuria_rhizophila", "Acidovorax_defluvii", "Pseudomonas_veronii", "Delftia_acidovorans", "Afipia_broomeae", "Pedomicrobium_australicum", "Rhizobium_radiobacter", "Bosea_vestrisii", "Bradyrhizobium_daqingense", "Staphylococcus_epidermidis", "Sphingomonas_yabuuchiae", "Delftia_tsuruhatensis", "Brevundimonas_vesicularis", "Bradyrhizobium_denitrificans", "Escherichia_flexneri", "Shigella_flexneri", "Halomonas_axialensis", "Methylobacterium_radiotolerans", "Enhydrobacter_aerosaccus", "Ralstonia_pickettii", "Halomonas_aquamarina", "Micrococcus_luteus", "Paracoccus_aminovorans", "Halomonas_sulfidaeris", "Bradyrhizobium_japonicum", "Pseudomonas_fluorescens", "Halomonas_meridiana", "Shewanella_algae", "Brevundimonas_diminuta", "Acinetobacter_johnsonii", "Bradyrhizobium_elkanii", "Sphingomonas_echinoides", "Acinetobacter_guillouiae", "Stenotrophomonas_maltophilia", "Variovorax_paradoxus", "Propionibacterium_acnes", "Pelomonas_puraquae")

# Species in this list will be removed from the OTU table before any analysis.
# This is useful for removing known, persistent lab contaminants or host DNA.
hardcoded_removal_list <- c(
    "Homo sapiens", # Example: remove human reads
    "Homo_sapiens", # Example: remove human reads
    "PhiX"          # Example: remove Phix spike-in
)

force_interactive <- FALSE

# --- Argument Parsing ---
if(interactive() || force_interactive) {
    args <- list(
        otu_table = '/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/analysis/CRC/Manuscript/Results/TCGA/consensus.taxa/tumor_normal_tcga.crc.consensus.species.v2.txt',
        metadata = "/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/analysis/CRC/Manuscript/Tables/TCGA/tcga.crc.metadata.v2.txt",
        prefix = "/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/02.Results/251107_TCGA_CRC_preprocessing_posthoc/251107_5p_5read_filtered_strict/decontam_result",
        threshold = 0.1, # default 0.1, strict 0.5
        min_prevalence = 0.05,
        min_abundance = 5,
        min_batches = 2,
        batch_column = "shipment_batch"
        )
} else{
    parser <- ArgumentParser(description="Perform decontam-based decontamination on an OTU table.")
    parser$add_argument("--otu_table", type="character", required=TRUE, help="Path to the input OTU table file.")
    parser$add_argument("--metadata", type="character", required=TRUE, help="Path to the metadata file.")
    parser$add_argument("--prefix", type="character", required=TRUE, help="Prefix for the output files.")
    parser$add_argument("--threshold", type="double", default=0.5, help="Threshold for decontam package [default: 0.5]")
    parser$add_argument("--min_prevalence", type="double", default=0.05, help="Minimum prevalence threshold to keep an OTU [default: 0.05]")
    parser$add_argument("--min_abundance", type="double", default=10, help="Minimum total read count to keep an OTU [default: 10]")
    parser$add_argument("--min_batches", type="integer", default=2, help="Minimum number of batches in which a taxon must be identified as contaminant [default: 2]")
    parser$add_argument("--batch_column", type="character", default="shipment_batch", help="Name of the column in metadata that specifies batch [default: shipment_batch]")
    parser$add_argument("--type_column", type="character", default="cohort", help="Name of the column in metadata that specifies sample type [default: cohort]")
    parser$add_argument("--tumor_values", type="character", default="Cancer,Precancer,Tumor,Primary Tumor", help="Comma-separated values to be classified as Tumor [default: Cancer,Precancer,Tumor,Primary Tumor]")
    parser$add_argument("--control_values", type="character", default="Normal,Healthy,Control,Blood", help="Comma-separated values to be classified as Control [default: Normal,Healthy,Control,Blood]")
    args <- parser$parse_args()
}                             
                                     
# --- Create Directories ---
# Create a subdirectory for plots in the current working directory
# (Nextflow will handle the output directory)
plots_dir <- "decontamination_plots"
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# --- Set Global Options for Type Mapping ---
# These options are used by load_and_prepare_data and create_or_expand_metadata
options(decontam_type_column = args$type_column)
options(decontam_tumor_values = strsplit(args$tumor_values, ",")[[1]])
options(decontam_control_values = strsplit(args$control_values, ",")[[1]])

cat("\n=== TYPE MAPPING CONFIGURATION ===\n")
cat("Type column:", args$type_column, "\n")
cat("Tumor values:", args$tumor_values, "\n")
cat("Control values:", args$control_values, "\n")
cat("==================================\n\n")

# --- Initial Data Preparation ---
# Use the helper function to load and preprocess data one time
prepared_data <- load_and_prepare_data(args$otu_table, args$metadata, batch_column = args$batch_column)
otu_original <- prepared_data$otu
meta_original <- prepared_data$meta
original_taxon_col_name <- prepared_data$name_col

# --- Pre-filtering Step: Hard-coded Removal ---
cat("\nINFO: Applying hard-coded removal filter...\n")
initial_taxa_count <- nrow(otu_original)
removal_pattern <- paste(hardcoded_removal_list, collapse = "|")
otu_pre_filtered <- otu_original[!grepl(removal_pattern, rownames(otu_original), ignore.case = TRUE), ]
cat("INFO: Removed", initial_taxa_count - nrow(otu_pre_filtered), "species based on the hard-coded list.\n")

# --- Filter Rare OTUs ---
cat("\nINFO: Filtering rare OTUs with a prevalence threshold of", args$min_prevalence * 100, "%...\n")

otu_pre_filtered_numeric <- otu_pre_filtered %>% mutate(across(everything(), as.numeric))
total_samples <- ncol(otu_pre_filtered_numeric)
min_sample_count <- total_samples * args$min_prevalence
presence_counts <- rowSums(otu_pre_filtered_numeric > 0)
                                     
total_abundance <- rowSums(otu_pre_filtered_numeric, na.rm = TRUE)                                     
rare_taxa <- rownames(otu_pre_filtered)[presence_counts < min_sample_count]
                                     
otu_clean <- filter_rare_otu(otu_pre_filtered, prevalence_threshold = args$min_prevalence, min_abundance = args$min_abundance, plots_dir = plots_dir)

# --- Step 1: Batch-specific Contaminant Identification ---
contaminants_decontam_raw <- identify_contaminants_decontam(
  otu_clean,
  meta_original,
  blacklist_species,
  plots_dir,
  threshold = args$threshold,
  batch_column = args$batch_column
)

# --- Step 2: Filter by Batch Frequency ---
contaminants_filtered <- filter_contaminants_by_frequency(
  contaminants_decontam_raw, 
  min_batches = args$min_batches
)

# --- Step 3: Generate Overall Cohort Scatter Plot ---
plot_overall_contaminants_scatter(
  otu_clean = otu_clean,
  meta_df = meta_original,
  master_contaminant_table_raw = contaminants_decontam_raw,
  final_contaminants_df = contaminants_filtered,
  blacklist = blacklist_species,
  plot_path = file.path(plots_dir, "scatter_overall_cohort.png")
)

# --- Step 4: Global Removal ---
otu_filtered <- remove_contaminants_global(otu_clean, contaminants_filtered)

# --- Step 5: Save Decontaminated Table ---
save_decontaminated_table(otu_filtered, original_taxon_col_name, args$prefix, "decontam_pkg")

# --- Step 6: Generate Proportion Comparison Plots ---
plot_contaminant_proportions_comparison(
  original_otu = otu_original,
  master_contaminant_table_raw = contaminants_decontam_raw,
  final_contaminants_df = contaminants_filtered,
  rare_taxa_list = rare_taxa,
  meta_df = meta_original,
  plots_dir = plots_dir,
  batch_column = args$batch_column
)

# --- Step 7: Generate Summary Report ---
generate_summary_report(
  otu_original = otu_original,
  otu_pre_filtered = otu_pre_filtered,
  otu_clean = otu_clean,
  otu_filtered = otu_filtered,
  meta_df = meta_original,
  master_contaminant_table_raw = contaminants_decontam_raw,
  final_contaminants_df = contaminants_filtered,
  output_prefix = args$prefix,
  min_prevalence = args$min_prevalence,
  min_abundance = args$min_abundance,
  threshold = args$threshold,
  min_batches = args$min_batches,
  batch_column = args$batch_column
)


# --- End of Script ---
cat("\n=============================================================\n")
cat("INFO: Decontamination script finished successfully.\n")
cat("=============================================================\n")








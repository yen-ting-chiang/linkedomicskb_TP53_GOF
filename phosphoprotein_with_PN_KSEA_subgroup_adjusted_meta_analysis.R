################################################################################
# TP53 Phosphoprotein KSEA (with Protein Normalization)
# Cross-Dataset Meta-Analysis
#
# Method: Stouffer's Weighted Z-score Method (kinase-level)
#   1. Read per-cancer KSEA results (z.score, p.value) from
#      phosphoprotein_KSEA_subgroup_adjusted/
#      {cancer_type}/
#   2. Use KSEA z.score directly as directional z-score (KSEA already
#      provides z-scores from KSEAapp), but also verify direction with
#      p.value for robustness:
#      z = sign(z.score) * abs(qnorm(p.value / 2))
#   3. Combine z-scores across cancers using Stouffer's weighted Z:
#      Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
#      where w_i = sqrt(n_i), n_i = sample size of study i
#   4. Derive meta p-value from Z_meta and apply BH FDR correction
#
# This approach corresponds to the MAPE_P (pathway-level meta-analysis)
# strategy, where enrichment statistics from individual studies are combined
# directly at the kinase level rather than first aggregating gene-level
# statistics.
#
# Methodology references:
#   - Stouffer SA, et al. The American Soldier. 1949.
#   - Shen K, Tseng GC. Meta-analysis for pathway enrichment analysis when
#     combining multiple genomic studies. Bioinformatics. 2010;26(10):1316-23.
#     PMID: 20410053
#   - Rhodes DR, et al. Large-scale meta-analysis of cancer microarray data
#     identifies common transcriptional profiles of neoplastic transformation
#     and progression. Proc Natl Acad Sci USA. 2004;101(24):9309-14.
#     PMID: 15184698
#   - Zaykin DV. Optimally weighted Z-test is a powerful method for combining
#     probabilities in meta-analysis. J Evol Biol. 2011;24(8):1836-41.
#     PMID: 21605215
#   - Wiredja DD, Koyuturk M, Chance MR. The KSEA App: a web-based tool for
#     kinase activity inference from quantitative phosphoproteomics.
#     Bioinformatics. 2017;33(21):3489-3491. PMID: 28655153
#
# Kinase-substrate data source: PhosphoSitePlus (via KSEAapp KSData)
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspot vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspot vs TP53wt,
#   DN vs TP53wt, NonDN vs TP53wt, DN vs NonDN
#
# Input:  phosphoprotein_KSEA_subgroup_adjusted/
# Output: phosphoprotein_KSEA_with_PN_meta/
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(openxlsx)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

ksea_input_dir <- file.path(base_path, "phosphoprotein_KSEA_subgroup_adjusted")

output_dir <- file.path(base_path, "phosphoprotein_KSEA_with_PN_meta")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 9 comparisons with file name patterns
comparisons <- data.frame(
  name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
           "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
           "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
  label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
            "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
            "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
  stringsAsFactors = FALSE
)

# Minimum number of cancer types required for meta-analysis
min_studies <- 2

cat("====================================================================\n")
cat("TP53 Phosphoprotein KSEA (with Protein Normalization)\n")
cat("Cross-Dataset Meta-Analysis\n")
cat("Method: Stouffer's Weighted Z-score (weight = sqrt(sample_size))\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Retrieve Sample Sizes Per Cancer Type Per Comparison
# ==============================================================================

cat("Retrieving sample sizes from TP53 classification...\n")

tp53_file <- file.path(base_path, "TP53_mutation_classification",
                       "all_CPTAC_TP53_classification.csv")
tp53_all <- read.csv(tp53_file, stringsAsFactors = FALSE)

# For each cancer type and comparison, compute the total sample size
get_comparison_sample_size <- function(tp53_ds, comp_name) {
  tp53_mut <- tp53_ds[tp53_ds$mt == 1, ]
  
  n <- switch(comp_name,
              "TP53mt_vs_TP53wt" = {
                n_mt <- sum(tp53_ds$mt == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_mt + n_wt
              },
              "MUT_GOF_vs_MUT_LOF" = {
                n_gof <- sum(tp53_mut$GOF == 1, na.rm = TRUE)
                n_lof <- sum(tp53_mut$LOF == 1, na.rm = TRUE)
                n_gof + n_lof
              },
              "Hotspot_vs_MUT_LOF" = {
                n_hot <- sum(tp53_mut$hotspot == 1, na.rm = TRUE)
                n_lof <- sum(tp53_mut$LOF == 1, na.rm = TRUE)
                n_hot + n_lof
              },
              "MUT_GOF_vs_TP53wt" = {
                n_gof <- sum(tp53_mut$GOF == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_gof + n_wt
              },
              "MUT_LOF_vs_TP53wt" = {
                n_lof <- sum(tp53_mut$LOF == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_lof + n_wt
              },
              "Hotspot_vs_TP53wt" = {
                n_hot <- sum(tp53_mut$hotspot == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_hot + n_wt
              },
              "DN_vs_TP53wt" = {
                n_dn <- sum(tp53_mut$DN == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_dn + n_wt
              },
              "NonDN_vs_TP53wt" = {
                n_ndn <- sum(tp53_mut$non_DN == 1, na.rm = TRUE)
                n_wt <- sum(tp53_ds$wt == 1, na.rm = TRUE)
                n_ndn + n_wt
              },
              "DN_vs_NonDN" = {
                n_dn <- sum(tp53_mut$DN == 1, na.rm = TRUE)
                n_ndn <- sum(tp53_mut$non_DN == 1, na.rm = TRUE)
                n_dn + n_ndn
              },
              NA_integer_
  )
  return(n)
}

# Build sample size lookup table
sample_size_table <- data.frame(
  cancer_type = character(),
  comparison = character(),
  n = integer(),
  stringsAsFactors = FALSE
)

for (ct in cancer_types) {
  tp53_ds <- tp53_all[tp53_all$cancer_type == ct, ]
  for (j in seq_len(nrow(comparisons))) {
    n <- get_comparison_sample_size(tp53_ds, comparisons$name[j])
    sample_size_table <- rbind(sample_size_table, data.frame(
      cancer_type = ct,
      comparison = comparisons$name[j],
      n = n,
      stringsAsFactors = FALSE
    ))
  }
}

cat("  Sample size table built:", nrow(sample_size_table), "entries\n\n")

# ==============================================================================
# Section 3: Utility Functions
# ==============================================================================

#' Convert KSEA z.score direction and p-value to a directional z-score
#'
#' Uses: z = sign(z.score) * abs(qnorm(p.value / 2))
#' This preserves the direction of kinase activity change while mapping
#' to a standard N(0,1) scale for Stouffer's method
#'
#' @param z_score Numeric vector of KSEA z-scores
#' @param p_val Numeric vector of p-values (from KSEAapp)
#' @return Numeric vector of standardized directional z-scores
ksea_to_z <- function(z_score, p_val) {
  # Clamp p-values to avoid Inf z-scores
  p_val[p_val < .Machine$double.xmin] <- .Machine$double.xmin
  p_val[p_val > 1] <- 1
  
  z <- sign(z_score) * abs(qnorm(p_val / 2))
  
  # Handle edge cases
  z[!is.finite(z)] <- NA_real_
  z
}

#' Stouffer's weighted Z-score meta-analysis for a single kinase
#'
#' Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
#' Weight w_i = sqrt(n_i)
#'
#' @param z_vec Numeric vector of z-scores (one per study)
#' @param n_vec Numeric vector of sample sizes (one per study)
#' @return Named list: Z_meta, p_meta, n_studies, direction, total_n
stouffer_weighted_z <- function(z_vec, n_vec) {
  # Remove NA entries
  valid <- is.finite(z_vec) & is.finite(n_vec) & n_vec > 0
  z_vec <- z_vec[valid]
  n_vec <- n_vec[valid]
  
  k <- length(z_vec)
  if (k < min_studies) {
    return(list(Z_meta = NA_real_, p_meta = NA_real_,
                n_studies = k, direction = NA_character_,
                total_n = NA_integer_))
  }
  
  w <- sqrt(n_vec)
  Z_meta <- sum(w * z_vec) / sqrt(sum(w^2))
  p_meta <- 2 * pnorm(-abs(Z_meta))
  
  direction <- ifelse(Z_meta > 0, "up", ifelse(Z_meta < 0, "down", "none"))
  
  list(Z_meta = Z_meta, p_meta = p_meta,
       n_studies = k, direction = direction,
       total_n = sum(n_vec))
}

# ==============================================================================
# Section 4: Main Meta-Analysis Loop (per comparison)
# ==============================================================================

# Storage for cross-comparison summary
all_summaries <- list()

cat("====================================================================\n")
cat("Kinase-Substrate Data Source: PhosphoSitePlus (via KSEAapp)\n")
cat("====================================================================\n")

for (j in seq_len(nrow(comparisons))) {
  comp_name <- comparisons$name[j]
  comp_label <- comparisons$label[j]
  
  cat("\n  --- Comparison:", comp_name, "---\n")
  
  # --- Step 1: Load all cancer-type KSEA results for this comparison ---
  all_ksea <- list()
  all_n <- c()
  included_cancers <- c()
  
  for (ct in cancer_types) {
    # KSEA output files are directly under {cancer_type}/ (no collection subdirectory)
    ksea_file <- file.path(ksea_input_dir, ct,
                           paste0("KSEA_", comp_name, ".csv"))
    if (!file.exists(ksea_file)) {
      cat("    [SKIP]", ct, ": KSEA file not found\n")
      next
    }
    
    ksea <- tryCatch(fread(ksea_file), error = function(e) NULL)
    if (is.null(ksea) || nrow(ksea) == 0) {
      cat("    [SKIP]", ct, ": empty or unreadable\n")
      next
    }
    
    # Verify required columns (KSEA output: Kinase.Gene, z.score, p.value)
    if (!all(c("Kinase.Gene", "z.score", "p.value") %in% colnames(ksea))) {
      cat("    [SKIP]", ct, ": missing required columns (Kinase.Gene, z.score, p.value)\n")
      next
    }
    
    # Get sample size for this cancer/comparison
    n_val <- sample_size_table$n[
      sample_size_table$cancer_type == ct &
        sample_size_table$comparison == comp_name
    ]
    if (length(n_val) == 0 || is.na(n_val) || n_val == 0) {
      cat("    [SKIP]", ct, ": sample size = 0 or unavailable\n")
      next
    }
    
    # Convert KSEA z.score + p.value to directional z-score for Stouffer's
    ksea$z <- ksea_to_z(ksea$z.score, ksea$p.value)
    
    # Rename Kinase.Gene to kinase for consistent downstream processing
    setnames(ksea, "Kinase.Gene", "kinase")
    
    # Keep relevant columns
    keep_cols <- intersect(c("kinase", "mS", "Enrichment", "m",
                             "z.score", "p.value", "FDR", "z"),
                           colnames(ksea))
    ksea <- ksea[, ..keep_cols]
    
    all_ksea[[ct]] <- ksea
    all_n[ct] <- n_val
    included_cancers <- c(included_cancers, ct)
    
    cat("    [OK]", ct, ":", nrow(ksea), "kinases, n =", n_val, "\n")
  }
  
  n_included <- length(included_cancers)
  cat("    Included cancers:", n_included, "/", length(cancer_types), "\n")
  
  if (n_included < min_studies) {
    cat("    [SKIP] Fewer than", min_studies, "studies available. Skipping.\n")
    all_summaries[[comp_name]] <- data.frame(
      comparison = comp_name,
      n_cancers_included = n_included,
      n_kinases_tested = NA,
      n_sig_005 = NA,
      n_up = NA,
      n_down = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # --- Step 2: Identify the union of all kinases across included cancers ---
  all_kinases <- unique(unlist(lapply(all_ksea, function(d) d$kinase)))
  cat("    Total unique kinases across studies:", length(all_kinases), "\n")
  
  # --- Step 3: Build per-cancer matrices (z-score, KSEA z.score, p.value) ---
  per_cancer_z <- matrix(NA_real_, nrow = length(all_kinases), ncol = n_included,
                         dimnames = list(all_kinases, included_cancers))
  per_cancer_zscore <- matrix(NA_real_, nrow = length(all_kinases), ncol = n_included,
                              dimnames = list(all_kinases, included_cancers))
  per_cancer_pval <- matrix(NA_real_, nrow = length(all_kinases), ncol = n_included,
                            dimnames = list(all_kinases, included_cancers))
  
  # Fill per-cancer matrices
  for (ct in included_cancers) {
    d <- all_ksea[[ct]]
    idx <- match(d$kinase, all_kinases)
    valid <- !is.na(idx)
    per_cancer_z[idx[valid], ct] <- d$z[valid]
    per_cancer_zscore[idx[valid], ct] <- d$z.score[valid]
    per_cancer_pval[idx[valid], ct] <- d$p.value[valid]
  }
  
  # --- Step 4: Stouffer's weighted Z for each kinase ---
  cat("    Running Stouffer's weighted Z-score meta-analysis...\n")
  
  n_vec <- all_n[included_cancers]
  
  meta_results <- data.frame(
    kinase = all_kinases,
    n_studies = integer(length(all_kinases)),
    total_n = integer(length(all_kinases)),
    Z_meta = numeric(length(all_kinases)),
    p_meta = numeric(length(all_kinases)),
    direction = character(length(all_kinases)),
    stringsAsFactors = FALSE
  )
  
  for (g in seq_along(all_kinases)) {
    z_row <- per_cancer_z[g, ]
    res <- stouffer_weighted_z(z_row, n_vec)
    meta_results$n_studies[g] <- res$n_studies
    meta_results$total_n[g] <- res$total_n
    meta_results$Z_meta[g] <- res$Z_meta
    meta_results$p_meta[g] <- res$p_meta
    meta_results$direction[g] <- res$direction
  }
  
  # --- Step 5: FDR correction ---
  meta_results$padj <- p.adjust(meta_results$p_meta, method = "BH")
  
  # Compute weighted mean KSEA z.score (weight = sqrt(n))
  w_mat <- matrix(sqrt(n_vec), nrow = length(all_kinases), ncol = n_included,
                  byrow = TRUE)
  w_mat[is.na(per_cancer_zscore)] <- NA_real_
  meta_results$mean_z.score <- rowSums(per_cancer_zscore * w_mat, na.rm = TRUE) /
    rowSums(w_mat, na.rm = TRUE)
  
  # Add per-cancer z-score columns (Stouffer's z)
  z_df <- as.data.frame(per_cancer_z)
  colnames(z_df) <- paste0("z_", included_cancers)
  
  # Add per-cancer KSEA z.score columns
  zscore_df <- as.data.frame(per_cancer_zscore)
  colnames(zscore_df) <- paste0("z.score_", included_cancers)
  
  # Add per-cancer pval columns
  pval_df <- as.data.frame(per_cancer_pval)
  colnames(pval_df) <- paste0("pval_", included_cancers)
  
  meta_results <- cbind(meta_results, z_df, zscore_df, pval_df)
  
  # --- Step 5b: Compute per-kinase up/down dataset counts and names ---
  # A dataset is "significantly up" if its per-cancer p.value < 0.05 and z.score > 0
  # A dataset is "significantly down" if its per-cancer p.value < 0.05 and z.score < 0
  up_datasets_num <- integer(nrow(meta_results))
  down_datasets_num <- integer(nrow(meta_results))
  up_datasets_name <- character(nrow(meta_results))
  down_datasets_name <- character(nrow(meta_results))
  
  for (g in seq_len(nrow(meta_results))) {
    kn <- meta_results$kinase[g]
    kn_idx <- match(kn, all_kinases)
    if (is.na(kn_idx)) next
    
    pvals_row <- per_cancer_pval[kn_idx, ]
    zscore_row <- per_cancer_zscore[kn_idx, ]
    
    sig_up <- included_cancers[!is.na(pvals_row) & !is.na(zscore_row) &
                                 pvals_row < 0.05 & zscore_row > 0]
    sig_down <- included_cancers[!is.na(pvals_row) & !is.na(zscore_row) &
                                   pvals_row < 0.05 & zscore_row < 0]
    
    up_datasets_num[g] <- length(sig_up)
    down_datasets_num[g] <- length(sig_down)
    up_datasets_name[g] <- paste(sig_up, collapse = ";")
    down_datasets_name[g] <- paste(sig_down, collapse = ";")
  }
  
  meta_results$up_datasets_num <- up_datasets_num
  meta_results$down_datasets_num <- down_datasets_num
  meta_results$up_datasets_name <- up_datasets_name
  meta_results$down_datasets_name <- down_datasets_name
  
  # --- Step 5c: Reorder columns ---
  # First 7 columns: kinase, Z_meta, padj, up_datasets_num,
  #   down_datasets_num, up_datasets_name, down_datasets_name
  priority_cols <- c("kinase", "Z_meta", "padj",
                     "up_datasets_num", "down_datasets_num",
                     "up_datasets_name", "down_datasets_name")
  remaining_cols <- setdiff(colnames(meta_results), priority_cols)
  meta_results <- meta_results[, c(priority_cols, remaining_cols)]
  
  # Sort by Z_meta descending
  meta_results <- meta_results[order(meta_results$Z_meta, decreasing = TRUE), ]
  
  # --- Step 6: Report ---
  n_sig <- sum(meta_results$padj < 0.05, na.rm = TRUE)
  n_up <- sum(meta_results$padj < 0.05 & meta_results$direction == "up", na.rm = TRUE)
  n_down <- sum(meta_results$padj < 0.05 & meta_results$direction == "down", na.rm = TRUE)
  
  cat("    Kinases tested:", nrow(meta_results), "\n")
  cat("    Significant (padj < 0.05):", n_sig, "\n")
  cat("      Up-regulated:", n_up, "\n")
  cat("      Down-regulated:", n_down, "\n")
  cat("    Included cancers:", paste(included_cancers, collapse = ", "), "\n")
  
  # --- Step 7: Save output ---
  out_csv <- file.path(output_dir, paste0("META_KSEA_", comp_name, ".csv"))
  fwrite(meta_results, out_csv)
  
  # XLSX with formatting
  wb <- createWorkbook()
  addWorksheet(wb, comp_label)
  writeData(wb, 1, meta_results)
  
  # Base style: Arial size 7, left-aligned horizontally, vertically centered
  base_style <- createStyle(
    fontName = "Arial", fontSize = 7,
    halign = "left", valign = "center"
  )
  addStyle(wb, 1, base_style,
           rows = 1:(nrow(meta_results) + 1),
           cols = 1:ncol(meta_results),
           gridExpand = TRUE, stack = FALSE)
  
  # Header style: bold, left-aligned, blue background (on top of base)
  header_style <- createStyle(
    fontName = "Arial", fontSize = 7,
    textDecoration = "bold", halign = "left", valign = "center",
    border = "bottom", fgFill = "#4472C4", fontColour = "white"
  )
  addStyle(wb, 1, header_style,
           rows = 1, cols = 1:ncol(meta_results), gridExpand = TRUE,
           stack = FALSE)
  setColWidths(wb, 1, cols = 1:ncol(meta_results), widths = "auto")
  
  # Highlight significant rows
  sig_rows <- which(meta_results$padj < 0.05) + 1  # +1 for header
  if (length(sig_rows) > 0) {
    sig_style <- createStyle(
      fontName = "Arial", fontSize = 7,
      halign = "left", valign = "center",
      fgFill = "#FFFFCC"
    )
    addStyle(wb, 1, sig_style,
             rows = sig_rows, cols = 1:ncol(meta_results), gridExpand = TRUE,
             stack = FALSE)
  }
  
  out_xlsx <- file.path(output_dir, paste0("META_KSEA_", comp_name, ".xlsx"))
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  
  cat("    Saved:", basename(out_csv), "&", basename(out_xlsx), "\n")
  
  # Store summary
  all_summaries[[comp_name]] <- data.frame(
    comparison = comp_name,
    n_cancers_included = n_included,
    cancers_included = paste(included_cancers, collapse = ";"),
    n_kinases_tested = nrow(meta_results),
    n_sig_005 = n_sig,
    n_up = n_up,
    n_down = n_down,
    stringsAsFactors = FALSE
  )
}

cat("\n")

# ==============================================================================
# Section 5: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Meta-Analysis Summary Statistics\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(all_summaries)

# CSV
fwrite(summary_df, file.path(output_dir, "META_KSEA_summary.csv"))

# XLSX
wb_sum <- createWorkbook()
addWorksheet(wb_sum, "Summary")
writeData(wb_sum, "Summary", summary_df)

# Base style: Arial size 7, left-aligned, vertically centered
base_style_sum <- createStyle(
  fontName = "Arial", fontSize = 7,
  halign = "left", valign = "center"
)
addStyle(wb_sum, "Summary", base_style_sum,
         rows = 1:(nrow(summary_df) + 1),
         cols = 1:ncol(summary_df),
         gridExpand = TRUE, stack = FALSE)

# Header style
header_style_sum <- createStyle(
  fontName = "Arial", fontSize = 7,
  textDecoration = "bold", halign = "left", valign = "center",
  border = "bottom", fgFill = "#4472C4", fontColour = "white"
)
addStyle(wb_sum, "Summary", header_style_sum,
         rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE,
         stack = FALSE)
setColWidths(wb_sum, "Summary", cols = 1:ncol(summary_df), widths = "auto")
saveWorkbook(wb_sum, file.path(output_dir, "META_KSEA_summary.xlsx"),
             overwrite = TRUE)

# Print
print(as.data.frame(summary_df), row.names = FALSE)
cat("\n")

# ==============================================================================
# Section 6: Sample Size Detail Table
# ==============================================================================

cat("\n=== Sample Size Detail (per cancer x comparison) ===\n\n")

ss_wide <- sample_size_table %>%
  pivot_wider(names_from = comparison, values_from = n)
print(as.data.frame(ss_wide), row.names = FALSE)

# Save sample size table
fwrite(ss_wide, file.path(output_dir, "sample_size_per_cancer_comparison.csv"))

cat("\n====================================================================\n")
cat("Phosphoprotein KSEA (with Protein Normalization)\n")
cat("Cross-Dataset Meta-Analysis Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

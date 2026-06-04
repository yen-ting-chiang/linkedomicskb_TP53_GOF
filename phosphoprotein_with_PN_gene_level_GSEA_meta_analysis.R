################################################################################
# TP53 Phosphoprotein (with Protein Normalization) Gene-Level GSEA
# Cross-Dataset Meta-Analysis
#
# Method: Stouffer's Weighted Z-score Method (pathway-level)
#   1. Read per-cancer gene-level GSEA results (NES, pval) from
#      phosphoprotein_with_PN_gene_level_GSEA/{agg_method}/{cancer_type}/{collection}/
#   2. Convert NES direction + pval to directional z-score:
#      z = sign(NES) * abs(qnorm(pval / 2))
#   3. Combine z-scores across cancers using Stouffer's weighted Z:
#      Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
#      where w_i = sqrt(n_i), n_i = sample size of study i
#   4. Derive meta p-value from Z_meta and apply BH FDR correction
#
# This approach corresponds to the MAPE_P (pathway-level meta-analysis)
# strategy, where enrichment statistics from individual studies are combined
# directly at the pathway level rather than first aggregating gene-level
# statistics.
#
# Two aggregation methods are processed separately:
#   - max_abs_t: gene-level t-statistic = max absolute t across phosphosites
#   - median_t:  gene-level t-statistic = median t across phosphosites
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
#
# 19 MSigDB Collections:
#   Hallmark, C2 (CGP, BIOCARTA, KEGG_LEGACY, KEGG_MEDICUS, PID, REACTOME,
#   WIKIPATHWAYS), C3 (MIR_MIRDB, MIR_MIR_LEGACY, TFT_GTRD, TFT_TFT_LEGACY),
#   C4 (3CA, CGN, CM), C5 (GO_BP, GO_CC, GO_MF), C6 (Oncogenic)
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspot vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspot vs TP53wt,
#   DN vs TP53wt, NonDN vs TP53wt, DN vs NonDN
#
# Input:  phosphoprotein_with_PN_gene_level_GSEA/{max_abs_t, median_t}/
# Output: phosphoprotein_with_PN_gene_level_GSEA_meta/{max_abs_t, median_t}/
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

# Two aggregation methods to process separately
agg_methods <- c("max_abs_t", "median_t")

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

# 19 MSigDB collections (matching the GSEA pipeline)
collection_names <- c(
  "Hallmark",
  "C2_CGP", "C2_CP_BIOCARTA", "C2_CP_KEGG_LEGACY", "C2_CP_KEGG_MEDICUS",
  "C2_CP_PID", "C2_CP_REACTOME", "C2_CP_WIKIPATHWAYS",
  "C3_MIR_MIRDB", "C3_MIR_MIR_LEGACY",
  "C3_TFT_GTRD", "C3_TFT_TFT_LEGACY",
  "C4_3CA", "C4_CGN", "C4_CM",
  "C5_GO_BP", "C5_GO_CC", "C5_GO_MF",
  "C6_Oncogenic"
)

# Minimum number of cancer types required for meta-analysis
min_studies <- 2

cat("====================================================================\n")
cat("TP53 Phosphoprotein (with PN) Gene-Level GSEA\n")
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

#' Convert NES direction and p-value to a directional z-score
#'
#' Uses: z = sign(NES) * abs(qnorm(pval / 2))
#' This preserves the direction of enrichment while mapping to N(0,1)
#'
#' @param nes Numeric vector of Normalized Enrichment Scores
#' @param p_val Numeric vector of p-values (from fgsea)
#' @return Numeric vector of z-scores
nes_to_z <- function(nes, p_val) {
  # Clamp p-values to avoid Inf z-scores
  p_val[p_val < .Machine$double.xmin] <- .Machine$double.xmin
  p_val[p_val > 1] <- 1
  
  z <- sign(nes) * abs(qnorm(p_val / 2))
  
  # Handle edge cases
  z[!is.finite(z)] <- NA_real_
  z
}

#' Stouffer's weighted Z-score meta-analysis for a single pathway
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
# Section 4: Main Meta-Analysis Loop
#   Outer loop: aggregation method (max_abs_t, median_t)
#   Inner loops: collection x comparison (same as template)
# ==============================================================================

for (agg_method in agg_methods) {
  
  cat("####################################################################\n")
  cat("Aggregation method:", agg_method, "\n")
  cat("####################################################################\n\n")
  
  # Input and output directories for this aggregation method
  gsea_input_dir <- file.path(base_path, "phosphoprotein_with_PN_gene_level_GSEA",
                              agg_method)
  output_dir <- file.path(base_path, "phosphoprotein_with_PN_gene_level_GSEA_meta",
                          agg_method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Storage for cross-comparison summary (one per collection)
  all_summaries <- list()
  for (coll_name in collection_names) {
    all_summaries[[coll_name]] <- list()
  }
  
  for (coll_name in collection_names) {
    cat("====================================================================\n")
    cat("Collection:", coll_name, " [", agg_method, "]\n")
    cat("====================================================================\n")
    
    # Create collection output subdirectory
    coll_out <- file.path(output_dir, coll_name)
    dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)
    
    for (j in seq_len(nrow(comparisons))) {
      comp_name <- comparisons$name[j]
      comp_label <- comparisons$label[j]
      
      cat("\n  --- Comparison:", comp_name, "---\n")
      
      # --- Step 1: Load all cancer-type GSEA results for this collection + comparison ---
      all_gsea <- list()
      all_n <- c()
      included_cancers <- c()
      
      for (ct in cancer_types) {
        gsea_file <- file.path(gsea_input_dir, ct, coll_name,
                               paste0("GSEA_", comp_name, ".csv"))
        if (!file.exists(gsea_file)) {
          cat("    [SKIP]", ct, ": GSEA file not found\n")
          next
        }
        
        gsea <- tryCatch(fread(gsea_file), error = function(e) NULL)
        if (is.null(gsea) || nrow(gsea) == 0) {
          cat("    [SKIP]", ct, ": empty or unreadable\n")
          next
        }
        
        # Verify required columns
        if (!all(c("pathway", "NES", "pval") %in% colnames(gsea))) {
          cat("    [SKIP]", ct, ": missing required columns (pathway, NES, pval)\n")
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
        
        # Convert NES + pval to directional z-score
        gsea$z <- nes_to_z(gsea$NES, gsea$pval)
        
        # Keep relevant columns
        keep_cols <- intersect(c("pathway", "pval", "padj", "ES", "NES", "size", "z"),
                               colnames(gsea))
        gsea <- gsea[, ..keep_cols]
        
        all_gsea[[ct]] <- gsea
        all_n[ct] <- n_val
        included_cancers <- c(included_cancers, ct)
        
        cat("    [OK]", ct, ":", nrow(gsea), "pathways, n =", n_val, "\n")
      }
      
      n_included <- length(included_cancers)
      cat("    Included cancers:", n_included, "/", length(cancer_types), "\n")
      
      if (n_included < min_studies) {
        cat("    [SKIP] Fewer than", min_studies, "studies available. Skipping.\n")
        all_summaries[[coll_name]][[comp_name]] <- data.frame(
          comparison = comp_name,
          n_cancers_included = n_included,
          n_pathways_tested = NA,
          n_sig_005 = NA,
          n_up = NA,
          n_down = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # --- Step 2: Identify the union of all pathways across included cancers ---
      all_pathways <- unique(unlist(lapply(all_gsea, function(d) d$pathway)))
      cat("    Total unique pathways across studies:", length(all_pathways), "\n")
      
      # --- Step 3: Build per-cancer matrices (z-score, NES, pval) ---
      per_cancer_z <- matrix(NA_real_, nrow = length(all_pathways), ncol = n_included,
                             dimnames = list(all_pathways, included_cancers))
      per_cancer_NES <- matrix(NA_real_, nrow = length(all_pathways), ncol = n_included,
                               dimnames = list(all_pathways, included_cancers))
      per_cancer_pval <- matrix(NA_real_, nrow = length(all_pathways), ncol = n_included,
                                dimnames = list(all_pathways, included_cancers))
      
      # Fill per-cancer matrices
      for (ct in included_cancers) {
        d <- all_gsea[[ct]]
        idx <- match(d$pathway, all_pathways)
        valid <- !is.na(idx)
        per_cancer_z[idx[valid], ct] <- d$z[valid]
        per_cancer_NES[idx[valid], ct] <- d$NES[valid]
        per_cancer_pval[idx[valid], ct] <- d$pval[valid]
      }
      
      # --- Step 4: Stouffer's weighted Z for each pathway ---
      cat("    Running Stouffer's weighted Z-score meta-analysis...\n")
      
      n_vec <- all_n[included_cancers]
      
      meta_results <- data.frame(
        pathway = all_pathways,
        n_studies = integer(length(all_pathways)),
        total_n = integer(length(all_pathways)),
        Z_meta = numeric(length(all_pathways)),
        p_meta = numeric(length(all_pathways)),
        direction = character(length(all_pathways)),
        stringsAsFactors = FALSE
      )
      
      for (g in seq_along(all_pathways)) {
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
      
      # Compute weighted mean NES (weight = sqrt(n))
      w_mat <- matrix(sqrt(n_vec), nrow = length(all_pathways), ncol = n_included,
                      byrow = TRUE)
      w_mat[is.na(per_cancer_NES)] <- NA_real_
      meta_results$mean_NES <- rowSums(per_cancer_NES * w_mat, na.rm = TRUE) /
        rowSums(w_mat, na.rm = TRUE)
      
      # Add per-cancer z-score columns
      z_df <- as.data.frame(per_cancer_z)
      colnames(z_df) <- paste0("z_", included_cancers)
      
      # Add per-cancer NES columns
      nes_df <- as.data.frame(per_cancer_NES)
      colnames(nes_df) <- paste0("NES_", included_cancers)
      
      # Add per-cancer pval columns
      pval_df <- as.data.frame(per_cancer_pval)
      colnames(pval_df) <- paste0("pval_", included_cancers)
      
      meta_results <- cbind(meta_results, z_df, nes_df, pval_df)
      
      # Sort by padj
      meta_results <- meta_results[order(meta_results$padj, meta_results$p_meta), ]
      
      # --- Step 6: Report ---
      n_sig <- sum(meta_results$padj < 0.05, na.rm = TRUE)
      n_up <- sum(meta_results$padj < 0.05 & meta_results$direction == "up", na.rm = TRUE)
      n_down <- sum(meta_results$padj < 0.05 & meta_results$direction == "down", na.rm = TRUE)
      
      cat("    Pathways tested:", nrow(meta_results), "\n")
      cat("    Significant (padj < 0.05):", n_sig, "\n")
      cat("      Up-regulated:", n_up, "\n")
      cat("      Down-regulated:", n_down, "\n")
      cat("    Included cancers:", paste(included_cancers, collapse = ", "), "\n")
      
      # --- Step 7: Save output ---
      out_csv <- file.path(coll_out, paste0("META_GSEA_", comp_name, ".csv"))
      fwrite(meta_results, out_csv)
      
      # XLSX with formatting
      wb <- createWorkbook()
      addWorksheet(wb, comp_label)
      writeData(wb, 1, meta_results)
      
      header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#4472C4", fontColour = "white"
      )
      addStyle(wb, 1, header_style,
               rows = 1, cols = 1:ncol(meta_results), gridExpand = TRUE)
      setColWidths(wb, 1, cols = 1:ncol(meta_results), widths = "auto")
      
      # Highlight significant rows
      sig_rows <- which(meta_results$padj < 0.05) + 1  # +1 for header
      if (length(sig_rows) > 0) {
        sig_style <- createStyle(fgFill = "#FFFFCC")
        addStyle(wb, 1, sig_style,
                 rows = sig_rows, cols = 1:ncol(meta_results), gridExpand = TRUE,
                 stack = TRUE)
      }
      
      out_xlsx <- file.path(coll_out, paste0("META_GSEA_", comp_name, ".xlsx"))
      saveWorkbook(wb, out_xlsx, overwrite = TRUE)
      
      cat("    Saved:", basename(out_csv), "&", basename(out_xlsx), "\n")
      
      # Store summary
      all_summaries[[coll_name]][[comp_name]] <- data.frame(
        comparison = comp_name,
        n_cancers_included = n_included,
        cancers_included = paste(included_cancers, collapse = ";"),
        n_pathways_tested = nrow(meta_results),
        n_sig_005 = n_sig,
        n_up = n_up,
        n_down = n_down,
        stringsAsFactors = FALSE
      )
    }
    cat("\n")
  }
  
  # ==========================================================================
  # Section 5: Summary Statistics (per collection, within this agg_method)
  # ==========================================================================
  
  cat("====================================================================\n")
  cat("Generating Meta-Analysis Summary Statistics [", agg_method, "]\n")
  cat("====================================================================\n\n")
  
  for (coll_name in collection_names) {
    cat("--- Summary for", coll_name, "[", agg_method, "] ---\n")
    
    summary_df <- bind_rows(all_summaries[[coll_name]])
    
    # CSV
    fwrite(summary_df, file.path(output_dir,
                                 paste0("META_GSEA_summary_", coll_name, ".csv")))
    
    # XLSX
    wb_sum <- createWorkbook()
    addWorksheet(wb_sum, "Summary")
    writeData(wb_sum, "Summary", summary_df)
    header_style <- createStyle(
      textDecoration = "bold", halign = "center",
      border = "bottom", fgFill = "#4472C4", fontColour = "white"
    )
    addStyle(wb_sum, "Summary", header_style,
             rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
    setColWidths(wb_sum, "Summary", cols = 1:ncol(summary_df), widths = "auto")
    saveWorkbook(wb_sum, file.path(output_dir,
                                   paste0("META_GSEA_summary_", coll_name, ".xlsx")), overwrite = TRUE)
    
    # Print
    print(as.data.frame(summary_df), row.names = FALSE)
    cat("\n")
  }
  
  # ==========================================================================
  # Section 6: Sample Size Detail Table (saved once per agg_method)
  # ==========================================================================
  
  cat("\n=== Sample Size Detail (per cancer x comparison) [", agg_method, "] ===\n\n")
  
  ss_wide <- sample_size_table %>%
    pivot_wider(names_from = comparison, values_from = n)
  print(as.data.frame(ss_wide), row.names = FALSE)
  
  # Save sample size table
  fwrite(ss_wide, file.path(output_dir, "sample_size_per_cancer_comparison.csv"))
  
  cat("\n====================================================================\n")
  cat("Phosphoprotein (with PN) Gene-Level GSEA Meta-Analysis Complete!\n")
  cat("Aggregation method:", agg_method, "\n")
  cat("All results saved to:", output_dir, "\n")
  cat("====================================================================\n\n")
}

cat("\n####################################################################\n")
cat("All aggregation methods processed. Meta-analysis fully complete!\n")
cat("####################################################################\n")

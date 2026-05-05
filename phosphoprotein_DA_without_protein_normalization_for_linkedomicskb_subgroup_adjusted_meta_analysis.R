################################################################################
# TP53 phosphoprotein-Level DPS Cross-Cancer Meta-Analysis
#
# Method: Stouffer's Weighted Z-score Method
#   1. Read per-cancer DPS results (t-statistic, P.Value) from
#      phosphoprotein_differential_analysis_subgroup_adjusted/
#   2. Convert t-statistic -> two-sided p-value -> directional z-score
#      z = sign(t) * abs(qnorm(P.Value / 2))
#   3. Combine z-scores across cancers using Stouffer's weighted Z:
#      Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
#      where w_i = sqrt(n_i), n_i = sample size of study i
#   4. Derive meta p-value from Z_meta and apply BH FDR correction
#
# Methodology references:
#   - Stouffer SA, et al. The American Soldier. 1949.
#   - Rhodes DR, et al. Large-scale meta-analysis of cancer microarray data
#     identifies common transcriptional profiles of neoplastic transformation
#     and progression. Proc Natl Acad Sci USA. 2004;101(24):9309-14.
#     PMID: 15184698
#   - Marot G, et al. Moderated effect size and P-value combinations for
#     microarray meta-analyses. Bioinformatics. 2009;25(20):2692-9.
#     PMID: 19628502
#   - Zaykin DV. Optimally weighted Z-test is a powerful method for combining
#     probabilities in meta-analysis. J Evol Biol. 2011;24(8):1836-41.
#     PMID: 21605215
#
# 9 Comparisons:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspot vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspot vs TP53wt,
#   DN vs TP53wt, NonDN vs TP53wt, DN vs NonDN
#
# Output: phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/phosphoprotein_DPS_without_protein_normalization_meta_analysis/
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(fgsea)
    library(readxl)
    library(KSEAapp)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

deg_dir <- file.path(base_path, "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted")

output_dir <- file.path(deg_dir, "phosphoprotein_DPS_without_protein_normalization_meta_analysis")
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
cat("TP53 phosphoprotein-Level DPS without protein normalization Cross-Cancer Meta-Analysis\n")
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
# (number of samples in both groups that were actually used in the comparison)
get_comparison_sample_size <- function(tp53_ds, comp_name) {
    # Determine which samples belong to each group
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

#' Convert t-statistic and two-sided p-value to a directional z-score
#'
#' Uses: z = sign(t) * abs(qnorm(p / 2))
#' This preserves the direction of effect while mapping to N(0,1)
#'
#' @param t_stat Numeric vector of t-statistics
#' @param p_val Numeric vector of two-sided p-values
#' @return Numeric vector of z-scores
t_to_z <- function(t_stat, p_val) {
    # Clamp p-values to avoid Inf z-scores
    p_val[p_val < .Machine$double.xmin] <- .Machine$double.xmin
    p_val[p_val > 1] <- 1

    z <- sign(t_stat) * abs(qnorm(p_val / 2))

    # Handle edge cases
    z[!is.finite(z)] <- NA_real_
    z
}

#' Stouffer's weighted Z-score meta-analysis for a single phosphosite
#'
#' Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
#' Weight w_i = sqrt(n_i)
#'
#' @param z_vec Numeric vector of z-scores (one per study)
#' @param n_vec Numeric vector of sample sizes (one per study)
#' @return Named list: Z_meta, p_meta, n_studies, direction
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
meta_summary <- list()

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Meta-analysis for:", comp_name, "\n")
    cat("==================================================================\n")

    # --- Step 1: Load all cancer-type DPS results for this comparison ---
    all_deg <- list()
    all_n <- c()
    included_cancers <- c()

    for (ct in cancer_types) {
        deg_file <- file.path(deg_dir, ct, paste0("DPS_", comp_name, ".csv"))
        if (!file.exists(deg_file)) {
            cat("  [SKIP]", ct, ": DPS file not found\n")
            next
        }

        deg <- tryCatch(fread(deg_file), error = function(e) NULL)
        if (is.null(deg) || nrow(deg) == 0) {
            cat("  [SKIP]", ct, ": empty or unreadable\n")
            next
        }

        # Verify required columns
        if (!all(c("phosphosite", "t", "P.Value") %in% colnames(deg))) {
            cat("  [SKIP]", ct, ": missing required columns (phosphosite, t, P.Value)\n")
            next
        }

        # Get sample size for this cancer/comparison
        n_val <- sample_size_table$n[
            sample_size_table$cancer_type == ct &
            sample_size_table$comparison == comp_name
        ]
        if (length(n_val) == 0 || is.na(n_val) || n_val == 0) {
            cat("  [SKIP]", ct, ": sample size = 0 or unavailable\n")
            next
        }

        # Convert t -> z
        deg$z <- t_to_z(deg$t, deg$P.Value)

        # Keep phosphosite, gene_symbol_phosphosite, logFC, z
        keep_cols <- intersect(c("gene_symbol_phosphosite", "phosphosite", "logFC", "t", "P.Value", "adj.P.Val", "z"), colnames(deg))
        deg <- deg[, ..keep_cols]

        all_deg[[ct]] <- deg
        all_n[ct] <- n_val
        included_cancers <- c(included_cancers, ct)

        cat("  [OK]", ct, ": ", nrow(deg), "phosphosites, n =", n_val,
            ", z range = [", round(min(deg$z, na.rm = TRUE), 2), ",",
            round(max(deg$z, na.rm = TRUE), 2), "]\n")
    }

    n_included <- length(included_cancers)
    cat("\n  Included cancers:", n_included, "/", length(cancer_types), "\n")

    if (n_included < min_studies) {
        cat("  [SKIP] Fewer than", min_studies, "studies available. Skipping.\n\n")
        meta_summary[[comp_name]] <- data.frame(
            comparison = comp_name,
            n_cancers_included = n_included,
            n_phosphosites_tested = NA,
            n_sig_005 = NA,
            n_up = NA,
            n_down = NA,
            stringsAsFactors = FALSE
        )
        next
    }

    # --- Step 2: Identify the union of all phosphosites across included cancers ---
    all_phosphosites <- unique(unlist(lapply(all_deg, function(d) d$phosphosite)))
    cat("  Total unique phosphosites across studies:", length(all_phosphosites), "\n")

    # --- Step 3: Build gene_symbol_phosphosite lookup (take the first non-empty mapping) ---
    gene_symbol_phosphosite_lookup <- setNames(rep(NA_character_, length(all_phosphosites)), all_phosphosites)
    for (ct in included_cancers) {
        d <- all_deg[[ct]]
        if ("gene_symbol_phosphosite" %in% colnames(d)) {
            mapped <- d$gene_symbol_phosphosite[match(all_phosphosites, d$phosphosite)]
            fill <- !is.na(mapped) & mapped != "" & is.na(gene_symbol_phosphosite_lookup)
            gene_symbol_phosphosite_lookup[fill] <- mapped[fill]
        }
    }

    # --- Step 4: Stouffer's weighted Z for each phosphosite ---
    cat("  Running Stouffer's weighted Z-score meta-analysis...\n")

    meta_results <- data.frame(
        gene_symbol_phosphosite = character(length(all_phosphosites)),
        phosphosite = all_phosphosites,
        n_studies = integer(length(all_phosphosites)),
        total_n = integer(length(all_phosphosites)),
        Z_meta = numeric(length(all_phosphosites)),
        p_meta = numeric(length(all_phosphosites)),
        direction = character(length(all_phosphosites)),
        stringsAsFactors = FALSE
    )

    # Also store per-cancer z-scores and logFC for the output
    per_cancer_z <- matrix(NA_real_, nrow = length(all_phosphosites), ncol = n_included,
                           dimnames = list(all_phosphosites, included_cancers))
    per_cancer_logFC <- matrix(NA_real_, nrow = length(all_phosphosites), ncol = n_included,
                               dimnames = list(all_phosphosites, included_cancers))

    # Fill per-cancer matrices
    for (ct in included_cancers) {
        d <- all_deg[[ct]]
        idx <- match(d$phosphosite, all_phosphosites)
        valid <- !is.na(idx)
        per_cancer_z[idx[valid], ct] <- d$z[valid]
        per_cancer_logFC[idx[valid], ct] <- d$logFC[valid]
    }

    # Compute meta Z for each phosphosite
    n_vec <- all_n[included_cancers]

    for (g in seq_along(all_phosphosites)) {
        z_row <- per_cancer_z[g, ]
        res <- stouffer_weighted_z(z_row, n_vec)
        meta_results$gene_symbol_phosphosite[g] <- ifelse(is.na(gene_symbol_phosphosite_lookup[all_phosphosites[g]]),
                                              "", gene_symbol_phosphosite_lookup[all_phosphosites[g]])
        meta_results$n_studies[g] <- res$n_studies
        meta_results$total_n[g] <- res$total_n
        meta_results$Z_meta[g] <- res$Z_meta
        meta_results$p_meta[g] <- res$p_meta
        meta_results$direction[g] <- res$direction
    }

    # --- Step 5: FDR correction ---
    meta_results$padj <- p.adjust(meta_results$p_meta, method = "BH")

    # Compute weighted mean logFC (weight = sqrt(n))
    w_mat <- matrix(sqrt(n_vec), nrow = length(all_phosphosites), ncol = n_included, byrow = TRUE)
    w_mat[is.na(per_cancer_logFC)] <- NA_real_
    meta_results$mean_logFC <- rowSums(per_cancer_logFC * w_mat, na.rm = TRUE) /
                               rowSums(w_mat, na.rm = TRUE)

    # Add per-cancer z-score columns
    z_df <- as.data.frame(per_cancer_z)
    colnames(z_df) <- paste0("z_", included_cancers)
    logFC_df <- as.data.frame(per_cancer_logFC)
    colnames(logFC_df) <- paste0("logFC_", included_cancers)

    meta_results <- cbind(meta_results, z_df, logFC_df)

    # Sort by padj
    meta_results <- meta_results[order(meta_results$padj, meta_results$p_meta), ]

    # --- Step 6: Report ---
    n_sig <- sum(meta_results$padj < 0.05, na.rm = TRUE)
    n_up <- sum(meta_results$padj < 0.05 & meta_results$direction == "up", na.rm = TRUE)
    n_down <- sum(meta_results$padj < 0.05 & meta_results$direction == "down", na.rm = TRUE)

    cat("  phosphosites tested:", nrow(meta_results), "\n")
    cat("  Significant (padj < 0.05):", n_sig, "\n")
    cat("    Up-regulated:", n_up, "\n")
    cat("    Down-regulated:", n_down, "\n")
    cat("  Included cancers:", paste(included_cancers, collapse = ", "), "\n")

    # --- Step 7: Save output ---
    out_csv <- file.path(output_dir, paste0("META_DPS_", comp_name, ".csv"))
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

    out_xlsx <- file.path(output_dir, paste0("META_DPS_", comp_name, ".xlsx"))
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    cat("  Saved:", basename(out_csv), "&", basename(out_xlsx), "\n\n")

    # Store summary
    meta_summary[[comp_name]] <- data.frame(
        comparison = comp_name,
        n_cancers_included = n_included,
        cancers_included = paste(included_cancers, collapse = ";"),
        n_phosphosites_tested = nrow(meta_results),
        n_sig_005 = n_sig,
        n_up = n_up,
        n_down = n_down,
        stringsAsFactors = FALSE
    )
}

# ==============================================================================
# Section 5: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Meta-Analysis Summary\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(meta_summary)

# CSV
fwrite(summary_df, file.path(output_dir, "META_DPS_summary.csv"))

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
saveWorkbook(wb_sum, file.path(output_dir, "META_DPS_summary.xlsx"), overwrite = TRUE)

# Print summary
cat("\n=== Cross-Cancer Meta-Analysis Summary ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

# ==============================================================================
# Section 6: Sample Size Detail Table
# ==============================================================================

cat("\n\n=== Sample Size Detail (per cancer x comparison) ===\n\n")

# Pivot sample size table for display
ss_wide <- sample_size_table %>%
    pivot_wider(names_from = comparison, values_from = n)
print(as.data.frame(ss_wide), row.names = FALSE)

# Save sample size table
fwrite(ss_wide, file.path(output_dir, "sample_size_per_cancer_comparison.csv"))

cat("\n====================================================================\n")
cat("DPS Meta-Analysis Complete. Proceeding to Meta-GSEA...\n")
cat("====================================================================\n")

# ==============================================================================
# Section 7: Load Gene Sets from PTMsigDB (for Meta-GSEA)
# ==============================================================================
# Methodology reference:
#   - Shen K, Tseng GC. Meta-analysis for pathway enrichment analysis when
#     combining multiple genomic studies. Bioinformatics. 2010;26(10):1316-23.
#     PMID: 20410053
#   This approach corresponds to MAPE_G (gene-level meta-analysis followed by
#   GSEA), where gene-level statistics are first combined across studies, and
#   the resulting meta-statistics are used as ranking for preranked GSEA.

cat("\n====================================================================\n")
cat("Meta-GSEA: Preranked GSEA on Cross-Cancer Z_meta Scores (PTMsigDB)\n")
cat("====================================================================\n\n")

cat("Loading PTMsigDB gene sets...\n")

#' Read PTMsigDB and build pathway lists keyed by GENE_SITE identifier
#' Also build a flanking-sequence-to-gene_site lookup table for ID mapping
#'
#' @param fp Path to PTMsigDB xlsx file
#' @return A list with two elements:
#'   - pathways: named list of character vectors (signature -> gene_site IDs)
#'   - flank_lookup: named character vector (flanking_seq -> gene_site)
read_ptmsigdb <- function(fp) {
    if (!file.exists(fp)) {
        stop("[PTMsigDB] File does not exist: ", fp)
    }

    df <- readxl::read_xlsx(fp)
    req <- c("signature", "site.annotation", "site.flanking")
    if (!all(req %in% names(df))) {
        stop("[PTMsigDB] xlsx is missing necessary fields: ",
             paste(setdiff(req, names(df)), collapse = ", "))
    }

    # Extract gene_site from site.annotation (e.g., "PPP1R12A_T696:15226371" -> "PPP1R12A_T696")
    gene_site <- toupper(sub(":.*$", "", trimws(df$site.annotation)))

    # Extract flanking sequence (15-mer)
    flanking <- toupper(trimws(df$site.flanking))

    # Build pathway lists keyed by gene_site
    by_sig <- split(gene_site, df$signature)
    pathways <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
    pathways[lengths(pathways) == 0] <- NULL

    # Build flanking-to-gene_site lookup
    # Use unique flanking sequences; for duplicates, keep first occurrence
    valid <- nzchar(flanking) & nzchar(gene_site) & nchar(flanking) == 15
    flank_lookup <- setNames(gene_site[valid], flanking[valid])
    flank_lookup <- flank_lookup[!duplicated(names(flank_lookup))]

    list(pathways = pathways, flank_lookup = flank_lookup)
}

ptmsigdb_file <- file.path(base_path, "data_PTMsigDB_all_sites_v2.0.0.xlsx")
ptmsigdb_data <- read_ptmsigdb(ptmsigdb_file)
pathways_ptm <- ptmsigdb_data$pathways
flank_lookup <- ptmsigdb_data$flank_lookup

cat("  PTMsigDB:", length(pathways_ptm), "phosphosite sets\n")
cat("  Flanking lookup table:", length(flank_lookup), "unique 15-mer entries\n\n")

# Named list for iteration (single collection for phosphoproteomics)
gsea_collections <- list(
    PTMsigDB = pathways_ptm
)

# fgsea parameters (same as per-cancer phosphoprotein without protein normalization GSEA pipeline)
gsea_minSize <- 5
gsea_maxSize <- 500

# ==============================================================================
# Section 8: Meta-GSEA Function (Phosphoprotein without protein normalization)
# ==============================================================================

#' Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
#'
#' DPS ID format: ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
#' This function extracts the 15-mer flanking sequence (field 4) and uses the
#' prebuilt lookup table to find the corresponding PTMsigDB gene_site ID
#' (e.g., "PPP1R12A_T696").
#'
#' @param dps_ids Character vector of DPS phosphosite IDs
#' @param lookup Named character vector (flanking -> gene_site)
#' @return Character vector of PTMsigDB gene_site IDs (NA if no match)
map_dps_to_ptmsigdb <- function(dps_ids, lookup) {
    parts <- strsplit(dps_ids, "|", fixed = TRUE)

    # Extract flanking sequence (field 4) and convert to uppercase
    flanking <- vapply(parts, function(x) {
        if (length(x) >= 4) toupper(x[4]) else NA_character_
    }, character(1))

    # Look up gene_site via flanking sequence
    mapped <- lookup[flanking]

    # Report mapping statistics
    n_total <- length(dps_ids)
    n_mapped <- sum(!is.na(mapped))
    cat(sprintf("    ID mapping: %d of %d phosphosites mapped to PTMsigDB (%.1f%%)\n",
                n_mapped, n_total, 100 * n_mapped / n_total))

    mapped
}

#' Run fgsea preranked GSEA from a META_DPS file using Z_meta as ranking
#' @param meta_file Path to META_DPS CSV file
#' @param pathways Named list of phosphosite sets (gene_site IDs)
#' @param lookup Flanking-to-gene_site lookup table
#' @return fgsea result data.frame or NULL
run_meta_gsea <- function(meta_file, pathways, lookup) {
    if (!file.exists(meta_file)) return(NULL)

    meta <- fread(meta_file)
    if (!("phosphosite" %in% colnames(meta)) || !("Z_meta" %in% colnames(meta))) {
        cat("    [WARN] Missing phosphosite/Z_meta columns in:", basename(meta_file), "\n")
        return(NULL)
    }

    # Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
    gene_site_ids <- map_dps_to_ptmsigdb(meta$phosphosite, lookup)

    # Keep only successfully mapped phosphosites
    mapped_mask <- !is.na(gene_site_ids)
    if (sum(mapped_mask) < 50) {
        cat("    [WARN] Too few mapped phosphosites:", sum(mapped_mask), "\n")
        return(NULL)
    }

    mapped_ids <- gene_site_ids[mapped_mask]
    mapped_z <- meta$Z_meta[mapped_mask]

    # Build named vector of Z_meta statistics
    stats_vec <- setNames(mapped_z, mapped_ids)

    # Remove duplicates (keep first = most significant by padj ordering)
    stats_vec <- stats_vec[!duplicated(names(stats_vec))]

    # Remove non-finite values
    stats_vec <- stats_vec[is.finite(stats_vec)]

    if (length(stats_vec) < 50) {
        cat("    [WARN] Too few unique mapped phosphosites after dedup:", length(stats_vec), "\n")
        return(NULL)
    }

    # Sort descending (required by fgsea)
    stats_vec <- sort(stats_vec, decreasing = TRUE)

    # Run fgsea
    res <- tryCatch({
        suppressWarnings(fgsea::fgseaMultilevel(
            pathways = pathways,
            stats = stats_vec,
            minSize = gsea_minSize,
            maxSize = gsea_maxSize,
            eps = 0
        ))
    }, error = function(e) {
        cat("    [ERROR] fgsea failed:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(res) || nrow(res) == 0) return(NULL)

    # Convert leadingEdge list to semicolon-separated string for CSV/XLSX output
    res$leadingEdge <- vapply(res$leadingEdge,
                              function(x) paste(x, collapse = ";"),
                              character(1))

    # Sort by padj
    res <- res[order(res$padj), ]
    as.data.frame(res)
}

# ==============================================================================
# Section 9: Run Meta-GSEA for Each Comparison
# ==============================================================================

gsea_output_dir <- file.path(output_dir, "meta_GSEA")
dir.create(gsea_output_dir, recursive = TRUE, showWarnings = FALSE)

# Storage for GSEA summary statistics (one per collection)
all_gsea_summaries <- list()
for (coll_name in names(gsea_collections)) {
    all_gsea_summaries[[coll_name]] <- list()
}

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Meta-GSEA for:", comp_name, "\n")
    cat("==================================================================\n")

    meta_file <- file.path(output_dir, paste0("META_DPS_", comp_name, ".csv"))
    if (!file.exists(meta_file)) {
        cat("  [SKIP] META_DPS file not found:", basename(meta_file), "\n\n")
        for (coll_name in names(gsea_collections)) {
            all_gsea_summaries[[coll_name]][[comp_name]] <- data.frame(
                comparison = comp_name,
                stringsAsFactors = FALSE
            )
        }
        next
    }

    for (coll_name in names(gsea_collections)) {
        pathways <- gsea_collections[[coll_name]]
        cat("\n  --- Collection:", coll_name, "---\n")

        # Create collection output subdirectory
        coll_out <- file.path(gsea_output_dir, coll_name)
        dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)

        # Summary row for this comparison + collection
        summary_row <- list(comparison = comp_name)

        cat("    ", comp_label, ":", sep = "")

        res <- run_meta_gsea(meta_file, pathways, flank_lookup)

        if (!is.null(res)) {
            n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
            n_up <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
            n_down <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
            cat(" ", nrow(res), "sets tested,", n_sig, "significant\n")

            # Save CSV
            fwrite(res, file.path(coll_out, paste0("META_GSEA_", comp_name, ".csv")))

            # Save XLSX
            wb <- createWorkbook()
            addWorksheet(wb, comp_label)
            writeData(wb, 1, res)
            saveWorkbook(wb, file.path(coll_out, paste0("META_GSEA_", comp_name, ".xlsx")),
                         overwrite = TRUE)

            summary_row[[paste0(comp_name, "_total")]] <- nrow(res)
            summary_row[[paste0(comp_name, "_sig")]] <- n_sig
            summary_row[[paste0(comp_name, "_up")]] <- n_up
            summary_row[[paste0(comp_name, "_down")]] <- n_down
        } else {
            cat(" no META_DPS file or insufficient data\n")
            summary_row[[paste0(comp_name, "_total")]] <- NA
            summary_row[[paste0(comp_name, "_sig")]] <- NA
            summary_row[[paste0(comp_name, "_up")]] <- NA
            summary_row[[paste0(comp_name, "_down")]] <- NA
        }

        all_gsea_summaries[[coll_name]][[comp_name]] <-
            as.data.frame(summary_row, stringsAsFactors = FALSE)
    }
    cat("\n")
}

# ==============================================================================
# Section 10: Meta-GSEA Summary Statistics (per collection)
# ==============================================================================

cat("====================================================================\n")
cat("Generating Meta-GSEA Summary Statistics\n")
cat("====================================================================\n\n")

for (coll_name in names(all_gsea_summaries)) {
    cat("--- Summary for", coll_name, "---\n")

    gsea_summary_df <- bind_rows(all_gsea_summaries[[coll_name]])

    # CSV
    fwrite(gsea_summary_df, file.path(gsea_output_dir,
        paste0("META_GSEA_summary_", coll_name, ".csv")))

    # XLSX
    wb_sum <- createWorkbook()
    addWorksheet(wb_sum, "Summary")
    writeData(wb_sum, "Summary", gsea_summary_df)
    header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#4472C4", fontColour = "white"
    )
    addStyle(wb_sum, "Summary", header_style,
             rows = 1, cols = 1:ncol(gsea_summary_df), gridExpand = TRUE)
    setColWidths(wb_sum, "Summary", cols = 1:ncol(gsea_summary_df), widths = "auto")
    saveWorkbook(wb_sum, file.path(gsea_output_dir,
        paste0("META_GSEA_summary_", coll_name, ".xlsx")), overwrite = TRUE)

    # Print significant counts only
    sig_cols <- grep("_sig$", colnames(gsea_summary_df), value = TRUE)
    if (length(sig_cols) > 0) {
        print_df <- gsea_summary_df[, c("comparison", sig_cols), drop = FALSE]
        print(as.data.frame(print_df), row.names = FALSE)
    }
    cat("\n")
}

cat("====================================================================\n")
cat("Meta-GSEA Complete. Proceeding to Meta-KSEA...\n")
cat("====================================================================\n")

# ==============================================================================
# Section 11: Meta-KSEA Setup (Kinase-Substrate Enrichment Analysis)
# ==============================================================================
# Uses KSEAapp with Z_meta as the input statistic.
# The PTMsigDB flanking-sequence lookup (loaded in Section 7) is reused to map
# DPS phosphosite IDs to Gene + Residue format required by KSEAapp.

cat("\n====================================================================\n")
cat("Meta-KSEA: Kinase-Substrate Enrichment Analysis on Z_meta Scores\n")
cat("====================================================================\n\n")

# Load KSData from KSEAapp namespace
data("KSData", package = "KSEAapp", envir = environment())
cat("  Loaded KSData from KSEAapp package\n")

# KSEA Parameters (same as per-cancer KSEA pipeline)
use_networkin <- FALSE
networkin_cutoff <- 5

# ==============================================================================
# Section 12: Meta-KSEA Function and Processing Loop
# ==============================================================================

#' Run KSEA from a META_DPS file using Z_meta as input statistic
#' @param meta_file Path to META_DPS CSV file
#' @param lookup Flanking-to-Gene_Site lookup table (from PTMsigDB)
#' @return KSEA result data.frame or NULL
run_meta_ksea <- function(meta_file, lookup) {
    if (!file.exists(meta_file)) return(NULL)

    meta <- fread(meta_file)
    if (!("phosphosite" %in% colnames(meta)) || !("Z_meta" %in% colnames(meta)) || !("p_meta" %in% colnames(meta))) {
        cat("    [WARN] Missing phosphosite/Z_meta/p_meta columns in:", basename(meta_file), "\n")
        return(NULL)
    }

    # Filter valid rows
    meta <- meta[is.finite(meta$Z_meta) & !is.na(meta$phosphosite), ]

    # Extract flanking sequence from DPS phosphosite IDs
    parts <- strsplit(meta$phosphosite, "|", fixed = TRUE)
    flanking_dps <- vapply(parts, function(x) if(length(x) >= 4) toupper(x[4]) else NA_character_, character(1))

    # Map to Gene_Site via PTMsigDB flanking lookup
    gene_sites <- lookup[flanking_dps]

    # Keep only successfully mapped phosphosites
    valid_mask <- !is.na(gene_sites)
    n_mapped <- sum(valid_mask)
    cat(sprintf("    ID mapping: %d of %d phosphosites mapped (%.1f%%)\n",
                n_mapped, length(gene_sites), 100 * n_mapped / length(gene_sites)))

    if (n_mapped < 50) {
        cat("    [WARN] Too few mapped phosphosites:", n_mapped, "\n")
        return(NULL)
    }

    mapped_gs <- gene_sites[valid_mask]
    mapped_z <- meta$Z_meta[valid_mask]
    mapped_p <- meta$p_meta[valid_mask]

    # Split Gene_Site into Gene and Residue
    # Format is typically GENE_S260 or GENE_S260-P
    genes <- sub("_.*$", "", mapped_gs)
    residues <- sub("^[^_]+_", "", mapped_gs)
    residues <- sub("-P$", "", residues)  # Remove -P suffix if present

    # Build PX dataframe for KSEAapp
    PX <- data.frame(
        Protein = genes,
        Gene = genes,
        Peptide = "",
        Residue.Both = residues,
        p = mapped_p,
        FC = 2^mapped_z,  # Convert Z_meta to fold change scale
        stringsAsFactors = FALSE
    )

    # Dedup PX to ensure one entry per Gene_Residue (keep max absolute Z_meta)
    PX$abs_z <- abs(mapped_z)
    PX <- PX %>%
        group_by(Gene, Residue.Both) %>%
        slice_max(abs_z, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(-abs_z)

    if (nrow(PX) < 50) {
        cat("    [WARN] Too few unique mapped sites:", nrow(PX), "\n")
        return(NULL)
    }

    # Run KSEA.Scores
    res <- tryCatch({
        KSEA.Scores(KSData, PX, NetworKIN = use_networkin, NetworKIN.cutoff = networkin_cutoff)
    }, error = function(e) {
        cat("    [ERROR] KSEA failed:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(res) || nrow(res) == 0) return(NULL)

    # Sort by FDR
    res <- res[order(res$FDR), ]
    as.data.frame(res)
}

# KSEA output directory
ksea_output_dir <- file.path(output_dir, "meta_KSEA")
dir.create(ksea_output_dir, recursive = TRUE, showWarnings = FALSE)

# Storage for KSEA summary
all_ksea_summaries <- list()

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Meta-KSEA for:", comp_name, "\n")
    cat("==================================================================\n")

    meta_file <- file.path(output_dir, paste0("META_DPS_", comp_name, ".csv"))
    if (!file.exists(meta_file)) {
        cat("  [SKIP] META_DPS file not found:", basename(meta_file), "\n\n")
        all_ksea_summaries[[comp_name]] <- data.frame(
            comparison = comp_name,
            stringsAsFactors = FALSE
        )
        next
    }

    cat("    ", comp_label, ":", sep = "")

    res <- run_meta_ksea(meta_file, flank_lookup)

    summary_row <- list(comparison = comp_name)

    if (!is.null(res)) {
        n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
        n_up <- sum(res$FDR < 0.05 & res$z.score > 0, na.rm = TRUE)
        n_down <- sum(res$FDR < 0.05 & res$z.score < 0, na.rm = TRUE)
        cat(" ", nrow(res), "kinases tested,", n_sig, "significant\n")

        # Save CSV
        fwrite(res, file.path(ksea_output_dir, paste0("META_KSEA_", comp_name, ".csv")))

        # Save XLSX
        wb <- createWorkbook()
        addWorksheet(wb, comp_label)
        writeData(wb, 1, res)
        saveWorkbook(wb, file.path(ksea_output_dir, paste0("META_KSEA_", comp_name, ".xlsx")),
                     overwrite = TRUE)

        summary_row[[paste0(comp_name, "_total")]] <- nrow(res)
        summary_row[[paste0(comp_name, "_sig")]] <- n_sig
        summary_row[[paste0(comp_name, "_up")]] <- n_up
        summary_row[[paste0(comp_name, "_down")]] <- n_down
    } else {
        cat(" no META_DPS file or insufficient data\n")
        summary_row[[paste0(comp_name, "_total")]] <- NA
        summary_row[[paste0(comp_name, "_sig")]] <- NA
        summary_row[[paste0(comp_name, "_up")]] <- NA
        summary_row[[paste0(comp_name, "_down")]] <- NA
    }

    all_ksea_summaries[[comp_name]] <- as.data.frame(summary_row, stringsAsFactors = FALSE)
    cat("\n")
}

# ==============================================================================
# Section 13: Meta-KSEA Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Meta-KSEA Summary Statistics\n")
cat("====================================================================\n\n")

ksea_summary_df <- bind_rows(all_ksea_summaries)

# CSV
fwrite(ksea_summary_df, file.path(ksea_output_dir, "META_KSEA_summary_PhosphoSitePlus.csv"))

# XLSX
wb_ksea_sum <- createWorkbook()
addWorksheet(wb_ksea_sum, "Summary")
writeData(wb_ksea_sum, "Summary", ksea_summary_df)
header_style <- createStyle(
    textDecoration = "bold", halign = "center",
    border = "bottom", fgFill = "#4472C4", fontColour = "white"
)
addStyle(wb_ksea_sum, "Summary", header_style,
         rows = 1, cols = 1:ncol(ksea_summary_df), gridExpand = TRUE)
setColWidths(wb_ksea_sum, "Summary", cols = 1:ncol(ksea_summary_df), widths = "auto")
saveWorkbook(wb_ksea_sum, file.path(ksea_output_dir, "META_KSEA_summary_PhosphoSitePlus.xlsx"),
             overwrite = TRUE)

# Print summary
cat("--- Meta-KSEA Summary ---\n")
sig_cols <- grep("_sig$", colnames(ksea_summary_df), value = TRUE)
if (length(sig_cols) > 0) {
    print_df <- ksea_summary_df[, c("comparison", sig_cols), drop = FALSE]
    print(as.data.frame(print_df), row.names = FALSE)
}
cat("\n")

cat("====================================================================\n")
cat("Full Pipeline Complete (DPS Meta-Analysis + Meta-GSEA + Meta-KSEA)!\n")
cat("DPS Meta-Analysis results saved to:", output_dir, "\n")
cat("Meta-GSEA results saved to:", gsea_output_dir, "\n")
cat("Meta-KSEA results saved to:", ksea_output_dir, "\n")
cat("====================================================================\n")

################################################################################
# TP53 protein-Level DEG Cross-Cancer Meta-Analysis
#
# Method: Stouffer's Weighted Z-score Method
#   1. Read per-cancer DEG results (t-statistic, P.Value) from
#      protein_differential_analysis_subgroup_adjusted/
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
# Output: protein_differential_analysis_subgroup_adjusted/protein_DEG_meta_analysis/
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(msigdbr)
    library(fgsea)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

deg_dir <- file.path(base_path, "protein_differential_analysis_subgroup_adjusted")

output_dir <- file.path(deg_dir, "protein_DEG_meta_analysis")
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
cat("TP53 protein-Level DEG Cross-Cancer Meta-Analysis\n")
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

#' Stouffer's weighted Z-score meta-analysis for a single gene
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

    # --- Step 1: Load all cancer-type DEG results for this comparison ---
    all_deg <- list()
    all_n <- c()
    included_cancers <- c()

    for (ct in cancer_types) {
        deg_file <- file.path(deg_dir, ct, paste0("DEG_", comp_name, ".csv"))
        if (!file.exists(deg_file)) {
            cat("  [SKIP]", ct, ": DEG file not found\n")
            next
        }

        deg <- tryCatch(fread(deg_file), error = function(e) NULL)
        if (is.null(deg) || nrow(deg) == 0) {
            cat("  [SKIP]", ct, ": empty or unreadable\n")
            next
        }

        # Verify required columns
        if (!all(c("gene", "t", "P.Value") %in% colnames(deg))) {
            cat("  [SKIP]", ct, ": missing required columns (gene, t, P.Value)\n")
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

        # Keep gene, gene_symbol, logFC, z
        keep_cols <- intersect(c("gene_symbol", "gene", "logFC", "t", "P.Value", "adj.P.Val", "z"), colnames(deg))
        deg <- deg[, ..keep_cols]

        all_deg[[ct]] <- deg
        all_n[ct] <- n_val
        included_cancers <- c(included_cancers, ct)

        cat("  [OK]", ct, ": ", nrow(deg), "genes, n =", n_val,
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
            n_genes_tested = NA,
            n_sig_005 = NA,
            n_up = NA,
            n_down = NA,
            stringsAsFactors = FALSE
        )
        next
    }

    # --- Step 2: Identify the union of all genes across included cancers ---
    all_genes <- unique(unlist(lapply(all_deg, function(d) d$gene)))
    cat("  Total unique genes across studies:", length(all_genes), "\n")

    # --- Step 3: Build gene_symbol lookup (take the first non-empty mapping) ---
    gene_symbol_lookup <- setNames(rep(NA_character_, length(all_genes)), all_genes)
    for (ct in included_cancers) {
        d <- all_deg[[ct]]
        if ("gene_symbol" %in% colnames(d)) {
            mapped <- d$gene_symbol[match(all_genes, d$gene)]
            fill <- !is.na(mapped) & mapped != "" & is.na(gene_symbol_lookup)
            gene_symbol_lookup[fill] <- mapped[fill]
        }
    }

    # --- Step 4: Stouffer's weighted Z for each gene ---
    cat("  Running Stouffer's weighted Z-score meta-analysis...\n")

    meta_results <- data.frame(
        gene_symbol = character(length(all_genes)),
        gene = all_genes,
        n_studies = integer(length(all_genes)),
        total_n = integer(length(all_genes)),
        Z_meta = numeric(length(all_genes)),
        p_meta = numeric(length(all_genes)),
        direction = character(length(all_genes)),
        stringsAsFactors = FALSE
    )

    # Also store per-cancer z-scores and logFC for the output
    per_cancer_z <- matrix(NA_real_, nrow = length(all_genes), ncol = n_included,
                           dimnames = list(all_genes, included_cancers))
    per_cancer_logFC <- matrix(NA_real_, nrow = length(all_genes), ncol = n_included,
                               dimnames = list(all_genes, included_cancers))

    # Fill per-cancer matrices
    for (ct in included_cancers) {
        d <- all_deg[[ct]]
        idx <- match(d$gene, all_genes)
        valid <- !is.na(idx)
        per_cancer_z[idx[valid], ct] <- d$z[valid]
        per_cancer_logFC[idx[valid], ct] <- d$logFC[valid]
    }

    # Compute meta Z for each gene
    n_vec <- all_n[included_cancers]

    for (g in seq_along(all_genes)) {
        z_row <- per_cancer_z[g, ]
        res <- stouffer_weighted_z(z_row, n_vec)
        meta_results$gene_symbol[g] <- ifelse(is.na(gene_symbol_lookup[all_genes[g]]),
                                              "", gene_symbol_lookup[all_genes[g]])
        meta_results$n_studies[g] <- res$n_studies
        meta_results$total_n[g] <- res$total_n
        meta_results$Z_meta[g] <- res$Z_meta
        meta_results$p_meta[g] <- res$p_meta
        meta_results$direction[g] <- res$direction
    }

    # --- Step 5: FDR correction ---
    meta_results$padj <- p.adjust(meta_results$p_meta, method = "BH")

    # Compute weighted mean logFC (weight = sqrt(n))
    w_mat <- matrix(sqrt(n_vec), nrow = length(all_genes), ncol = n_included, byrow = TRUE)
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

    cat("  Genes tested:", nrow(meta_results), "\n")
    cat("  Significant (padj < 0.05):", n_sig, "\n")
    cat("    Up-regulated:", n_up, "\n")
    cat("    Down-regulated:", n_down, "\n")
    cat("  Included cancers:", paste(included_cancers, collapse = ", "), "\n")

    # --- Step 7: Save output ---
    out_csv <- file.path(output_dir, paste0("META_DEG_", comp_name, ".csv"))
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

    out_xlsx <- file.path(output_dir, paste0("META_DEG_", comp_name, ".xlsx"))
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    cat("  Saved:", basename(out_csv), "&", basename(out_xlsx), "\n\n")

    # Store summary
    meta_summary[[comp_name]] <- data.frame(
        comparison = comp_name,
        n_cancers_included = n_included,
        cancers_included = paste(included_cancers, collapse = ";"),
        n_genes_tested = nrow(meta_results),
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
fwrite(summary_df, file.path(output_dir, "META_DEG_summary.csv"))

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
saveWorkbook(wb_sum, file.path(output_dir, "META_DEG_summary.xlsx"), overwrite = TRUE)

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
cat("DEG Meta-Analysis Complete. Proceeding to Meta-GSEA...\n")
cat("====================================================================\n")

# ==============================================================================
# Section 7: Load Gene Sets from MSigDB (for Meta-GSEA)
# ==============================================================================
# Methodology reference:
#   - Shen K, Tseng GC. Meta-analysis for pathway enrichment analysis when
#     combining multiple genomic studies. Bioinformatics. 2010;26(10):1316-23.
#     PMID: 20410053
#   This approach corresponds to MAPE_G (gene-level meta-analysis followed by
#   GSEA), where gene-level statistics are first combined across studies, and
#   the resulting meta-statistics are used as ranking for preranked GSEA.

cat("\n====================================================================\n")
cat("Meta-GSEA: Preranked GSEA on Cross-Cancer Z_meta Scores\n")
cat("====================================================================\n\n")

cat("Loading MSigDB gene sets...\n")

# Build a global mapping from ensembl_gene to gene_symbol
cat("Building Ensembl to Gene Symbol mapping...\n")
all_genes_df <- suppressMessages(msigdbr(species = "Homo sapiens"))
ens2sym_df <- all_genes_df[!is.na(all_genes_df$ensembl_gene) & all_genes_df$ensembl_gene != "", ]
ens2sym_df <- ens2sym_df[!duplicated(ens2sym_df$ensembl_gene), ]
ens2sym_dict <- setNames(ens2sym_df$gene_symbol, ens2sym_df$ensembl_gene)
rm(all_genes_df, ens2sym_df) # free memory

# Helper: convert msigdbr output to named list for fgsea
msigdbr_to_list <- function(df) {
    df <- df[!is.na(df$ensembl_gene) & df$ensembl_gene != "", ]
    split(df$ensembl_gene, df$gs_name)
}

collections_to_fetch <- list(
    list(name="Hallmark",               cat="H",  sub=NULL),
    list(name="C2_CGP",                 cat="C2", sub="CGP"),
    list(name="C2_CP_BIOCARTA",         cat="C2", sub="CP:BIOCARTA"),
    list(name="C2_CP_KEGG_LEGACY",      cat="C2", sub="CP:KEGG_LEGACY"),
    list(name="C2_CP_KEGG_MEDICUS",     cat="C2", sub="CP:KEGG_MEDICUS"),
    list(name="C2_CP_PID",              cat="C2", sub="CP:PID"),
    list(name="C2_CP_REACTOME",         cat="C2", sub="CP:REACTOME"),
    list(name="C2_CP_WIKIPATHWAYS",     cat="C2", sub="CP:WIKIPATHWAYS"),
    list(name="C3_MIR_MIRDB",           cat="C3", sub="MIR:MIRDB"),
    list(name="C3_MIR_MIR_LEGACY",      cat="C3", sub="MIR:MIR_LEGACY"),
    list(name="C3_TFT_GTRD",            cat="C3", sub="TFT:GTRD"),
    list(name="C3_TFT_TFT_LEGACY",      cat="C3", sub="TFT:TFT_LEGACY"),
    list(name="C4_3CA",                 cat="C4", sub="3CA"),
    list(name="C4_CGN",                 cat="C4", sub="CGN"),
    list(name="C4_CM",                  cat="C4", sub="CM"),
    list(name="C5_GO_BP",               cat="C5", sub="GO:BP"),
    list(name="C5_GO_CC",               cat="C5", sub="GO:CC"),
    list(name="C5_GO_MF",               cat="C5", sub="GO:MF"),
    list(name="C6_Oncogenic",           cat="C6", sub=NULL)#,
    #list(name="C9_CellType",            cat="C9", sub=NULL)
)

collections <- list()

for (info in collections_to_fetch) {
    if (is.null(info$sub)) {
        df <- suppressMessages(msigdbr(species = "Homo sapiens", category = info$cat))
    } else {
        df <- suppressMessages(msigdbr(species = "Homo sapiens", category = info$cat, subcategory = info$sub))
    }
    df <- df %>% dplyr::select(gs_name, ensembl_gene)
    collections[[info$name]] <- msigdbr_to_list(df)
    cat(sprintf("  %-22s : %d gene sets\n", info$name, length(collections[[info$name]])))
}
cat("\n")

# fgsea parameters (same as per-cancer GSEA pipeline)
gsea_minSize <- 15
gsea_maxSize <- 500

# ==============================================================================
# Section 8: Meta-GSEA Function
# ==============================================================================

#' Run fgsea preranked GSEA from a META_DEG file using Z_meta as ranking
#' @param meta_file Path to META_DEG CSV file
#' @param pathways Named list of gene sets (Ensembl IDs)
#' @param mapping_dict Named vector mapping Ensembl IDs to Gene Symbols
#' @return fgsea result data.frame or NULL
run_meta_gsea <- function(meta_file, pathways, mapping_dict = ens2sym_dict) {
    if (!file.exists(meta_file)) return(NULL)

    meta <- fread(meta_file)
    if (!("gene" %in% colnames(meta)) || !("Z_meta" %in% colnames(meta))) {
        cat("    [WARN] Missing gene/Z_meta columns in:", basename(meta_file), "\n")
        return(NULL)
    }

    # Strip Ensembl version suffix (e.g., ENSG00000112742.10 -> ENSG00000112742)
    meta$gene <- sub("\\.[0-9]+$", "", meta$gene)

    # Build named vector of Z_meta statistics
    stats_vec <- setNames(meta$Z_meta, meta$gene)

    # Remove duplicates (keep first = most significant by padj ordering)
    stats_vec <- stats_vec[!duplicated(names(stats_vec))]

    # Remove non-finite values
    stats_vec <- stats_vec[is.finite(stats_vec)]

    if (length(stats_vec) < 100) {
        cat("    [WARN] Too few genes:", length(stats_vec), "\n")
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
    # And map Ensembl IDs to Gene Symbols
    res$leadingEdge <- vapply(res$leadingEdge, function(x) {
        symbols <- mapping_dict[x]
        symbols[is.na(symbols) | symbols == ""] <- x[is.na(symbols) | symbols == ""]
        paste(symbols, collapse = ";")
    }, character(1))

    # Sort by padj
    res <- res[order(res$padj), ]
    as.data.frame(res)
}

# ==============================================================================
# Section 9: Run Meta-GSEA for Each Comparison
# ==============================================================================

gsea_output_dir <- file.path(output_dir, "protein_meta_GSEA")
dir.create(gsea_output_dir, recursive = TRUE, showWarnings = FALSE)

# Storage for GSEA summary statistics (one per collection)
all_gsea_summaries <- list()
for (coll_name in names(collections)) {
    all_gsea_summaries[[coll_name]] <- list()
}

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Meta-GSEA for:", comp_name, "\n")
    cat("==================================================================\n")

    meta_file <- file.path(output_dir, paste0("META_DEG_", comp_name, ".csv"))
    if (!file.exists(meta_file)) {
        cat("  [SKIP] META_DEG file not found:", basename(meta_file), "\n\n")
        for (coll_name in names(collections)) {
            all_gsea_summaries[[coll_name]][[comp_name]] <- data.frame(
                comparison = comp_name,
                stringsAsFactors = FALSE
            )
        }
        next
    }

    for (coll_name in names(collections)) {
        pathways <- collections[[coll_name]]
        cat("\n  --- Collection:", coll_name, "---\n")

        # Create collection output subdirectory
        coll_out <- file.path(gsea_output_dir, coll_name)
        dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)

        # Summary row for this comparison + collection
        summary_row <- list(comparison = comp_name)

        cat("    ", comp_label, ":", sep = "")

        res <- run_meta_gsea(meta_file, pathways)

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

            header_style <- createStyle(
                textDecoration = "bold", halign = "center",
                border = "bottom", fgFill = "#4472C4", fontColour = "white"
            )
            addStyle(wb, 1, header_style,
                     rows = 1, cols = 1:ncol(res), gridExpand = TRUE)
            setColWidths(wb, 1, cols = 1:ncol(res), widths = "auto")

            saveWorkbook(wb, file.path(coll_out, paste0("META_GSEA_", comp_name, ".xlsx")),
                         overwrite = TRUE)

            summary_row[[paste0(comp_name, "_total")]] <- nrow(res)
            summary_row[[paste0(comp_name, "_sig")]] <- n_sig
            summary_row[[paste0(comp_name, "_up")]] <- n_up
            summary_row[[paste0(comp_name, "_down")]] <- n_down
        } else {
            cat(" no META_DEG file or insufficient data\n")
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
cat("Full Pipeline Complete (DEG Meta-Analysis + Meta-GSEA)!\n")
cat("DEG Meta-Analysis results saved to:", output_dir, "\n")
cat("Meta-GSEA results saved to:", gsea_output_dir, "\n")
cat("====================================================================\n")

################################################################################
# Phos without PN Gene-Level GSEA Meta vs Phos with PN Gene-Level GSEA Meta
# Integration
#
# Purpose:
#   Integrate phosphoprotein gene-level GSEA meta-analysis results from
#   without-protein-normalization and with-protein-normalization pipelines
#   by merging on the "pathway" column for each aggregation method, collection,
#   and comparison pair.
#
# Input:
#   - phosphoprotein_without_PN_gene_level_GSEA_meta/{agg_method}/{collection}/META_GSEA_{comparison}.csv
#     (from phosphoprotein_without_PN_gene_level_GSEA_meta_analysis.R)
#   - phosphoprotein_with_PN_gene_level_GSEA_meta/{agg_method}/{collection}/META_GSEA_{comparison}.csv
#     (from phosphoprotein_with_PN_gene_level_GSEA_meta_analysis.R)
#
# Output:
#   - phos_without_PN_gene_level_GSEA_meta_vs_phos_with_PN_gene_level_GSEA_meta/
#       {agg_method}/{collection}/META_GSEA_{comparison}.csv and .xlsx
#       {agg_method}/integration_summary.csv and .xlsx
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(ggrepel)
    library(msigdbr)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Input directories (base; agg_method subdirectory will be appended)
without_PN_meta_base <- file.path(base_path, "phosphoprotein_without_PN_gene_level_GSEA_meta")
with_PN_meta_base    <- file.path(base_path, "phosphoprotein_with_PN_gene_level_GSEA_meta")

# Output directory (base)
output_base <- file.path(base_path, "phos_without_PN_gene_level_GSEA_meta_vs_phos_with_PN_gene_level_GSEA_meta")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# Two aggregation methods
agg_methods <- c("max_abs_t", "median_t")

# 9 comparisons
comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

# 19 MSigDB collections (matching the gene-level GSEA meta-analysis pipelines)
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

cat("====================================================================\n")
cat("Phos without PN Gene-Level GSEA Meta vs\n")
cat("Phos with PN Gene-Level GSEA Meta Integration\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Integration Loop
#   Outer loop: aggregation method (max_abs_t, median_t)
#   Inner loops: collection x comparison
# ==============================================================================

for (agg_method in agg_methods) {

    cat("####################################################################\n")
    cat("Aggregation method:", agg_method, "\n")
    cat("####################################################################\n\n")

    # Input directories for this aggregation method
    without_PN_meta_dir <- file.path(without_PN_meta_base, agg_method)
    with_PN_meta_dir    <- file.path(with_PN_meta_base, agg_method)

    # Output directory for this aggregation method
    output_dir <- file.path(output_base, agg_method)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    integration_summary <- list()

    for (coll_name in collection_names) {
        cat("==================================================================\n")
        cat("Collection:", coll_name, " [", agg_method, "]\n")
        cat("==================================================================\n")

        # Check if collection directories exist in at least one input
        without_PN_coll_dir <- file.path(without_PN_meta_dir, coll_name)
        with_PN_coll_dir    <- file.path(with_PN_meta_dir, coll_name)

        if (!dir.exists(without_PN_coll_dir) && !dir.exists(with_PN_coll_dir)) {
            cat("  [SKIP] Collection directory not found in either input\n\n")
            next
        }

        coll_out_dir <- file.path(output_dir, coll_name)
        dir.create(coll_out_dir, recursive = TRUE, showWarnings = FALSE)

        for (j in seq_len(nrow(comparisons))) {
            comp_name  <- comparisons$name[j]
            comp_label <- comparisons$label[j]

            cat("  Integration for:", comp_name, "\n")

            without_PN_file <- file.path(without_PN_coll_dir, paste0("META_GSEA_", comp_name, ".csv"))
            with_PN_file    <- file.path(with_PN_coll_dir, paste0("META_GSEA_", comp_name, ".csv"))

            # Check if either file exists
            if (!file.exists(without_PN_file) && !file.exists(with_PN_file)) {
                cat("    [SKIP] Missing both input files\n")
                integration_summary[[paste(coll_name, comp_name, sep = "_")]] <- data.frame(
                    collection = coll_name,
                    comparison = comp_name,
                    n_phos_gl_GSEA_without_PN_pathways = NA,
                    n_phos_gl_GSEA_with_PN_pathways = NA,
                    n_merged_pathways = NA,
                    n_both_sig_005 = NA,
                    stringsAsFactors = FALSE
                )
                next
            }

            # --- Load phos_without_PN gene-level GSEA meta data ---
            without_PN_data <- data.frame(pathway = character())
            if (file.exists(without_PN_file)) {
                tmp_without_PN <- fread(without_PN_file)
                if (nrow(tmp_without_PN) > 0 && "pathway" %in% colnames(tmp_without_PN)) {
                    # Rename all columns except "pathway" with _phos_gl_GSEA_without_PN suffix
                    cols_to_rename <- colnames(tmp_without_PN)[colnames(tmp_without_PN) != "pathway"]
                    setnames(tmp_without_PN, cols_to_rename, paste0(cols_to_rename, "_phos_gl_GSEA_without_PN"))
                    # Remove duplicate pathways
                    tmp_without_PN <- tmp_without_PN[!duplicated(tmp_without_PN$pathway), ]
                    without_PN_data <- tmp_without_PN
                }
            }

            # --- Load phos_with_PN gene-level GSEA meta data ---
            with_PN_data <- data.frame(pathway = character())
            if (file.exists(with_PN_file)) {
                tmp_with_PN <- fread(with_PN_file)
                if (nrow(tmp_with_PN) > 0 && "pathway" %in% colnames(tmp_with_PN)) {
                    # Rename all columns except "pathway" with _phos_gl_GSEA_with_PN suffix
                    cols_to_rename <- colnames(tmp_with_PN)[colnames(tmp_with_PN) != "pathway"]
                    setnames(tmp_with_PN, cols_to_rename, paste0(cols_to_rename, "_phos_gl_GSEA_with_PN"))
                    # Remove duplicate pathways
                    tmp_with_PN <- tmp_with_PN[!duplicated(tmp_with_PN$pathway), ]
                    with_PN_data <- tmp_with_PN
                }
            }

            if (nrow(without_PN_data) == 0 && nrow(with_PN_data) == 0) {
                cat("    [SKIP] No data in either file\n")
                integration_summary[[paste(coll_name, comp_name, sep = "_")]] <- data.frame(
                    collection = coll_name,
                    comparison = comp_name,
                    n_phos_gl_GSEA_without_PN_pathways = 0,
                    n_phos_gl_GSEA_with_PN_pathways = 0,
                    n_merged_pathways = 0,
                    n_both_sig_005 = NA,
                    stringsAsFactors = FALSE
                )
                next
            }

            # --- Merge by pathway (full join to keep all pathways) ---
            merged <- merge(without_PN_data, with_PN_data, by = "pathway", all = TRUE)

            # --- Compute per-dataset significance counts ---
            all_cols <- colnames(merged)

            nes_without_PN_cols  <- grep("^NES_.*_phos_gl_GSEA_without_PN$",  all_cols, value = TRUE)
            pval_without_PN_cols <- grep("^pval_.*_phos_gl_GSEA_without_PN$", all_cols, value = TRUE)
            nes_with_PN_cols     <- grep("^NES_.*_phos_gl_GSEA_with_PN$",  all_cols, value = TRUE)
            pval_with_PN_cols    <- grep("^pval_.*_phos_gl_GSEA_with_PN$", all_cols, value = TRUE)

            # Extract cancer type names from NES columns to ensure matched pairing
            without_PN_cancers <- sub("^NES_(.*)_phos_gl_GSEA_without_PN$", "\\1", nes_without_PN_cols)
            with_PN_cancers    <- sub("^NES_(.*)_phos_gl_GSEA_with_PN$", "\\1", nes_with_PN_cols)

            # phos_without_PN: count datasets with significant up (NES > 0 & pval < 0.05) and down
            merged$without_PN_up_datasets_num   <- 0L
            merged$without_PN_down_datasets_num <- 0L
            for (ct in without_PN_cancers) {
                nes_col  <- paste0("NES_", ct, "_phos_gl_GSEA_without_PN")
                pval_col <- paste0("pval_", ct, "_phos_gl_GSEA_without_PN")
                if (nes_col %in% all_cols && pval_col %in% all_cols) {
                    nes_vals  <- merged[[nes_col]]
                    pval_vals <- merged[[pval_col]]
                    is_up   <- !is.na(nes_vals) & !is.na(pval_vals) & nes_vals > 0 & pval_vals < 0.05
                    is_down <- !is.na(nes_vals) & !is.na(pval_vals) & nes_vals < 0 & pval_vals < 0.05
                    merged$without_PN_up_datasets_num   <- merged$without_PN_up_datasets_num + as.integer(is_up)
                    merged$without_PN_down_datasets_num <- merged$without_PN_down_datasets_num + as.integer(is_down)
                }
            }

            # phos_with_PN: count datasets with significant up and down
            merged$with_PN_up_datasets_num   <- 0L
            merged$with_PN_down_datasets_num <- 0L
            for (ct in with_PN_cancers) {
                nes_col  <- paste0("NES_", ct, "_phos_gl_GSEA_with_PN")
                pval_col <- paste0("pval_", ct, "_phos_gl_GSEA_with_PN")
                if (nes_col %in% all_cols && pval_col %in% all_cols) {
                    nes_vals  <- merged[[nes_col]]
                    pval_vals <- merged[[pval_col]]
                    is_up   <- !is.na(nes_vals) & !is.na(pval_vals) & nes_vals > 0 & pval_vals < 0.05
                    is_down <- !is.na(nes_vals) & !is.na(pval_vals) & nes_vals < 0 & pval_vals < 0.05
                    merged$with_PN_up_datasets_num   <- merged$with_PN_up_datasets_num + as.integer(is_up)
                    merged$with_PN_down_datasets_num <- merged$with_PN_down_datasets_num + as.integer(is_down)
                }
            }

            # --- Collect per-dataset significant dataset names ---
            # phos_without_PN: names of significantly up/down datasets
            without_PN_up_names_list   <- vector("list", nrow(merged))
            without_PN_down_names_list <- vector("list", nrow(merged))
            for (i in seq_len(nrow(merged))) {
                up_names   <- character()
                down_names <- character()
                for (ct in without_PN_cancers) {
                    nes_val  <- merged[[paste0("NES_", ct, "_phos_gl_GSEA_without_PN")]][i]
                    pval_val <- merged[[paste0("pval_", ct, "_phos_gl_GSEA_without_PN")]][i]
                    if (!is.na(nes_val) && !is.na(pval_val) && pval_val < 0.05) {
                        if (nes_val > 0) up_names   <- c(up_names, ct)
                        if (nes_val < 0) down_names <- c(down_names, ct)
                    }
                }
                without_PN_up_names_list[[i]]   <- paste(up_names, collapse = "; ")
                without_PN_down_names_list[[i]] <- paste(down_names, collapse = "; ")
            }
            merged$without_PN_up_datasets_name   <- unlist(without_PN_up_names_list)
            merged$without_PN_down_datasets_name <- unlist(without_PN_down_names_list)

            # phos_with_PN: names of significantly up/down datasets
            with_PN_up_names_list   <- vector("list", nrow(merged))
            with_PN_down_names_list <- vector("list", nrow(merged))
            for (i in seq_len(nrow(merged))) {
                up_names   <- character()
                down_names <- character()
                for (ct in with_PN_cancers) {
                    nes_val  <- merged[[paste0("NES_", ct, "_phos_gl_GSEA_with_PN")]][i]
                    pval_val <- merged[[paste0("pval_", ct, "_phos_gl_GSEA_with_PN")]][i]
                    if (!is.na(nes_val) && !is.na(pval_val) && pval_val < 0.05) {
                        if (nes_val > 0) up_names   <- c(up_names, ct)
                        if (nes_val < 0) down_names <- c(down_names, ct)
                    }
                }
                with_PN_up_names_list[[i]]   <- paste(up_names, collapse = "; ")
                with_PN_down_names_list[[i]] <- paste(down_names, collapse = "; ")
            }
            merged$with_PN_up_datasets_name   <- unlist(with_PN_up_names_list)
            merged$with_PN_down_datasets_name <- unlist(with_PN_down_names_list)

            # Reorder columns: priority columns first, then count/name columns, then rest
            priority_cols <- c("pathway", "Z_meta_phos_gl_GSEA_with_PN", "padj_phos_gl_GSEA_with_PN",
                               "Z_meta_phos_gl_GSEA_without_PN", "padj_phos_gl_GSEA_without_PN")
            count_name_cols <- c("without_PN_up_datasets_num", "without_PN_down_datasets_num",
                                 "without_PN_up_datasets_name", "without_PN_down_datasets_name",
                                 "with_PN_up_datasets_num", "with_PN_down_datasets_num",
                                 "with_PN_up_datasets_name", "with_PN_down_datasets_name")
            # Keep only columns that actually exist in merged
            priority_cols   <- priority_cols[priority_cols %in% colnames(merged)]
            count_name_cols <- count_name_cols[count_name_cols %in% colnames(merged)]
            other_cols <- setdiff(colnames(merged), c(priority_cols, count_name_cols))
            setcolorder(merged, c(priority_cols, count_name_cols, other_cols))

            # Compute n_both_sig (both padj < 0.05)
            n_both_sig <- NA
            if ("padj_phos_gl_GSEA_without_PN" %in% colnames(merged) && "padj_phos_gl_GSEA_with_PN" %in% colnames(merged)) {
                n_both_sig <- sum(
                    merged$padj_phos_gl_GSEA_without_PN < 0.05 & merged$padj_phos_gl_GSEA_with_PN < 0.05,
                    na.rm = TRUE
                )
            }

            # Sort by Z_meta_phos_gl_GSEA_with_PN descending
            if ("Z_meta_phos_gl_GSEA_with_PN" %in% colnames(merged)) {
                merged <- merged[order(-merged$Z_meta_phos_gl_GSEA_with_PN, na.last = TRUE), ]
            }

            n_without_PN <- sum(!is.na(without_PN_data$pathway) & without_PN_data$pathway != "")
            n_with_PN    <- sum(!is.na(with_PN_data$pathway) & with_PN_data$pathway != "")

            cat("    Phos gene-level GSEA without PN pathways:", n_without_PN, "\n")
            cat("    Phos gene-level GSEA with PN pathways:", n_with_PN, "\n")
            cat("    Merged pathways:", nrow(merged), "\n")
            if (!is.na(n_both_sig)) {
                cat("    Both significant (padj < 0.05):", n_both_sig, "\n")
            }

            # --- Save CSV ---
            out_csv <- file.path(coll_out_dir, paste0("META_GSEA_", comp_name, ".csv"))
            fwrite(merged, out_csv)

            # --- Save XLSX ---
            wb <- createWorkbook()
            addWorksheet(wb, comp_label)
            writeData(wb, 1, merged)

            nc <- ncol(merged)
            nr <- nrow(merged)

            # Header style: Arial 7, bold, blue background, white text, centered
            header_style <- createStyle(
                fontName = "Arial", fontSize = 7,
                textDecoration = "bold", halign = "center", valign = "center",
                border = "bottom", fgFill = "#4472C4", fontColour = "white"
            )
            addStyle(wb, 1, header_style,
                     rows = 1, cols = 1:nc, gridExpand = TRUE)

            # Body style: Arial 7, left-aligned, vertically centered
            body_style <- createStyle(
                fontName = "Arial", fontSize = 7,
                halign = "left", valign = "center"
            )
            if (nr > 0) {
                addStyle(wb, 1, body_style,
                         rows = 2:(nr + 1), cols = 1:nc, gridExpand = TRUE, stack = FALSE)
            }

            setColWidths(wb, 1, cols = 1:nc, widths = "auto")

            # Color rows using the same logic as scatter plots
            has_both <- "Z_meta_phos_gl_GSEA_without_PN" %in% colnames(merged) &
                        "Z_meta_phos_gl_GSEA_with_PN" %in% colnames(merged) &
                        "padj_phos_gl_GSEA_without_PN" %in% colnames(merged) &
                        "padj_phos_gl_GSEA_with_PN" %in% colnames(merged)

            if (has_both) {
                z_wo  <- merged$Z_meta_phos_gl_GSEA_without_PN
                z_w   <- merged$Z_meta_phos_gl_GSEA_with_PN
                pa_wo <- merged$padj_phos_gl_GSEA_without_PN
                pa_w  <- merged$padj_phos_gl_GSEA_with_PN

                # Red: both positive & both significant
                red_rows <- which(z_wo > 0 & z_w > 0 & pa_wo < 0.05 & pa_w < 0.05) + 1
                # Blue: both negative & both significant
                blue_rows <- which(z_wo < 0 & z_w < 0 & pa_wo < 0.05 & pa_w < 0.05) + 1
                # Orange: with_PN positive & with_PN sig, but NOT (without_PN positive & without_PN sig)
                orange_rows <- which(!(z_wo > 0 & pa_wo < 0.05) & z_w > 0 & pa_w < 0.05) + 1
                # Green: with_PN negative & with_PN sig, but NOT (without_PN negative & without_PN sig)
                green_rows <- which(!(z_wo < 0 & pa_wo < 0.05) & z_w < 0 & pa_w < 0.05) + 1

                nc <- ncol(merged)
                if (length(red_rows) > 0) {
                    addStyle(wb, 1, createStyle(fgFill = "#FFCCCC"),
                             rows = red_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
                }
                if (length(blue_rows) > 0) {
                    addStyle(wb, 1, createStyle(fgFill = "#CCE5FF"),
                             rows = blue_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
                }
                if (length(orange_rows) > 0) {
                    addStyle(wb, 1, createStyle(fgFill = "#FFE0B2"),
                             rows = orange_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
                }
                if (length(green_rows) > 0) {
                    addStyle(wb, 1, createStyle(fgFill = "#C8E6C9"),
                             rows = green_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
                }
            }

            out_xlsx <- file.path(coll_out_dir, paste0("META_GSEA_", comp_name, ".xlsx"))
            saveWorkbook(wb, out_xlsx, overwrite = TRUE)

            cat("    Saved:", basename(out_csv), "&", basename(out_xlsx), "\n")

            # Store summary
            integration_summary[[paste(coll_name, comp_name, sep = "_")]] <- data.frame(
                collection = coll_name,
                comparison = comp_name,
                n_phos_gl_GSEA_without_PN_pathways = n_without_PN,
                n_phos_gl_GSEA_with_PN_pathways = n_with_PN,
                n_merged_pathways = nrow(merged),
                n_both_sig_005 = n_both_sig,
                stringsAsFactors = FALSE
            )
        }
        cat("\n")
    }

    # ==========================================================================
    # Section 3: Integration Summary (per agg_method)
    # ==========================================================================

    cat("====================================================================\n")
    cat("Generating Integration Summary [", agg_method, "]\n")
    cat("====================================================================\n\n")

    summary_df <- bind_rows(integration_summary)

    if (nrow(summary_df) > 0) {
        # CSV
        fwrite(summary_df, file.path(output_dir, "integration_summary.csv"))

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
        saveWorkbook(wb_sum, file.path(output_dir, "integration_summary.xlsx"), overwrite = TRUE)

        # Print summary table
        cat("=== Integration Summary [", agg_method, "] ===\n\n")
        print(as.data.frame(summary_df), row.names = FALSE)
        cat("\n")
    } else {
        cat("No data merged.\n")
    }

    cat("====================================================================\n")
    cat("Integration Complete! [", agg_method, "]\n")
    cat("Results saved to:", output_dir, "\n")
    cat("====================================================================\n\n")
}

# ==============================================================================
# Section 4: Scatter Plots for Gene-Level GSEA Meta Integration
# ==============================================================================

cat("\n====================================================================\n")
cat("Generating Scatter Plots\n")
cat("====================================================================\n\n")

# Build pathway -> gene set size lookup from MSigDB
cat("Loading MSigDB gene set sizes...\n")
all_gs <- suppressMessages(msigdbr(species = "Homo sapiens"))
pathway_size_df <- all_gs %>%
    group_by(gs_name) %>%
    summarise(gene_set_size = n_distinct(gene_symbol), .groups = "drop")
pathway_size_lookup <- setNames(pathway_size_df$gene_set_size, pathway_size_df$gs_name)
rm(all_gs, pathway_size_df)
cat("  Loaded sizes for", length(pathway_size_lookup), "gene sets\n\n")

# User can modify this vector to select which collections to plot
collections_to_plot <- c("Hallmark")
# "Hallmark",
# "C2_CGP", "C2_CP_BIOCARTA", "C2_CP_KEGG_LEGACY", "C2_CP_KEGG_MEDICUS",
# "C2_CP_PID", "C2_CP_REACTOME", "C2_CP_WIKIPATHWAYS",
# "C3_MIR_MIRDB", "C3_MIR_MIR_LEGACY",
# "C3_TFT_GTRD", "C3_TFT_TFT_LEGACY",
# "C4_3CA", "C4_CGN", "C4_CM",
# "C5_GO_BP", "C5_GO_CC", "C5_GO_MF",
# "C6_Oncogenic"

for (agg_method in agg_methods) {

    cat("####################################################################\n")
    cat("Scatter plots for aggregation method:", agg_method, "\n")
    cat("####################################################################\n\n")

    output_dir <- file.path(output_base, agg_method)
    plot_dir <- file.path(output_dir, "scatter_plot")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    for (coll_name in collections_to_plot) {
        cat("Plotting collection:", coll_name, " [", agg_method, "]\n")

        for (j in seq_len(nrow(comparisons))) {
            comp_name <- comparisons$name[j]

            int_file <- file.path(output_dir, coll_name, paste0("META_GSEA_", comp_name, ".csv"))
            if (!file.exists(int_file)) {
                next
            }

            df <- fread(int_file)

            # Keep only rows where both Z_meta are finite
            df <- df[is.finite(Z_meta_phos_gl_GSEA_without_PN) & is.finite(Z_meta_phos_gl_GSEA_with_PN)]

            if (nrow(df) < 5) {
                next
            }

            # Assign colors based on rules
            df$fill_color <- "medium_grey" # Rule 3
            df$fill_color[df$padj_phos_gl_GSEA_without_PN > 0.05 & df$padj_phos_gl_GSEA_with_PN > 0.05] <- "light_grey" # Rule 2

            # Identify highlighted pathways (Rule 1)
            # Red: both Z_meta > 0 & both padj < 0.05
            path_red <- df$pathway[df$Z_meta_phos_gl_GSEA_without_PN > 0 & df$Z_meta_phos_gl_GSEA_with_PN > 0 & df$padj_phos_gl_GSEA_without_PN < 0.05 & df$padj_phos_gl_GSEA_with_PN < 0.05]

            # Blue: both Z_meta < 0 & both padj < 0.05
            path_blue <- df$pathway[df$Z_meta_phos_gl_GSEA_without_PN < 0 & df$Z_meta_phos_gl_GSEA_with_PN < 0 & df$padj_phos_gl_GSEA_without_PN < 0.05 & df$padj_phos_gl_GSEA_with_PN < 0.05]

            # Orange: with_PN positive & sig, excluding (without_PN positive & sig)
            path_orange <- df$pathway[!(df$Z_meta_phos_gl_GSEA_without_PN > 0 & df$padj_phos_gl_GSEA_without_PN < 0.05) & df$Z_meta_phos_gl_GSEA_with_PN > 0 & df$padj_phos_gl_GSEA_with_PN < 0.05]

            # Green: with_PN negative & sig, excluding (without_PN negative & sig)
            path_green <- df$pathway[!(df$Z_meta_phos_gl_GSEA_without_PN < 0 & df$padj_phos_gl_GSEA_without_PN < 0.05) & df$Z_meta_phos_gl_GSEA_with_PN < 0 & df$padj_phos_gl_GSEA_with_PN < 0.05]

            # Apply Rule 1 colors (highest priority, applied last)
            df$fill_color[df$pathway %in% path_green] <- "green"
            df$fill_color[df$pathway %in% path_orange] <- "orange"
            df$fill_color[df$pathway %in% path_blue]  <- "blue"
            df$fill_color[df$pathway %in% path_red]  <- "red"

            # Prepare data for label filtering
            df$Z_sum <- df$Z_meta_phos_gl_GSEA_without_PN + df$Z_meta_phos_gl_GSEA_with_PN

            # Filter red labels (top 3 Z_sum if > 3)
            if (length(path_red) > 3) {
                red_df <- df[df$pathway %in% path_red, ]
                red_df <- red_df[order(-Z_sum)]
                label_red <- red_df$pathway[1:3]
            } else {
                label_red <- path_red
            }

            # Filter blue labels (bottom 3 Z_sum if > 3)
            if (length(path_blue) > 3) {
                blue_df <- df[df$pathway %in% path_blue, ]
                blue_df <- blue_df[order(Z_sum)]
                label_blue <- blue_df$pathway[1:3]
            } else {
                label_blue <- path_blue
            }

            # Filter orange labels (top 3 Z_meta_with_PN if > 3)
            if (length(path_orange) > 3) {
                orange_df <- df[df$pathway %in% path_orange, ]
                orange_df <- orange_df[order(-Z_meta_phos_gl_GSEA_with_PN)]
                label_orange <- orange_df$pathway[1:3]
            } else {
                label_orange <- path_orange
            }

            # Filter green labels (bottom 3 Z_meta_with_PN if > 3)
            if (length(path_green) > 3) {
                green_df <- df[df$pathway %in% path_green, ]
                green_df <- green_df[order(Z_meta_phos_gl_GSEA_with_PN)]
                label_green <- green_df$pathway[1:3]
            } else {
                label_green <- path_green
            }

            # Rule 1 union for labeling
            rule1_label_paths <- unique(c(label_red, label_blue, label_orange, label_green))
            df$show_label <- df$pathway %in% rule1_label_paths

            # Format labels based on collection
            df$label <- df$pathway
            if (coll_name == "Hallmark") {
                df$label <- gsub("^HALLMARK_", "", df$label)
            } else if (coll_name == "C2_CP_KEGG_LEGACY") {
                df$label <- gsub("^KEGG_", "", df$label)
            }

            # Define color palette
            fill_values <- c(
                "red"         = "#D62728",
                "blue"        = "#1F77B4",
                "orange"      = "#FF7F0E",
                "green"       = "#2CA02C",
                "light_grey"  = "#D9D9D9",
                "medium_grey" = "#A0A0A0"
            )

            # Map gene set size from MSigDB lookup
            df$gene_set_size <- pathway_size_lookup[df$pathway]
            # Fill NA sizes with median to avoid dropping points
            median_size <- median(df$gene_set_size, na.rm = TRUE)
            if (is.na(median_size)) median_size <- 50  # fallback default
            df$gene_set_size[is.na(df$gene_set_size)] <- median_size

            # Split data for layered plotting (grey in background, colored in foreground)
            df_grey <- df[fill_color %in% c("light_grey", "medium_grey")]
            df_colored <- df[fill_color %in% c("red", "blue", "orange", "green")]
            df_label <- df[show_label == TRUE]

            # Let each axis adapt to its own data range (with 5% padding)
            x_range <- range(df$Z_meta_phos_gl_GSEA_with_PN, na.rm = TRUE)
            y_range <- range(df$Z_meta_phos_gl_GSEA_without_PN, na.rm = TRUE)
            x_pad <- diff(x_range) * 0.05
            y_pad <- diff(y_range) * 0.05

            p <- ggplot(df, aes(x = Z_meta_phos_gl_GSEA_with_PN, y = Z_meta_phos_gl_GSEA_without_PN)) +
                # y=x reference line
                geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
                # Background grey points
                geom_point(data = df_grey, aes(color = fill_color, size = gene_set_size), shape = 16, alpha = 0.5) +
                # Foreground colored points
                geom_point(data = df_colored, aes(color = fill_color, size = gene_set_size), shape = 16, alpha = 0.85) +
                # Pathway labels
                geom_text_repel(
                    data = df_label,
                    aes(label = label),
                    size = 2.2,
                    segment.size = 0.25,
                    segment.color = "grey40",
                    segment.alpha = 0.7,
                    min.segment.length = 0,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    force = 40,
                    force_pull = 0.2,
                    max.overlaps = Inf,
                    max.iter = 100000,
                    max.time = 10,
                    seed = 42
                ) +
                scale_color_manual(values = fill_values, guide = "none") +
                scale_size_continuous(range = c(0.5, 3.5), name = "Gene Set Size") +
                coord_cartesian(
                    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
                    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
                    clip = "off"
                ) +
                labs(
                    x = "Phos Gene-Level GSEA with PN Z_meta",
                    y = "Phos Gene-Level GSEA without PN Z_meta",
                    title = paste(gsub("_", " ", comp_name), "-", coll_name, "-", agg_method)
                ) +
                theme_classic(base_size = 10) +
                theme(
                    plot.background  = element_rect(fill = "white", color = NA),
                    panel.background = element_rect(fill = "white", color = NA),
                    panel.border     = element_blank(),
                    axis.line        = element_line(color = "black", linewidth = 0.4),
                    axis.ticks       = element_line(color = "black", linewidth = 0.3),
                    axis.text        = element_text(color = "black", size = 8),
                    axis.title       = element_text(color = "black", size = 10),
                    plot.title       = element_text(hjust = 0.5, face = "bold", size = 11),
                    plot.margin      = margin(15, 40, 15, 15)
                )

            # Save TIFF
            tiff_file <- file.path(plot_dir, paste0("scatter_", coll_name, "_", comp_name, ".tiff"))
            tiff(tiff_file, width = 5.25, height = 4.95, units = "in", res = 300, compression = "lzw")
            print(p)
            dev.off()

            # Save PDF
            pdf_file <- file.path(plot_dir, paste0("scatter_", coll_name, "_", comp_name, ".pdf"))
            pdf(pdf_file, width = 5.25, height = 4.95)
            print(p)
            dev.off()

            cat("    Saved plots for", comp_name, "\n")
        }
    }
}

cat("====================================================================\n")
cat("Scatter plots completed!\n")
cat("====================================================================\n")

cat("\n####################################################################\n")
cat("All aggregation methods processed. Integration fully complete!\n")
cat("####################################################################\n")

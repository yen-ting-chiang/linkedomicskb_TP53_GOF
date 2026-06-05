################################################################################
# Phosphoprotein DPS (with Protein Normalization) vs
# Phosphoprotein DPS (without Protein Normalization) Meta-Analysis Integration
#
# Purpose:
#   Integrate phosphoprotein-level meta-analysis results from two normalization
#   approaches:
#     (1) Phosphoprotein DPS with protein normalization (standard approach)
#     (2) Phosphoprotein DPS without protein normalization
#   by merging on the shared "phosphosite" column (exact match).
#
#   Because both datasets are at the phosphosite level, the merge key is
#   simply the "phosphosite" column (e.g.,
#   ENSG00000076003.5|ENSP00000264156.2|S762|EIESEIDSEEELINK|1).
#   An inner join is performed: only phosphosites present in both datasets
#   are retained.
#
#   Computed derived columns:
#     - Z_meta_sum  : with_PN Z_meta + without_PN Z_meta (concordant regulation)
#     - Z_meta_diff : without_PN Z_meta - with_PN Z_meta (discordant regulation)
#     - padj_max    : max(with_PN padj, without_PN padj) (conservative threshold)
#
# Input:
#   - phosphoprotein_differential_analysis_subgroup_adjusted/
#       phosphoprotein_DPS_meta_analysis/
#       META_DPS_{comparison}.csv  (9 comparisons, excluding meta_GSEA/meta_KSEA
#       subdirectories)
#   - phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/
#       phosphoprotein_DPS_without_protein_normalization_meta_analysis/
#       META_DPS_{comparison}.csv  (9 comparisons, excluding meta_GSEA/meta_KSEA
#       subdirectories)
#
# Output:
#   - phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta/
#       META_integrated_{comparison}.csv and .xlsx  (9 comparisons)
#       integration_summary.csv and .xlsx
#   - phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta/scatter_plot/
#       scatter_{comparison}.tiff and .pdf  (9 comparisons)
#
# Rationale:
#   Z_meta_sum captures concordant regulation across both normalization
#   approaches; phosphosites with large |Z_meta_sum| are consistently altered
#   regardless of normalization. Z_meta_diff captures discordance (differences
#   attributable to the normalization method); phosphosites with large
#   |Z_meta_diff| behave differently depending on whether protein abundance
#   is regressed out. padj_max provides a conservative significance threshold
#   requiring both approaches to yield individually significant results.
#
# References:
#   - Zhang B, et al. Proteogenomic characterization of human colon and
#     rectal cancer. Nature. 2014;513(7518):382-7. PMID: 25043054
#   - Mertins P, et al. Proteogenomics connects somatic mutations to
#     signalling in breast cancer. Nature. 2016;534(7605):55-62.
#     PMID: 27251275
#   - Vasaikar S, et al. Proteogenomic Analysis of Human Colon Cancer
#     Reveals New Therapeutic Vulnerabilities. Cell. 2019;177(4):1035-1049.e19.
#     PMID: 31031003
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(ggrepel)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Input directories (meta-analysis output from the two normalization pipelines)
# (1) Phosphoprotein DPS with protein normalization
phos_with_PN_meta_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_subgroup_adjusted",
    "phosphoprotein_DPS_meta_analysis")

# (2) Phosphoprotein DPS without protein normalization
phos_without_PN_meta_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted",
    "phosphoprotein_DPS_without_protein_normalization_meta_analysis")

# Output directory
output_dir <- file.path(base_path, "phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 9 comparisons (must match the naming convention of META_DPS files)
comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

cat("====================================================================\n")
cat("Phospho DPS (with PN) vs Phospho DPS (w/o PN) Meta Integration\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Integration Loop (per comparison)
# ==============================================================================

# Storage for summary statistics
integration_summary <- list()

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Integration for:", comp_name, "\n")
    cat("==================================================================\n")

    # --- Step 1: Locate input files ---
    with_PN_file <- file.path(phos_with_PN_meta_dir,
                              paste0("META_DPS_", comp_name, ".csv"))
    without_PN_file <- file.path(phos_without_PN_meta_dir,
                                 paste0("META_DPS_", comp_name, ".csv"))

    if (!file.exists(with_PN_file)) {
        cat("  [SKIP] Phospho (with PN) META_DPS file not found:",
            basename(with_PN_file), "\n\n")
        integration_summary[[comp_name]] <- data.frame(
            comparison = comp_name,
            n_phos_with_PN = NA, n_phos_without_PN = NA, n_merged_rows = NA,
            n_unique_phosphosites_merged = NA,
            n_both_sig_005 = NA, n_concordant_sig = NA, n_discordant_sig = NA,
            stringsAsFactors = FALSE
        )
        next
    }
    if (!file.exists(without_PN_file)) {
        cat("  [SKIP] Phospho (w/o PN) META_DPS file not found:",
            basename(without_PN_file), "\n\n")
        integration_summary[[comp_name]] <- data.frame(
            comparison = comp_name,
            n_phos_with_PN = NA, n_phos_without_PN = NA, n_merged_rows = NA,
            n_unique_phosphosites_merged = NA,
            n_both_sig_005 = NA, n_concordant_sig = NA, n_discordant_sig = NA,
            stringsAsFactors = FALSE
        )
        next
    }

    # --- Step 2: Read and prepare data ---
    with_PN_data <- fread(with_PN_file)
    without_PN_data <- fread(without_PN_file)

    cat("  Phosphosites (with PN):", nrow(with_PN_data), "\n")
    cat("  Phosphosites (w/o PN):", nrow(without_PN_data), "\n")

    # Select core columns for merging, prefixed by normalization approach
    # With protein normalization core columns
    with_PN_core <- with_PN_data[, .(
        phosphosite,
        gene_symbol_phosphosite,
        with_PN_n_studies   = n_studies,
        with_PN_total_n     = total_n,
        with_PN_Z_meta      = Z_meta,
        with_PN_p_meta      = p_meta,
        with_PN_padj        = padj,
        with_PN_direction   = direction,
        with_PN_mean_logFC  = mean_logFC
    )]

    # Without protein normalization core columns
    without_PN_core <- without_PN_data[, .(
        phosphosite,
        gene_symbol_phosphosite,
        without_PN_n_studies   = n_studies,
        without_PN_total_n     = total_n,
        without_PN_Z_meta      = Z_meta,
        without_PN_p_meta      = p_meta,
        without_PN_padj        = padj,
        without_PN_direction   = direction,
        without_PN_mean_logFC  = mean_logFC
    )]

    # Remove duplicate phosphosite entries (keep first occurrence,
    # which is sorted by padj from the meta-analysis pipeline)
    with_PN_core <- with_PN_core[!duplicated(phosphosite)]
    without_PN_core <- without_PN_core[!duplicated(phosphosite)]

    # --- Step 3: Inner join by phosphosite ---
    # Both datasets share the same phosphosite identifier
    merged <- merge(with_PN_core, without_PN_core,
                    by = "phosphosite",
                    all.x = FALSE, all.y = FALSE,  # inner join
                    suffixes = c("_with_PN", "_without_PN"))

    # Resolve gene_symbol_phosphosite: prefer with_PN version
    merged$gene_symbol_phosphosite <- merged$gene_symbol_phosphosite_with_PN
    merged[, c("gene_symbol_phosphosite_with_PN",
               "gene_symbol_phosphosite_without_PN") := NULL]

    cat("  Merged rows (inner join):", nrow(merged), "\n")
    cat("  Unique phosphosites in merged data:",
        length(unique(merged$phosphosite)), "\n")

    if (nrow(merged) == 0) {
        cat("  [WARN] No overlapping phosphosites found. Skipping.\n\n")
        integration_summary[[comp_name]] <- data.frame(
            comparison = comp_name,
            n_phos_with_PN = nrow(with_PN_data),
            n_phos_without_PN = nrow(without_PN_data),
            n_merged_rows = 0,
            n_unique_phosphosites_merged = 0,
            n_both_sig_005 = 0, n_concordant_sig = 0, n_discordant_sig = 0,
            stringsAsFactors = FALSE
        )
        next
    }

    # --- Step 4: Compute derived columns ---
    # Z_meta_sum: concordant signal (with_PN Z_meta + without_PN Z_meta)
    merged$Z_meta_sum <- merged$with_PN_Z_meta + merged$without_PN_Z_meta

    # Z_meta_diff: discordance signal (without_PN Z_meta - with_PN Z_meta)
    merged$Z_meta_diff <- merged$without_PN_Z_meta - merged$with_PN_Z_meta

    # padj_max: conservative threshold = max(with_PN padj, without_PN padj)
    merged$padj_max <- pmax(merged$with_PN_padj, merged$without_PN_padj,
                            na.rm = TRUE)

    # --- Step 5: Arrange output column order ---
    col_order <- c(
        "gene_symbol_phosphosite", "phosphosite",
        "with_PN_Z_meta", "without_PN_Z_meta", "Z_meta_sum", "Z_meta_diff",
        "with_PN_padj", "without_PN_padj", "padj_max",
        "with_PN_direction", "without_PN_direction",
        "with_PN_mean_logFC", "without_PN_mean_logFC",
        "with_PN_n_studies", "without_PN_n_studies",
        "with_PN_total_n", "without_PN_total_n",
        "with_PN_p_meta", "without_PN_p_meta"
    )
    setcolorder(merged, col_order)

    # Sort by padj_max (most significant first)
    merged <- merged[order(padj_max,
                           pmax(with_PN_p_meta, without_PN_p_meta,
                                na.rm = TRUE))]

    # --- Step 6: Summary statistics ---
    n_both_sig <- sum(merged$padj_max < 0.05, na.rm = TRUE)

    # Concordant: same direction and padj_max < 0.05
    n_concordant <- sum(
        merged$padj_max < 0.05 &
        merged$with_PN_direction == merged$without_PN_direction,
        na.rm = TRUE
    )

    # Discordant: different direction and both individually significant
    n_discordant <- sum(
        merged$with_PN_padj < 0.05 &
        merged$without_PN_padj < 0.05 &
        merged$with_PN_direction != merged$without_PN_direction,
        na.rm = TRUE
    )

    cat("  Both significant (padj_max < 0.05):", n_both_sig, "\n")
    cat("    Concordant (same direction):", n_concordant, "\n")
    cat("    Discordant (opposite direction):", n_discordant, "\n")

    # --- Step 7: Save CSV ---
    out_csv <- file.path(output_dir,
                         paste0("META_integrated_", comp_name, ".csv"))
    fwrite(merged, out_csv)

    # --- Step 8: Save XLSX with formatting ---
    wb <- createWorkbook()
    addWorksheet(wb, comp_label)
    writeData(wb, 1, merged)

    header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#4472C4", fontColour = "white"
    )
    addStyle(wb, 1, header_style,
             rows = 1, cols = 1:ncol(merged), gridExpand = TRUE)
    setColWidths(wb, 1, cols = 1:ncol(merged), widths = "auto")

    # Highlight rows where padj_max < 0.05 (both approaches significant)
    sig_rows <- which(merged$padj_max < 0.05) + 1  # +1 for header
    if (length(sig_rows) > 0) {
        sig_style <- createStyle(fgFill = "#FFFFCC")
        addStyle(wb, 1, sig_style,
                 rows = sig_rows, cols = 1:ncol(merged), gridExpand = TRUE,
                 stack = TRUE)
    }

    out_xlsx <- file.path(output_dir,
                          paste0("META_integrated_", comp_name, ".xlsx"))
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    cat("  Saved:", basename(out_csv), "&", basename(out_xlsx), "\n\n")

    # Store summary
    integration_summary[[comp_name]] <- data.frame(
        comparison                  = comp_name,
        n_phos_with_PN              = nrow(with_PN_data),
        n_phos_without_PN           = nrow(without_PN_data),
        n_merged_rows               = nrow(merged),
        n_unique_phosphosites_merged = length(unique(merged$phosphosite)),
        n_both_sig_005              = n_both_sig,
        n_concordant_sig            = n_concordant,
        n_discordant_sig            = n_discordant,
        stringsAsFactors            = FALSE
    )
}

# ==============================================================================
# Section 3: Integration Summary
# ==============================================================================

cat("====================================================================\n")
cat("Generating Integration Summary\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(integration_summary)

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
saveWorkbook(wb_sum, file.path(output_dir, "integration_summary.xlsx"),
             overwrite = TRUE)

# Print summary
cat("\n=== Phospho DPS (with PN) vs (w/o PN) Meta Integration Summary ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Integration tables saved. Proceeding to scatter plots...\n")
cat("====================================================================\n")

# ==============================================================================
# Section 4: Phospho (w/o PN) Z_meta vs Phospho (with PN) Z_meta Scatter Plots
# ==============================================================================
#
# Axis mapping (following the template's protein/phospho role assignment):
#   x-axis: Phospho (with PN) Z_meta  (replaces the phospho w/o PN role)
#   y-axis: Phospho (w/o PN) Z_meta   (replaces the protein role)
#
# Highlight rules use the y-axis variable (without_PN_Z_meta),
# matching the template's behavior of highlighting by the phospho axis.

cat("\n====================================================================\n")
cat("Generating Phospho (w/o PN) vs Phospho (with PN) Z_meta Scatter Plots\n")
cat("====================================================================\n\n")

plot_dir <- file.path(output_dir, "scatter_plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Number of top/bottom phosphosites to highlight per category
n_top <- 10

# Load reference gene list (TP53 target genes from supplementary data)
# Source: Cao S, et al. Estimation of tumor cell total mRNA expression in
#   15 cancer types predicts disease progression. Nat Biotechnol.
#   2022;40(11):1624-1633. Supplementary Table 5.
ref_gene_file <- file.path(base_path, "41420_2023_1413_MOESM5_ESM.csv")
ref_genes <- fread(ref_gene_file, select = "gene")$gene
ref_genes <- unique(ref_genes[!is.na(ref_genes) & ref_genes != ""])
cat("  Reference gene list loaded:", length(ref_genes), "genes\n\n")

# purple color for reference gene list outline
purple_color <- "#B452CD"

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("Scatter plot for:", comp_name, "\n")

    # Read the integrated output file
    int_file <- file.path(output_dir,
                          paste0("META_integrated_", comp_name, ".csv"))
    if (!file.exists(int_file)) {
        cat("  [SKIP] Integrated file not found\n\n")
        next
    }

    df <- fread(int_file)

    # Keep only rows where both Z_meta values are finite
    df <- df[is.finite(without_PN_Z_meta) & is.finite(with_PN_Z_meta)]

    if (nrow(df) < 5) {
        cat("  [SKIP] Too few phosphosites with both Z_meta values:",
            nrow(df), "\n\n")
        next
    }

    # ------------------------------------------------------------------
    # Identify highlighted phosphosites (coloring rules)
    # ------------------------------------------------------------------

    # Rule 1: Top 10 largest with_PN_Z_meta -> orange
    #         (prerequisite: both with_PN and without_PN significantly UP,
    #          i.e., padj < 0.05 and direction == "up" for both)
    #         Top 10 smallest with_PN_Z_meta -> teal
    #         (prerequisite: both with_PN and without_PN significantly DOWN,
    #          i.e., padj < 0.05 and direction == "down" for both)
    df_sig_up <- df[with_PN_padj < 0.05 & with_PN_direction == "up" &
                    without_PN_padj < 0.05 & without_PN_direction == "up"]
    df_sig_up <- df_sig_up[order(-with_PN_Z_meta)]
    sites_phospho_top <- df_sig_up$phosphosite[
        1:min(n_top, nrow(df_sig_up))]

    df_sig_down <- df[with_PN_padj < 0.05 & with_PN_direction == "down" &
                      without_PN_padj < 0.05 & without_PN_direction == "down"]
    df_sig_down <- df_sig_down[order(with_PN_Z_meta)]
    sites_phospho_bot <- df_sig_down$phosphosite[
        1:min(n_top, nrow(df_sig_down))]

    # Rule 1 highlight phosphosites (union)
    rule1_sites <- unique(c(sites_phospho_top, sites_phospho_bot))

    # Rule 2: Among Rule 1 points, check if gene_symbol is in reference list
    #   -> keep Rule 1 color but add thick purple outline
    # Extract gene_symbol from gene_symbol_phosphosite (part before underscore)
    df$gene_symbol <- sub("_[^_]*$", "", df$gene_symbol_phosphosite)
    df$is_ref_gene <- df$gene_symbol %in% ref_genes

    # Assign fill color (Rule 4 first, then overwrite upward):
    # Rule 4 default: medium grey
    df$fill_color <- "medium_grey"

    # Rule 3: both with_PN_padj > 0.05 AND without_PN_padj > 0.05 -> light grey
    df$fill_color[df$without_PN_padj > 0.05 &
                  df$with_PN_padj > 0.05] <- "light_grey"

    # Rule 1 (highest priority, overwrites Rule 3/4):
    df$fill_color[df$phosphosite %in% sites_phospho_bot] <- "teal"
    df$fill_color[df$phosphosite %in% sites_phospho_top] <- "orange"

    # Rule 2: purple outline for Rule 1 points whose gene is in reference list
    df$has_purple_outline <- (df$phosphosite %in% rule1_sites) & df$is_ref_gene

    # Collect all labeled phosphosites (Rule 1 highlight groups)
    labeled_sites <- rule1_sites
    df$show_label <- df$phosphosite %in% labeled_sites

    # ------------------------------------------------------------------
    # Prepare plot data subsets
    # ------------------------------------------------------------------

    # Define fill color palette
    fill_values <- c(
        "orange"      = "#FF7F0E",
        "teal"        = "#009688",
        "light_grey"  = "#D9D9D9",
        "medium_grey" = "#A0A0A0"
    )

    # Split data into layers for proper z-order
    df_grey   <- df[fill_color %in% c("light_grey", "medium_grey")]
    df_rule1_no_outline <- df[fill_color %in% c("orange", "teal") &
                              has_purple_outline == FALSE]
    df_rule1_with_outline <- df[fill_color %in% c("orange", "teal") &
                                has_purple_outline == TRUE]

    # Axis range (symmetric for clarity)
    axis_max <- max(abs(c(df$with_PN_Z_meta, df$without_PN_Z_meta)),
                    na.rm = TRUE) * 1.05

    # Subset for labels (use gene_symbol_phosphosite for clarity)
    df_label <- df[show_label == TRUE]

    # Report highlight counts
    n_ref_in_rule1 <- sum(df$has_purple_outline)
    cat("  Rule 1 highlighted:", length(rule1_sites),
        "(", n_ref_in_rule1, "with reference gene overlap)\n")

    # Build the plot
    # x-axis: Phospho (with PN) Z_meta
    # y-axis: Phospho (w/o PN) Z_meta
    p <- ggplot(df, aes(x = with_PN_Z_meta, y = without_PN_Z_meta)) +
        # y = x reference line (drawn first, behind points)
        geom_abline(intercept = 0, slope = 1,
                    linetype = "dashed", color = "grey50", linewidth = 0.4) +
        # Layer 1: grey background points (Rule 3 & Rule 4)
        geom_point(
            data = df_grey,
            aes(color = fill_color),
            shape = 16, size = 0.6, alpha = 0.5
        ) +
        # Layer 2: Rule 1 highlighted points WITHOUT purple outline
        geom_point(
            data = df_rule1_no_outline,
            aes(color = fill_color),
            shape = 16, size = 1.8, alpha = 0.85
        ) +
        # Layer 3: Rule 1 highlighted points WITH purple outline (Rule 2)
        #          shape 21 = circle with separate fill and border color
        geom_point(
            data = df_rule1_with_outline,
            aes(fill = fill_color),
            shape = 21, color = purple_color,
            size = 2.2, alpha = 0.85, stroke = 0.8
        ) +
        # Phosphosite labels with connecting lines (ggrepel)
        # Uses gene_symbol_phosphosite (e.g., MCM6_S762) for unique labels
        geom_text_repel(
            data = df_label,
            aes(label = gene_symbol_phosphosite),
            size = 2.2,
            fontface = "italic",
            segment.size = 0.3,
            segment.color = "grey40",
            segment.alpha = 0.7,
            min.segment.length = 0,
            box.padding = 0.5,
            point.padding = 0.3,
            force = 8,
            force_pull = 0.5,
            max.overlaps = Inf,
            max.iter = 20000,
            seed = 42
        ) +
        scale_color_manual(values = fill_values, guide = "none") +
        scale_fill_manual(values = fill_values, guide = "none") +
        coord_cartesian(xlim = c(-axis_max, axis_max),
                        ylim = c(-axis_max, axis_max)) +
        labs(
            x = expression("Phospho (with PN) " * italic(Z)[meta]),
            y = expression("Phospho (w/o PN) " * italic(Z)[meta]),
            title = gsub("_", " ", comp_name)
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
            plot.title       = element_text(hjust = 0.5, face = "bold",
                                            size = 11),
            plot.margin      = margin(10, 10, 10, 10)
        )

    # Save TIFF (300 dpi, journal submission quality)
    tiff_file <- file.path(plot_dir, paste0("scatter_", comp_name, ".tiff"))
    tiff(tiff_file, width = 3.5, height = 3.5, units = "in",
         res = 300, compression = "lzw")
    print(p)
    dev.off()

    # Save PDF (vector format)
    pdf_file <- file.path(plot_dir, paste0("scatter_", comp_name, ".pdf"))
    pdf(pdf_file, width = 3.5, height = 3.5)
    print(p)
    dev.off()

    cat("  Saved:", basename(tiff_file), "&", basename(pdf_file), "\n")
    cat("  Highlighted phosphosites:", length(labeled_sites), "\n\n")
}

cat("====================================================================\n")
cat("All tasks complete!\n")
cat("Integration results saved to:", output_dir, "\n")
cat("Scatter plots saved to:", plot_dir, "\n")
cat("====================================================================\n")

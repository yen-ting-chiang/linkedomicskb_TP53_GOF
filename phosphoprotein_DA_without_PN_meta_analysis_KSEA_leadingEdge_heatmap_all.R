################################################################################
# KSEA Leading Edge Substrates Heatmap (pheatmap)
# (Meta-analysis KSEA leading edges x per-dataset DPS without PN logFC)
#
# Purpose:
#   For each selected KSEA kinase, extract its leading edge substrates from the
#   meta-KSEA results (without protein normalization), then look up the logFC of
#   each substrate's corresponding gene_symbol_phosphosite in the per-dataset
#   DPS results (without protein normalization). The resulting heatmap shows:
#     - X-axis: datasets (alphabetically sorted)
#     - Y-axis: leading edge substrates (alphabetically sorted)
#     - Color : logFC value from per-dataset DPS results
#
# Input:
#   1. META_KSEA_*.csv from
#        phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/
#        phosphoprotein_DPS_without_protein_normalization_meta_analysis/meta_KSEA/
#      These contain a "leadingEdge_substrates" column with semicolon-separated
#      substrate IDs (GENE_SITE format, e.g., "MCM6_S762;TP53BP1_S1678;...").
#
#   2. DPS_*.csv from
#        phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/
#        {dataset}/
#      These contain per-dataset phosphosite differential analysis results with
#      columns: gene_symbol_phosphosite, logFC, adj.P.Val, etc.
#
# Output:
#   phosphoprotein_DA_without_PN_meta_analysis_KSEA_leadingEdge_heatmap/
#     {color_scheme}/{comparison}/{kinase_name}_heatmap.tiff
#     {color_scheme}/{comparison}/{kinase_name}_heatmap.pdf
#     data/{comparison}/{kinase_name}_heatmap_data.csv
#     data/{comparison}/{kinase_name}_heatmap_data.xlsx
#
# Methodology references:
#   - Casado P, et al. Kinase-substrate enrichment analysis provides insights
#     into the heterogeneity of signaling pathway activation in leukemia cells.
#     Sci Signal. 2013;6(268):rs6. PMID: 23532336
#   - Wiredja DD, et al. The KSEA App: a web-based tool for kinase activity
#     inference from quantitative phosphoproteomics. Bioinformatics.
#     2017;33(21):3489-3491. PMID: 28655153
#   - Kolde R. pheatmap: Pretty Heatmaps. R package. CRAN.
################################################################################

# ==============================================================================
# Section 1: Setup and User Configuration
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(pheatmap)
    library(RColorBrewer)
    library(grDevices)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Meta-KSEA input directory (without protein normalization meta-analysis results)
meta_ksea_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted",
    "phosphoprotein_DPS_without_protein_normalization_meta_analysis",
    "meta_KSEA")

# Per-dataset DPS input directory (without protein normalization)
dps_input_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted")

# Output directory
output_dir <- file.path(base_path,
    "phosphoprotein_DA_without_PN_meta_analysis_KSEA_leadingEdge_heatmap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Available datasets
all_datasets <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
                   "LSCC", "LUAD", "OV", "PDAC", "UCEC")

# ==============================================================================
# Section 1a: User-configurable Comparisons
# ==============================================================================

# All 9 available comparisons
comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

# USER: Specify which comparisons to generate heatmaps for.
# Default: only TP53mt_vs_TP53wt. Uncomment others as needed.
selected_comparisons <- c(
    "TP53mt_vs_TP53wt"#,
    #"MUT_GOF_vs_MUT_LOF",
    #"Hotspot_vs_MUT_LOF",
    #"MUT_GOF_vs_TP53wt",
    #"MUT_LOF_vs_TP53wt",
    #"Hotspot_vs_TP53wt",
    #"DN_vs_TP53wt",
    #"NonDN_vs_TP53wt",
    #"DN_vs_NonDN"
)

# ==============================================================================
# Section 1b: User-configurable Kinases (KSEA Pathways)
# ==============================================================================

# USER: Specify which KSEA kinases to generate heatmaps for.
# Each kinase is plotted SEPARATELY.
# Default: CDK1 and CDK2.
selected_kinases <- c(
    "CDK1",
    "CDK2"#,
    #"CSNK2A1",
    #"AKT1",
    #"MAPK1",
    #"CHEK2",
    #"CDK6",
    #"ATR",
    #"CDK4"
)

# ==============================================================================
# Section 1c: User-configurable Color Schemes
# ==============================================================================

# Multiple journal-common diverging color palettes for heatmap.
# USER: Select which color scheme(s) to use by modifying selected_color_schemes.
# All schemes use a diverging palette centered at 0 (white/neutral).

color_schemes <- list(
    # Palette 1: Red-Blue (RdBu reversed) - Nature, Cell common
    "RdBu" = list(
        name  = "RdBu",
        label = "Red-Blue (Nature/Cell style)",
        low   = "#2166AC",   # Blue (negative logFC)
        mid   = "#F7F7F7",   # White (zero)
        high  = "#B2182B"    # Red (positive logFC)
    ),
    # Palette 2: Red-Yellow-Blue (RdYlBu reversed) - Science common
    "RdYlBu" = list(
        name  = "RdYlBu",
        label = "Red-Yellow-Blue (Science style)",
        low   = "#313695",   # Dark blue
        mid   = "#FFFFBF",   # Yellow
        high  = "#A50026"    # Dark red
    ),
    # Palette 3: Purple-Orange (PuOr) - colorblind-friendly
    "PuOr" = list(
        name  = "PuOr",
        label = "Purple-Orange (colorblind-friendly)",
        low   = "#542788",   # Purple (negative)
        mid   = "#F7F7F7",   # White
        high  = "#E08214"    # Orange (positive)
    ),
    # Palette 4: Teal-Red (custom) - high-impact journal style
    "TealRed" = list(
        name  = "TealRed",
        label = "Teal-Red (Lancet/NEJM style)",
        low   = "#008080",   # Teal (negative)
        mid   = "#F5F5F5",   # Near-white
        high  = "#CD2626"    # Red (positive)
    ),
    # Palette 5: Blue-White-Red (BWR) - classic heatmap
    "BWR" = list(
        name  = "BWR",
        label = "Blue-White-Red (classic)",
        low   = "#0571B0",   # Blue
        mid   = "#FFFFFF",   # White
        high  = "#CA0020"    # Red
    ),
    # Palette 6: Green-Black-Red (GBR) - microarray legacy
    "GBR" = list(
        name  = "GBR",
        label = "Green-Black-Red (microarray legacy)",
        low   = "#008837",   # Green (negative)
        mid   = "#000000",   # Black (zero)
        high  = "#C51B7D"    # Magenta-Red (positive)
    ),
    # Palette 7: Viridis-inspired diverging
    "ViridisDiv" = list(
        name  = "ViridisDiv",
        label = "Teal-Ivory-Magenta (Viridis-inspired)",
        low   = "#21908C",   # Teal
        mid   = "#FCFDBF",   # Ivory
        high  = "#9C179E"    # Magenta
    ),
    # Palette 8: Spectral reversed
    "Spectral" = list(
        name  = "Spectral",
        label = "Blue-Yellow-Red (Spectral)",
        low   = "#3288BD",   # Blue
        mid   = "#FFFFBF",   # Yellow
        high  = "#D53E4F"    # Red
    )
)

# USER: Specify which color schemes to generate.
# Select one or more from: "RdBu", "RdYlBu", "PuOr", "TealRed",
#                           "BWR", "GBR", "ViridisDiv", "Spectral"
# To generate all schemes, use: names(color_schemes)
selected_color_schemes <- names(color_schemes)

# ==============================================================================
# Section 1d: Pheatmap Visual Settings
# ==============================================================================

# USER: Customize pheatmap appearance parameters

# Font sizes
fontsize_row   <- 8    # Y-axis label font size (leading edge substrates)
fontsize_col   <- 9    # X-axis label font size (datasets)
fontsize_main  <- 11   # Title font size
fontsize_legend <- 8   # Legend font size

# Cell dimensions (in points)
cellwidth  <- 30
cellheight <- 12

# Show cell values on heatmap?
display_values <- FALSE

# Number format for cell values (if display_values = TRUE)
number_format <- "%.2f"

# Border color for cells
cell_border_color <- "grey80"

# NA cell color
na_color <- "grey90"

# Annotation for significance: mark cells with adj.P.Val < 0.05
# Set to TRUE to overlay "*" on significant cells
show_significance <- FALSE

# ==============================================================================
# Section 1e: pheatmap Output Dimensions
# ==============================================================================

# USER: Customize output figure dimensions.
# Set to "auto" for automatic calculation based on matrix size,
# or specify numeric values (inches).
fig_width  <- "auto"
fig_height <- "auto"

# TIFF resolution (dpi) for publication
tiff_dpi <- 300

cat("====================================================================\n")
cat("KSEA Leading Edge Heatmap (pheatmap)\n")
cat("Meta-KSEA leading edges x per-dataset DPS without PN logFC\n")
cat("====================================================================\n\n")

cat("Selected kinases:", paste(selected_kinases, collapse = ", "), "\n")
cat("Selected comparisons:", paste(selected_comparisons, collapse = ", "), "\n")
cat("Selected color schemes:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

# ==============================================================================
# Section 2: Helper Functions
# ==============================================================================

#' Parse leading edge substrates from KSEA result
#'
#' The leadingEdge_substrates column contains semicolon-separated substrate IDs.
#' Returns a character vector of substrate IDs (GENE_SITE format).
#'
#' @param le_string A single character string (semicolon-separated)
#' @return Character vector of substrate IDs (or empty if input is NA/"_"/empty)
parse_leading_edge <- function(le_string) {
    if (is.na(le_string) || le_string == "" || le_string == "_") {
        return(character(0))
    }
    trimws(unlist(strsplit(le_string, ";")))
}

#' Generate a diverging color palette for pheatmap
#'
#' Creates a smooth color gradient from low -> mid -> high with n_colors steps.
#'
#' @param low  Hex color for negative extreme
#' @param mid  Hex color for zero/center
#' @param high Hex color for positive extreme
#' @param n_colors Total number of colors in the palette (default 100)
#' @return Character vector of hex colors
make_diverging_palette <- function(low, mid, high, n_colors = 100) {
    half <- floor(n_colors / 2)
    c(
        colorRampPalette(c(low, mid))(half),
        colorRampPalette(c(mid, high))(n_colors - half)
    )
}

# ==============================================================================
# Section 3: Main Processing Loop
# ==============================================================================

# Filter and validate selected comparisons
comparisons <- comparisons[comparisons$name %in% selected_comparisons, , drop = FALSE]
if (nrow(comparisons) == 0) {
    stop("No valid comparisons selected. Please choose from: ",
         paste(c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
                 "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
                 "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"), collapse = ", "))
}

# Validate selected color schemes
selected_color_schemes <- intersect(selected_color_schemes, names(color_schemes))
if (length(selected_color_schemes) == 0) {
    stop("No valid color schemes selected. Please choose from: ",
         paste(names(color_schemes), collapse = ", "))
}

# Create output subdirectories
for (cs_name in selected_color_schemes) {
    dir.create(file.path(output_dir, cs_name), recursive = TRUE, showWarnings = FALSE)
}
data_out_dir <- file.path(output_dir, "data")
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

for (j_comp in seq_len(nrow(comparisons))) {
    comp_name  <- comparisons$name[j_comp]
    comp_label <- comparisons$label[j_comp]

    cat("==================================================================\n")
    cat("Comparison:", comp_name, "\n")
    cat("==================================================================\n")

    # ------------------------------------------------------------------
    # Step 1: Read META_KSEA file to get leading edge substrates
    # ------------------------------------------------------------------
    ksea_file <- file.path(meta_ksea_dir, paste0("META_KSEA_", comp_name, ".csv"))
    if (!file.exists(ksea_file)) {
        cat("  [SKIP] META_KSEA file not found:", basename(ksea_file), "\n\n")
        next
    }

    ksea <- fread(ksea_file)

    # Verify required columns
    if (!all(c("Kinase.Gene", "leadingEdge_substrates") %in% colnames(ksea))) {
        cat("  [SKIP] Missing required columns (Kinase.Gene, leadingEdge_substrates)\n\n")
        next
    }

    # ------------------------------------------------------------------
    # Step 2: Load per-dataset DPS data (without PN) for this comparison
    # ------------------------------------------------------------------
    cat("  Loading per-dataset DPS data...\n")

    dps_all <- list()

    for (ds in all_datasets) {
        dps_file <- file.path(dps_input_dir, ds, paste0("DPS_", comp_name, ".csv"))
        if (!file.exists(dps_file)) {
            cat("    [SKIP] DPS file not found for", ds, "\n")
            next
        }

        dps_df <- tryCatch(
            fread(dps_file),
            error = function(e) {
                cat("    [ERROR] Failed to read", ds, ":", e$message, "\n")
                NULL
            }
        )

        if (is.null(dps_df) || nrow(dps_df) == 0) next

        # Ensure required columns exist
        required_cols <- c("gene_symbol_phosphosite", "logFC")
        if (!all(required_cols %in% colnames(dps_df))) {
            cat("    [SKIP]", ds, ": missing required columns\n")
            next
        }

        # Keep only needed columns and add dataset identifier
        keep_cols <- intersect(c("gene_symbol_phosphosite", "logFC", "adj.P.Val"),
                               colnames(dps_df))
        dps_df <- dps_df[, ..keep_cols]
        dps_df$dataset <- ds

        dps_all[[ds]] <- dps_df
    }

    if (length(dps_all) == 0) {
        cat("  [SKIP] No DPS data available for any dataset\n\n")
        next
    }

    # Combine all datasets
    dps_combined <- rbindlist(dps_all, fill = TRUE)
    datasets_available <- sort(unique(dps_combined$dataset))
    cat("  Datasets loaded:", paste(datasets_available, collapse = ", "), "\n")

    # ------------------------------------------------------------------
    # Step 3: Generate heatmap for each selected kinase
    # ------------------------------------------------------------------
    for (kinase_name in selected_kinases) {

        cat("\n  *** Kinase:", kinase_name, "***\n")

        # Check if this kinase exists in KSEA results
        ksea_row <- ksea[ksea$Kinase.Gene == kinase_name, ]

        if (nrow(ksea_row) == 0) {
            cat("    [SKIP] Kinase not found in KSEA results\n")
            next
        }

        # Parse leading edge substrates
        le_str <- as.character(ksea_row$leadingEdge_substrates[1])
        le_subs <- parse_leading_edge(le_str)
        cat("    Leading edge substrates:", length(le_subs), "\n")

        if (length(le_subs) == 0) {
            cat("    [SKIP] No leading edge substrates\n")
            next
        }

        # Sort leading edge substrates alphabetically (for y-axis)
        le_subs <- sort(le_subs)

        # Filter DPS data to leading edge substrates
        plot_data <- dps_combined[gene_symbol_phosphosite %in% le_subs]

        if (nrow(plot_data) == 0) {
            cat("    [SKIP] No matching phosphosites found in DPS data\n")
            next
        }

        # Handle duplicates (keep first occurrence per phosphosite per dataset)
        plot_data <- plot_data[!duplicated(
            plot_data[, .(gene_symbol_phosphosite, dataset)])]

        # Report coverage
        sites_found <- sort(unique(plot_data$gene_symbol_phosphosite))
        datasets_found <- sort(unique(plot_data$dataset))
        cat("    Substrates found in DPS:", length(sites_found), "/",
            length(le_subs), "\n")
        cat("    Datasets with data:", length(datasets_found), "\n")

        # Create wide-format logFC matrix (substrates x datasets)
        logfc_wide <- dcast(plot_data, gene_symbol_phosphosite ~ dataset,
                            value.var = "logFC")

        # Build numeric matrix
        logfc_mat <- as.data.frame(logfc_wide)
        rownames(logfc_mat) <- logfc_mat$gene_symbol_phosphosite
        logfc_mat$gene_symbol_phosphosite <- NULL
        logfc_mat <- as.matrix(logfc_mat)

        # Sort rows alphabetically (y-axis: leading edge substrates)
        logfc_mat <- logfc_mat[order(rownames(logfc_mat)), , drop = FALSE]

        # Sort columns alphabetically (x-axis: datasets)
        logfc_mat <- logfc_mat[, order(colnames(logfc_mat)), drop = FALSE]

        # Also create adj.P.Val matrix if significance marking is enabled
        sig_mat <- NULL
        if (show_significance && "adj.P.Val" %in% colnames(plot_data)) {
            pval_wide <- dcast(plot_data, gene_symbol_phosphosite ~ dataset,
                               value.var = "adj.P.Val")
            pval_df <- as.data.frame(pval_wide)
            rownames(pval_df) <- pval_df$gene_symbol_phosphosite
            pval_df$gene_symbol_phosphosite <- NULL
            pval_mat <- as.matrix(pval_df)
            pval_mat <- pval_mat[rownames(logfc_mat), colnames(logfc_mat),
                                 drop = FALSE]
            sig_mat <- matrix("", nrow = nrow(logfc_mat), ncol = ncol(logfc_mat),
                              dimnames = dimnames(logfc_mat))
            sig_mat[!is.na(pval_mat) & pval_mat < 0.05] <- "*"
        }

        # Determine symmetric color scale limits
        logfc_max <- max(abs(logfc_mat), na.rm = TRUE)
        logfc_lim <- ceiling(logfc_max * 10) / 10  # round up to 1 decimal

        # ------------------------------------------------------------------
        # Step 3a: Save data CSV and XLSX
        # ------------------------------------------------------------------
        comp_data_dir <- file.path(data_out_dir, comp_name)
        dir.create(comp_data_dir, recursive = TRUE, showWarnings = FALSE)

        # Prepare export data: add row names as first column
        export_df <- data.frame(
            leadingEdge_substrate = rownames(logfc_mat),
            logfc_mat,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        # Save CSV
        data_csv <- file.path(comp_data_dir,
                              paste0(kinase_name, "_heatmap_data.csv"))
        fwrite(export_df, data_csv)

        # Save XLSX with formatting
        data_xlsx <- file.path(comp_data_dir,
                               paste0(kinase_name, "_heatmap_data.xlsx"))
        wb_data <- createWorkbook()
        addWorksheet(wb_data, "Heatmap_Data")
        writeData(wb_data, 1, export_df)

        header_style <- createStyle(
            fontName = "Arial", fontSize = 8,
            textDecoration = "bold", halign = "center", valign = "center",
            border = "bottom", fgFill = "#4472C4", fontColour = "white"
        )
        body_style <- createStyle(
            fontName = "Arial", fontSize = 8,
            halign = "left", valign = "center"
        )
        addStyle(wb_data, 1, header_style,
                 rows = 1, cols = 1:ncol(export_df), gridExpand = TRUE)
        if (nrow(export_df) > 0) {
            addStyle(wb_data, 1, body_style,
                     rows = 2:(nrow(export_df) + 1), cols = 1:ncol(export_df),
                     gridExpand = TRUE, stack = FALSE)
        }
        setColWidths(wb_data, 1, cols = 1:ncol(export_df), widths = "auto")
        saveWorkbook(wb_data, data_xlsx, overwrite = TRUE)

        cat("    Data saved:", basename(data_csv), "&", basename(data_xlsx), "\n")

        # ------------------------------------------------------------------
        # Step 3b: Generate pheatmaps for each color scheme
        # ------------------------------------------------------------------
        for (cs_name in selected_color_schemes) {
            cs <- color_schemes[[cs_name]]

            # Generate diverging color palette
            heatmap_colors <- make_diverging_palette(cs$low, cs$mid, cs$high,
                                                     n_colors = 100)

            # Set symmetric breaks
            breaks <- seq(-logfc_lim, logfc_lim, length.out = 101)

            # Determine cell display values
            display_mat <- if (display_values) logfc_mat else FALSE

            # Calculate figure dimensions
            n_rows <- nrow(logfc_mat)
            n_cols <- ncol(logfc_mat)

            if (fig_width == "auto") {
                w <- max(4, 2.5 + n_cols * (cellwidth / 72))
            } else {
                w <- as.numeric(fig_width)
            }

            if (fig_height == "auto") {
                h <- max(4, 2.0 + n_rows * (cellheight / 72))
            } else {
                h <- as.numeric(fig_height)
            }

            # Construct annotation title
            main_title <- paste0(kinase_name, " KSEA Leading Edge Substrates\n(",
                                 gsub("_", " ", comp_name), ", ", cs$label, ")")

            # Create comparison subdirectory within color scheme directory
            comp_cs_dir <- file.path(output_dir, cs_name, comp_name)
            dir.create(comp_cs_dir, recursive = TRUE, showWarnings = FALSE)

            # --- Save TIFF ---
            tiff_file <- file.path(comp_cs_dir,
                                   paste0(kinase_name, "_heatmap.tiff"))
            tiff(tiff_file, width = w, height = h,
                 units = "in", res = tiff_dpi, compression = "lzw")

            pheatmap(
                mat               = logfc_mat,
                color             = heatmap_colors,
                breaks            = breaks,
                cluster_rows      = FALSE,     # Alphabetical order
                cluster_cols      = FALSE,     # Alphabetical order
                display_numbers   = display_mat,
                number_format     = number_format,
                number_color      = "black",
                fontsize_row      = fontsize_row,
                fontsize_col      = fontsize_col,
                fontsize          = fontsize_legend,
                cellwidth         = cellwidth,
                cellheight        = cellheight,
                border_color      = cell_border_color,
                na_col            = na_color,
                main              = main_title,
                angle_col         = 45,
                legend             = TRUE,
                legend_breaks      = pretty(c(-logfc_lim, logfc_lim), n = 5),
                legend_labels      = pretty(c(-logfc_lim, logfc_lim), n = 5)
            )

            dev.off()

            # --- Save PDF ---
            pdf_file <- file.path(comp_cs_dir,
                                  paste0(kinase_name, "_heatmap.pdf"))
            pdf(pdf_file, width = w, height = h)

            pheatmap(
                mat               = logfc_mat,
                color             = heatmap_colors,
                breaks            = breaks,
                cluster_rows      = FALSE,
                cluster_cols      = FALSE,
                display_numbers   = display_mat,
                number_format     = number_format,
                number_color      = "black",
                fontsize_row      = fontsize_row,
                fontsize_col      = fontsize_col,
                fontsize          = fontsize_legend,
                cellwidth         = cellwidth,
                cellheight        = cellheight,
                border_color      = cell_border_color,
                na_col            = na_color,
                main              = main_title,
                angle_col         = 45,
                legend             = TRUE,
                legend_breaks      = pretty(c(-logfc_lim, logfc_lim), n = 5),
                legend_labels      = pretty(c(-logfc_lim, logfc_lim), n = 5)
            )

            dev.off()
        }

        cat("    Heatmaps saved for all selected color schemes\n")
    }

    cat("\n")
}

cat("====================================================================\n")
cat("KSEA Leading Edge Heatmap pipeline completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

################################################################################
# Top Z_meta Phosphosites Heatmap (pheatmap)
# (Meta-analysis top Z_meta phosphosites x per-dataset logFC)
#
# Purpose:
#   For each selected comparison, extract the top N phosphosites ranked by
#   Z_meta from the meta-analysis results (without protein normalization),
#   then display the per-dataset logFC values as a heatmap.
#     - X-axis: datasets (alphabetically sorted)
#     - Y-axis: top N gene_symbol_phosphosite, ordered by Z_meta (descending)
#     - Color : logFC value from each dataset (logFC_* columns in META_DPS)
#
# Input:
#   META_DPS_*.csv from
#     phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/
#     phosphoprotein_DPS_without_protein_normalization_meta_analysis/
#   These contain per-phosphosite meta-analysis results with Z_meta scores
#   and per-dataset logFC values (logFC_BRCA, logFC_COAD, etc.).
#
# Output:
#   phosphoprotein_DA_without_PN_meta_analysis_total_heatmap/
#     {color_scheme}/{comparison}/top{N}_Z_meta_heatmap.tiff
#     {color_scheme}/{comparison}/top{N}_Z_meta_heatmap.pdf
#     data/{comparison}/top{N}_Z_meta_heatmap_data.csv
#     data/{comparison}/top{N}_Z_meta_heatmap_data.xlsx
#
# Methodology references:
#   - Kolde R. pheatmap: Pretty Heatmaps. R package. CRAN.
#   - Stouffer SA, et al. The American Soldier. Princeton University Press. 1949.
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
    library(msigdbr)
    library(KSEAapp)
    library(grid)
    library(gtable)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Meta-analysis input directory (contains META_DPS_*.csv files)
meta_dps_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted",
    "phosphoprotein_DPS_without_protein_normalization_meta_analysis")

# Output directory
output_dir <- file.path(base_path,
    "phosphoprotein_DA_without_PN_meta_analysis_total_heatmap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
# Section 1b: Top N Phosphosites by Z_meta (Y-axis limit)
# ==============================================================================

# USER: Specify how many top-ranked phosphosites to display.
# Phosphosites are ranked by Z_meta (descending) from the cross-cancer
# meta-analysis (META_DPS file). Only the top N phosphosites will be
# shown on the heatmap y-axis.
# Set to Inf to show all phosphosites without filtering.
top_n_sites <- 30

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
fontsize_row   <- 8    # Y-axis label font size (gene_symbol_phosphosite)
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

# ==============================================================================
# Section 1f: Row Annotation Settings
# ==============================================================================

# USER: Customize row annotation colors.
# 7 annotation columns are added to the left side of the heatmap:
#
# Columns 1-3: MSigDB Hallmark gene set membership (parent gene level)
#   1. E2F_TARGETS       : HALLMARK_E2F_TARGETS
#   2. G2M_CHECKPOINT    : HALLMARK_G2M_CHECKPOINT
#   3. MYC_TARGETS_V1    : HALLMARK_MYC_TARGETS_V1
#
# Columns 4-5: MSigDB C2:CP:KEGG_LEGACY pathway membership (parent gene level)
#   4. CELL_CYCLE        : KEGG_CELL_CYCLE
#   5. SPLICEOSOME       : KEGG_SPLICEOSOME
#
# Columns 6-7: KSEA kinase substrate membership (phosphosite level)
#   6. CDK1_substrate   : Whether the phosphosite is a known CDK1 substrate
#   7. CDK2_substrate   : Whether the phosphosite is a known CDK2 substrate
#
# Available color presets (journal-common palettes):
#   NPG (Nature):     "#E64B35" (red),    "#4DBBD5" (cyan)
#   Lancet:           "#00468B" (navy),    "#ED0000" (red)
#   NEJM:             "#BC3C29" (brick),   "#0072B5" (blue)
#   JAMA:             "#374E55" (charcoal),"#DF8F44" (amber)
#   JCO:              "#0073C2" (blue),    "#EFC000" (gold)
#   Cell:             "#7570B3" (purple),  "#D95F02" (orange)
#
# Set the fill color for "Yes" (member). "No" is always white (no fill).

# Hallmark annotations
annot_color_E2F        <- "#E64B35"     # NPG red
annot_color_G2M        <- "#4DBBD5"     # NPG cyan
annot_color_MYC_V1     <- "#00468B"     # Lancet navy

# KEGG_LEGACY annotations
annot_color_CELL_CYCLE    <- "#D95F02"  # Cell orange
annot_color_SPLICEOSOME   <- "#EFC000"  # JCO gold

# KSEA kinase substrate annotations
annot_color_CDK1       <- "#BC3C29"     # NEJM brick
annot_color_CDK2       <- "#ED0000"     # Lancet red

# ==============================================================================
# Section 2: Helper Functions
# ==============================================================================

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

cat("====================================================================\n")
cat("Top Z_meta Phosphosites Heatmap (pheatmap)\n")
cat("Meta-analysis top Z_meta phosphosites x per-dataset logFC\n")
cat("====================================================================\n\n")

cat("Selected comparisons:", paste(selected_comparisons, collapse = ", "), "\n")
cat("Top N phosphosites:", top_n_sites, "\n")
cat("Selected color schemes:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

# ==============================================================================
# Section 2a: Load MSigDB Gene Sets and KSEA Kinase-Substrate Data
# ==============================================================================

cat("Loading MSigDB gene sets for row annotations...\n")

# --- Hallmark gene sets ---
hallmark_df <- suppressMessages(
    msigdbr(species = "Homo sapiens", collection = "H")
)

e2f_genes <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_E2F_TARGETS"]
)
g2m_genes <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_G2M_CHECKPOINT"]
)
myc_v1_genes <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_MYC_TARGETS_V1"]
)

cat("  HALLMARK_E2F_TARGETS:", length(e2f_genes), "genes\n")
cat("  HALLMARK_G2M_CHECKPOINT:", length(g2m_genes), "genes\n")
cat("  HALLMARK_MYC_TARGETS_V1:", length(myc_v1_genes), "genes\n")

# --- KEGG_LEGACY gene sets (C2:CP:KEGG_LEGACY) ---
cat("\nLoading MSigDB C2:CP:KEGG_LEGACY gene sets...\n")

kegg_df <- suppressMessages(
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
)

cell_cycle_genes <- unique(
    kegg_df$gene_symbol[kegg_df$gs_name == "KEGG_CELL_CYCLE"]
)
spliceosome_genes <- unique(
    kegg_df$gene_symbol[kegg_df$gs_name == "KEGG_SPLICEOSOME"]
)

cat("  KEGG_CELL_CYCLE:", length(cell_cycle_genes), "genes\n")
cat("  KEGG_SPLICEOSOME:", length(spliceosome_genes), "genes\n")

# --- KSEA kinase-substrate data (CDK1 and CDK2) ---
cat("\nLoading KSEA kinase-substrate data from KSEAapp KSData...\n")

data("KSData", package = "KSEAapp", envir = environment())

# Build substrate IDs in GENE_SITE format (e.g., "TP53_S315")
# to match the gene_symbol_phosphosite column
cdk1_substrates <- unique(paste0(
    KSData$SUB_GENE[KSData$GENE == "CDK1"], "_",
    KSData$SUB_MOD_RSD[KSData$GENE == "CDK1"]
))
cdk2_substrates <- unique(paste0(
    KSData$SUB_GENE[KSData$GENE == "CDK2"], "_",
    KSData$SUB_MOD_RSD[KSData$GENE == "CDK2"]
))

cat("  CDK1 substrates:", length(cdk1_substrates), "\n")
cat("  CDK2 substrates:", length(cdk2_substrates), "\n\n")

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
    # Step 1: Read META_DPS file
    # ------------------------------------------------------------------
    meta_dps_file <- file.path(meta_dps_dir,
                               paste0("META_DPS_", comp_name, ".csv"))
    if (!file.exists(meta_dps_file)) {
        cat("  [SKIP] META_DPS file not found:", basename(meta_dps_file), "\n\n")
        next
    }

    meta_dps <- tryCatch(
        fread(meta_dps_file),
        error = function(e) {
            cat("  [ERROR] Failed to read META_DPS file:", e$message, "\n\n")
            NULL
        }
    )

    if (is.null(meta_dps) || nrow(meta_dps) == 0) {
        cat("  [SKIP] META_DPS file is empty\n\n")
        next
    }

    # Verify required columns
    if (!"gene_symbol_phosphosite" %in% colnames(meta_dps) ||
        !"Z_meta" %in% colnames(meta_dps)) {
        cat("  [SKIP] Missing required columns (gene_symbol_phosphosite, Z_meta)\n\n")
        next
    }

    cat("  META_DPS loaded:", nrow(meta_dps), "phosphosites\n")

    # ------------------------------------------------------------------
    # Step 2: Identify logFC columns (per-dataset logFC values)
    # ------------------------------------------------------------------
    logfc_cols <- grep("^logFC_", colnames(meta_dps), value = TRUE)

    if (length(logfc_cols) == 0) {
        cat("  [SKIP] No logFC_* columns found in META_DPS file\n\n")
        next
    }

    # Extract dataset names from logFC column names (e.g., "logFC_BRCA" -> "BRCA")
    dataset_names <- sub("^logFC_", "", logfc_cols)
    cat("  Datasets found:", paste(sort(dataset_names), collapse = ", "), "\n")

    # ------------------------------------------------------------------
    # Step 3: Select top N phosphosites by Z_meta (descending)
    # ------------------------------------------------------------------
    # Remove rows with NA Z_meta
    meta_dps <- meta_dps[!is.na(Z_meta)]

    # Sort by Z_meta descending
    meta_dps <- meta_dps[order(-Z_meta)]

    # Apply top N filter
    n_total <- nrow(meta_dps)
    if (is.finite(top_n_sites) && n_total > top_n_sites) {
        meta_dps_top <- meta_dps[1:top_n_sites]
    } else {
        meta_dps_top <- meta_dps
    }

    cat("  Top", nrow(meta_dps_top), "phosphosites selected (of", n_total, "total)\n")

    # ------------------------------------------------------------------
    # Step 4: Build logFC matrix (phosphosites x datasets)
    # ------------------------------------------------------------------
    logfc_mat <- as.data.frame(meta_dps_top[, ..logfc_cols])
    rownames(logfc_mat) <- meta_dps_top$gene_symbol_phosphosite

    # Rename columns from "logFC_BRCA" to "BRCA"
    colnames(logfc_mat) <- sub("^logFC_", "", colnames(logfc_mat))

    # Convert to numeric matrix
    logfc_mat <- as.matrix(logfc_mat)

    # Order rows by Z_meta descending (already sorted, highest at top)
    # Order columns alphabetically (x-axis: datasets)
    logfc_mat <- logfc_mat[, order(colnames(logfc_mat)), drop = FALSE]

    cat("  Matrix dimensions:", nrow(logfc_mat), "rows x",
        ncol(logfc_mat), "cols\n")

    # Determine symmetric color scale limits
    logfc_max <- max(abs(logfc_mat), na.rm = TRUE)
    logfc_lim <- ceiling(logfc_max * 10) / 10  # round up to 1 decimal

    cat("  logFC range: [", round(-logfc_lim, 2), ",",
        round(logfc_lim, 2), "]\n")

    # ------------------------------------------------------------------
    # Step 4b: Build row annotations (gene set & kinase substrate membership)
    # ------------------------------------------------------------------
    # Extract parent gene from phosphosite ID (e.g., "MCM2_S41" -> "MCM2")
    parent_genes <- sub("_[STY][0-9]+$", "", rownames(logfc_mat))

    annot_row <- data.frame(
        # KSEA kinase substrates (columns 6-7, phosphosite level)
        CDK2_substrate = ifelse(rownames(logfc_mat) %in% cdk2_substrates, "Yes", "No"),
        CDK1_substrate = ifelse(rownames(logfc_mat) %in% cdk1_substrates, "Yes", "No"),
        # KEGG_LEGACY pathways (columns 4-5, parent gene level)
        SPLICEOSOME         = ifelse(parent_genes %in% spliceosome_genes, "Yes", "No"),
        CELL_CYCLE          = ifelse(parent_genes %in% cell_cycle_genes, "Yes", "No"),
        # Hallmark gene sets (columns 1-3, parent gene level)
        MYC_TARGETS_V1 = ifelse(parent_genes %in% myc_v1_genes, "Yes", "No"),
        G2M_CHECKPOINT = ifelse(parent_genes %in% g2m_genes, "Yes", "No"),
        E2F_TARGETS    = ifelse(parent_genes %in% e2f_genes, "Yes", "No"),
        row.names = rownames(logfc_mat),
        stringsAsFactors = FALSE
    )

    # Convert all columns to factors for consistent pheatmap color mapping
    for (col_name in colnames(annot_row)) {
        annot_row[[col_name]] <- factor(annot_row[[col_name]], levels = c("Yes", "No"))
    }

    # Annotation color mapping
    annot_colors <- list(
        CDK2_substrate       = c("Yes" = annot_color_CDK2,       "No" = "white"),
        CDK1_substrate       = c("Yes" = annot_color_CDK1,       "No" = "white"),
        SPLICEOSOME          = c("Yes" = annot_color_SPLICEOSOME,  "No" = "white"),
        CELL_CYCLE           = c("Yes" = annot_color_CELL_CYCLE,    "No" = "white"),
        MYC_TARGETS_V1       = c("Yes" = annot_color_MYC_V1,     "No" = "white"),
        G2M_CHECKPOINT       = c("Yes" = annot_color_G2M,        "No" = "white"),
        E2F_TARGETS          = c("Yes" = annot_color_E2F,        "No" = "white")
    )

    # Report annotation counts
    cat("  Row annotation membership counts:\n")
    for (col_name in colnames(annot_row)) {
        n_yes <- sum(annot_row[[col_name]] == "Yes")
        cat("    ", col_name, ":", n_yes, "/", nrow(logfc_mat), "\n")
    }

    # ------------------------------------------------------------------
    # Step 5: Save data CSV and XLSX
    # ------------------------------------------------------------------
    comp_data_dir <- file.path(data_out_dir, comp_name)
    dir.create(comp_data_dir, recursive = TRUE, showWarnings = FALSE)

    # Prepare export data: add row names as first column plus Z_meta
    export_df <- data.frame(
        gene_symbol_phosphosite = rownames(logfc_mat),
        Z_meta = meta_dps_top$Z_meta,
        logfc_mat,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    # Save CSV
    data_csv <- file.path(comp_data_dir,
                          paste0("top", nrow(logfc_mat), "_Z_meta_heatmap_data.csv"))
    fwrite(export_df, data_csv)

    # Save XLSX with formatting
    data_xlsx <- file.path(comp_data_dir,
                           paste0("top", nrow(logfc_mat), "_Z_meta_heatmap_data.xlsx"))
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

    cat("  Data saved:", basename(data_csv), "&", basename(data_xlsx), "\n")

    # ------------------------------------------------------------------
    # Step 6: Generate pheatmaps for each color scheme
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
            # Add extra width for row annotations to prevent truncation on the left
            # Approx 0.5 inches per annotation column + 1.0 inch for annotation legends
            annot_width <- ncol(annot_row) * 0.5 + 1.0
            w <- max(4, 2.5 + n_cols * (cellwidth / 72) + annot_width)
        } else {
            w <- as.numeric(fig_width)
        }

        if (fig_height == "auto") {
            h <- max(4, 2.0 + n_rows * (cellheight / 72))
        } else {
            h <- as.numeric(fig_height)
        }

        # Construct title
        main_title <- paste0("Top ", nrow(logfc_mat),
                             " Phosphosites by Z_meta\n(",
                             gsub("_", " ", comp_name), ", ", cs$label, ")")

        # Create comparison subdirectory within color scheme directory
        comp_cs_dir <- file.path(output_dir, cs_name, comp_name)
        dir.create(comp_cs_dir, recursive = TRUE, showWarnings = FALSE)

        # ------------------------------------------------------------------
        # Step 6a: Generate the heatmap object silently
        # ------------------------------------------------------------------
        ph <- pheatmap(
            mat               = logfc_mat,
            color             = heatmap_colors,
            breaks            = breaks,
            cluster_rows      = FALSE,     # Z_meta rank order
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
            legend_labels      = pretty(c(-logfc_lim, logfc_lim), n = 5),
            annotation_row     = annot_row,
            annotation_colors  = annot_colors,
            silent             = TRUE      # Prevent drawing immediately
        )

        # ------------------------------------------------------------------
        # Step 6b: Manipulate gtable to move row names to the left
        # ------------------------------------------------------------------
        g <- ph$gtable
        row_names_idx <- which(g$layout$name == "row_names")

        if (length(row_names_idx) > 0) {
            # Right-align text and shift to the right edge of its cell
            # pheatmap might store text as textGrob or gTree depending on version
            if (inherits(g$grobs[[row_names_idx]], "text")) {
                g$grobs[[row_names_idx]]$hjust <- 1
                g$grobs[[row_names_idx]]$x <- rep(unit(1, "npc") - unit(3, "points"),
                                                  length(g$grobs[[row_names_idx]]$x))
            } else if (inherits(g$grobs[[row_names_idx]], "gTree")) {
                text_idx <- grep("text", names(g$grobs[[row_names_idx]]$children))
                if (length(text_idx) > 0) {
                    g$grobs[[row_names_idx]]$children[[text_idx]]$hjust <- 1
                    g$grobs[[row_names_idx]]$children[[text_idx]]$x <-
                        rep(unit(1, "npc") - unit(3, "points"),
                            length(g$grobs[[row_names_idx]]$children[[text_idx]]$x))
                }
            }

            # Get original column of row_names
            orig_col <- g$layout[row_names_idx, "l"]

            # Move to the first column (left side)
            g$layout[row_names_idx, c("l", "r")] <- 1

            # Swap widths between original column and column 1
            orig_width_1 <- g$widths[1]
            g$widths[1] <- g$widths[orig_col]
            g$widths[orig_col] <- orig_width_1
        }

        # ------------------------------------------------------------------
        # Step 6c: Draw to Output Devices
        # ------------------------------------------------------------------

        # --- Save TIFF ---
        tiff_file <- file.path(comp_cs_dir,
                               paste0("top", nrow(logfc_mat), "_Z_meta_heatmap.tiff"))
        tiff(tiff_file, width = w, height = h,
             units = "in", res = tiff_dpi, compression = "lzw")
        grid.newpage()
        grid.draw(g)
        dev.off()

        # --- Save PDF ---
        pdf_file <- file.path(comp_cs_dir,
                              paste0("top", nrow(logfc_mat), "_Z_meta_heatmap.pdf"))
        pdf(pdf_file, width = w, height = h)
        grid.newpage()
        grid.draw(g)
        dev.off()
    }

    cat("  Heatmaps saved for all selected color schemes\n")

    cat("\n")
}

cat("====================================================================\n")
cat("Top Z_meta Phosphosites Heatmap pipeline completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

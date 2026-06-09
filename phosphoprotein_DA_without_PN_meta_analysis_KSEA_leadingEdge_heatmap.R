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
#     - Y-axis: top N leading edge substrates, ordered by their KSEA rank
#               (ranked by absolute Z_meta from meta-analysis, highest first)
#     - Color : logFC value from per-dataset DPS results
#
# Input:
#   1. META_KSEA_*.csv from
#        phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/
#        phosphoprotein_DPS_without_protein_normalization_meta_analysis/meta_KSEA/
#      These contain a "leadingEdge_substrates" column with semicolon-separated
#      substrate IDs (GENE_SITE format, e.g., "MCM6_S762;TP53BP1_S1678;...").
#
#   2. META_DPS_*.csv from the same meta-analysis directory
#      These contain per-phosphosite meta-analysis results with Z_meta scores,
#      used to rank the leading edge substrates for y-axis ordering.
#
#   3. DPS_*.csv from
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
    library(msigdbr)
    library(grid)
    library(gtable)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Meta-analysis base directory (without protein normalization)
meta_base_dir <- file.path(base_path,
    "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted",
    "phosphoprotein_DPS_without_protein_normalization_meta_analysis")

# Meta-KSEA input directory (contains META_KSEA_*.csv with leading edge substrates)
meta_ksea_dir <- file.path(meta_base_dir, "meta_KSEA")

# Meta-DPS directory (contains META_DPS_*.csv with Z_meta for substrate ranking)
meta_dps_dir <- meta_base_dir

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
# Section 1b2: Top N Leading Edge Substrates (Y-axis limit)
# ==============================================================================

# USER: Specify how many top-ranked leading edge substrates to display.
# Substrates are ranked by their absolute Z_meta from the cross-cancer
# meta-analysis (META_DPS file). Only the top N substrates (by KSEA rank)
# will be shown on the heatmap y-axis.
# Set to Inf to show all leading edge substrates without filtering.
top_n_substrates <- 20

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

# ==============================================================================
# Section 1f: Row Annotation Settings (Hallmark Gene Set Membership)
# ==============================================================================

# USER: Customize row annotation colors for Hallmark gene set membership.
# Two annotation columns are added to the heatmap:
#   1. E2F_TARGETS: whether the substrate's parent gene is in HALLMARK_E2F_TARGETS
#   2. G2M_CHECKPOINT: whether the substrate's parent gene is in HALLMARK_G2M_CHECKPOINT
#
# Available color presets (journal-common palettes):
#   NPG (Nature):     "#E64B35" (red),    "#4DBBD5" (cyan)
#   Lancet:           "#00468B" (navy),    "#ED0000" (red)
#   NEJM:             "#BC3C29" (brick),   "#0072B5" (blue)
#   JAMA:             "#374E55" (charcoal),"#DF8F44" (amber)
#   JCO:              "#0073C2" (blue),    "#EFC000" (gold)
#   Cell:             "#7570B3" (purple),  "#D95F02" (orange)
#
# Set the fill color for "Yes" (gene is a member).
# "No" is always displayed as white (no fill).

annot_color_E2F <- "#E64B35"     # NPG red for E2F_TARGETS
annot_color_G2M <- "#4DBBD5"     # NPG cyan for G2M_CHECKPOINT

cat("====================================================================\n")
cat("KSEA Leading Edge Heatmap (pheatmap)\n")
cat("Meta-KSEA leading edges x per-dataset DPS without PN logFC\n")
cat("====================================================================\n\n")

cat("Selected kinases:", paste(selected_kinases, collapse = ", "), "\n")
cat("Selected comparisons:", paste(selected_comparisons, collapse = ", "), "\n")
cat("Selected color schemes:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n\n")

# ==============================================================================
# Section 2a: Load MSigDB Hallmark Gene Sets for Row Annotations
# ==============================================================================

cat("Loading MSigDB Hallmark gene sets for row annotations...\n")

hallmark_df <- suppressMessages(
    msigdbr(species = "Homo sapiens", collection = "H")
)

# Extract gene lists for E2F_TARGETS and G2M_CHECKPOINT
e2f_genes <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_E2F_TARGETS"]
)
g2m_genes <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_G2M_CHECKPOINT"]
)

cat("  HALLMARK_E2F_TARGETS:", length(e2f_genes), "genes\n")
cat("  HALLMARK_G2M_CHECKPOINT:", length(g2m_genes), "genes\n\n")

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
    # Step 2b: Load META_DPS data for substrate ranking (Z_meta)
    # ------------------------------------------------------------------
    meta_dps_file <- file.path(meta_dps_dir,
                               paste0("META_DPS_", comp_name, ".csv"))
    meta_dps <- NULL
    if (file.exists(meta_dps_file)) {
        meta_dps <- tryCatch(
            fread(meta_dps_file),
            error = function(e) {
                cat("  [WARN] Failed to read META_DPS file:", e$message, "\n")
                NULL
            }
        )
        if (!is.null(meta_dps)) {
            cat("  META_DPS loaded for substrate ranking:",
                nrow(meta_dps), "phosphosites\n")
        }
    } else {
        cat("  [WARN] META_DPS file not found:", basename(meta_dps_file),
            "\n         Substrates will be ordered alphabetically as fallback\n")
    }

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

        # Parse leading edge substrates (preserve original order from KSEA)
        le_str <- as.character(ksea_row$leadingEdge_substrates[1])
        le_subs <- parse_leading_edge(le_str)
        cat("    Total leading edge substrates:", length(le_subs), "\n")

        if (length(le_subs) == 0) {
            cat("    [SKIP] No leading edge substrates\n")
            next
        }

        # ------------------------------------------------------------------
        # Rank leading edge substrates by KSEA rank (Z_meta from META_DPS)
        # ------------------------------------------------------------------
        # The Z_meta from the meta-analysis represents the cross-cancer
        # combined effect size. Substrates with higher |Z_meta| are ranked
        # higher (more consistently differential across cancers).
        # The kinase z.score direction determines sorting direction:
        #   Positive kinase z.score -> rank substrates by Z_meta descending
        #   Negative kinase z.score -> rank substrates by Z_meta ascending

        kinase_z <- ksea_row$z.score[1]

        if (!is.null(meta_dps) &&
            all(c("gene_symbol_phosphosite", "Z_meta") %in% colnames(meta_dps))) {

            # Match leading edge substrates to their Z_meta values
            le_zmeta <- meta_dps[gene_symbol_phosphosite %in% le_subs,
                                 .(gene_symbol_phosphosite, Z_meta)]

            # Remove duplicates (keep first occurrence)
            le_zmeta <- le_zmeta[!duplicated(le_zmeta$gene_symbol_phosphosite)]

            if (nrow(le_zmeta) > 0) {
                # Sort by KSEA rank:
                # For activated kinase (z > 0): highest Z_meta first (most upregulated)
                # For deactivated kinase (z < 0): lowest Z_meta first (most downregulated)
                if (!is.na(kinase_z) && kinase_z > 0) {
                    le_zmeta <- le_zmeta[order(-Z_meta)]  # Descending
                } else {
                    le_zmeta <- le_zmeta[order(Z_meta)]   # Ascending
                }

                # Reorder le_subs by KSEA rank (matched substrates first,
                # then any unmatched substrates appended alphabetically)
                ranked_subs <- le_zmeta$gene_symbol_phosphosite
                unranked_subs <- sort(setdiff(le_subs, ranked_subs))
                le_subs <- c(ranked_subs, unranked_subs)

                cat("    Ranked by Z_meta (",
                    ifelse(!is.na(kinase_z) && kinase_z > 0,
                           "descending, kinase activated",
                           "ascending, kinase deactivated"),
                    ")\n", sep = "")
            } else {
                cat("    [NOTE] No Z_meta match found, using original order\n")
            }
        } else {
            cat("    [NOTE] META_DPS not available, using original order\n")
        }

        # Apply top N filter
        n_before_filter <- length(le_subs)
        if (is.finite(top_n_substrates) && length(le_subs) > top_n_substrates) {
            le_subs <- le_subs[1:top_n_substrates]
        }
        cat("    Substrates after top", min(top_n_substrates, n_before_filter),
            "filter:", length(le_subs), "/", n_before_filter, "\n")

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
        sites_found <- unique(plot_data$gene_symbol_phosphosite)
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

        # Order rows by KSEA rank (le_subs order, highest rank at top)
        # Only keep rows that exist in the matrix
        row_order <- intersect(le_subs, rownames(logfc_mat))
        logfc_mat <- logfc_mat[row_order, , drop = FALSE]

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
        # Step 3a0: Build row annotations (Hallmark gene set membership)
        # ------------------------------------------------------------------
        # Extract parent gene from substrate ID (e.g., "MCM2_S41" -> "MCM2")
        parent_genes <- sub("_[STY][0-9]+$", "", rownames(logfc_mat))

        annot_row <- data.frame(
            G2M_CHECKPOINT = ifelse(parent_genes %in% g2m_genes, "Yes", "No"),
            E2F_TARGETS    = ifelse(parent_genes %in% e2f_genes, "Yes", "No"),
            row.names = rownames(logfc_mat),
            stringsAsFactors = FALSE
        )

        # Convert to factors for consistent pheatmap color mapping
        annot_row$E2F_TARGETS    <- factor(annot_row$E2F_TARGETS,
                                            levels = c("Yes", "No"))
        annot_row$G2M_CHECKPOINT <- factor(annot_row$G2M_CHECKPOINT,
                                            levels = c("Yes", "No"))

        # Annotation color mapping
        annot_colors <- list(
            E2F_TARGETS    = c("Yes" = annot_color_E2F, "No" = "white"),
            G2M_CHECKPOINT = c("Yes" = annot_color_G2M, "No" = "white")
        )

        n_e2f_yes <- sum(annot_row$E2F_TARGETS == "Yes")
        n_g2m_yes <- sum(annot_row$G2M_CHECKPOINT == "Yes")
        cat("    E2F_TARGETS members:", n_e2f_yes, "/", nrow(logfc_mat), "\n")
        cat("    G2M_CHECKPOINT members:", n_g2m_yes, "/", nrow(logfc_mat), "\n")

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

            # Construct annotation title
            main_title <- paste0(kinase_name, " KSEA Leading Edge Top ",
                                 nrow(logfc_mat), " Substrates\n(",
                                 gsub("_", " ", comp_name), ", ", cs$label, ")")

            # Create comparison subdirectory within color scheme directory
            comp_cs_dir <- file.path(output_dir, cs_name, comp_name)
            dir.create(comp_cs_dir, recursive = TRUE, showWarnings = FALSE)

            # ------------------------------------------------------------------
            # Step 3b1: Generate the heatmap object silently
            # ------------------------------------------------------------------
            ph <- pheatmap(
                mat               = logfc_mat,
                color             = heatmap_colors,
                breaks            = breaks,
                cluster_rows      = FALSE,     # KSEA rank order
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
            # Step 3b2: Manipulate gtable to move row names to the left
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
            # Step 3b3: Draw to Output Devices
            # ------------------------------------------------------------------

            # --- Save TIFF ---
            tiff_file <- file.path(comp_cs_dir,
                                   paste0(kinase_name, "_heatmap.tiff"))
            tiff(tiff_file, width = w, height = h,
                 units = "in", res = tiff_dpi, compression = "lzw")
            grid.newpage()
            grid.draw(g)
            dev.off()

            # --- Save PDF ---
            pdf_file <- file.path(comp_cs_dir,
                                  paste0(kinase_name, "_heatmap.pdf"))
            pdf(pdf_file, width = w, height = h)
            grid.newpage()
            grid.draw(g)
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

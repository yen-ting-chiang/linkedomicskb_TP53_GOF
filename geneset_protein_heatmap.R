################################################################################
# Gene Set Protein Polar Heatmap (Circular Heatmap)
#
# Purpose:
#   Visualize protein-level DEG results (logFC) for genes within specified
#   MSigDB gene sets across multiple CPTAC datasets using a circular (polar)
#   heatmap. Each dataset is a concentric ring. Genes are arranged around the
#   circle. Color encodes logFC, and cells with adj.P.Val < 0.05 are marked
#   with a black "*".
#
# Input:
#   - protein_differential_analysis_subgroup_adjusted/{dataset}/DEG_{comparison}.csv
#   - MSigDB gene set definitions (via msigdbr package)
#
# Output:
#   - geneset_protein_heatmap/{color_scheme}/{comparison}/{gene_set_name}.tiff
#   - geneset_protein_heatmap/{color_scheme}/{comparison}/{gene_set_name}.pdf
#   - geneset_protein_heatmap/data/{comparison}/{gene_set_name}.csv
#   - geneset_protein_heatmap/data/{comparison}/{gene_set_name}.xlsx
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(msigdbr)
    library(circlize)
    library(ComplexHeatmap)
    library(grid)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Input directory (DEG results from protein_differential_analysis_subgroup_adjusted)
input_dir <- file.path(base_path, "protein_differential_analysis_subgroup_adjusted")

# Output directory
output_dir <- file.path(base_path, "geneset_protein_heatmap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Available datasets
all_datasets <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

# USER: Specify which datasets to include in the heatmap.
# Options: "BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"
# To use all available datasets, keep it as: all_datasets

#selected_datasets <- all_datasets
#selected_datasets <- c("BRCA", "COAD", "HNSCC", "LUAD", "PDAC", "UCEC") 
#selected_datasets <- c("BRCA", "GBM", "OV", "UCEC")
selected_datasets <- c("GBM", "HNSCC", "LSCC", "LUAD", "OV", "UCEC")

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

# USER: Specify which comparisons to generate heatmaps for.
# Uncomment the ones you want to include in the list below.
selected_comparisons <- c(
    #"TP53mt_vs_TP53wt",
    "MUT_GOF_vs_MUT_LOF"#,
    #"Hotspot_vs_MUT_LOF",
    #"MUT_GOF_vs_TP53wt",
    #"MUT_LOF_vs_TP53wt",
    #"Hotspot_vs_TP53wt",
    #"DN_vs_TP53wt",
    #"NonDN_vs_TP53wt",
    #"DN_vs_NonDN"
)

# ==============================================================================
# Section 1b: Clustering Settings
# ==============================================================================

# USER: Enable or disable hierarchical clustering for gene ordering (around the
#       circle) and dataset ordering (concentric ring order, outer to inner).
#       When disabled, axes fall back to alphabetical order.

cluster_genes    <- TRUE   # Cluster genes around the circle?
cluster_datasets <- TRUE   # Cluster dataset ring order?

# USER: Select distance metric for clustering.
# Options:
#   "euclidean"   : Standard Euclidean distance (default, most common)
#   "manhattan"   : Sum of absolute differences (robust to outliers)
#   "maximum"     : Maximum coordinate difference (Chebyshev distance)
#   "canberra"    : Weighted version of Manhattan (sensitive to values near 0)
#   "minkowski"   : Generalized distance (p = 2 equals Euclidean)

cluster_distance <- "euclidean"

# USER: Select linkage (agglomeration) method for clustering.
# Options:
#   "ward.D2"  : Ward's minimum variance (default, compact spherical clusters)
#   "complete" : Maximum distance between clusters (tends to produce compact clusters)
#   "average"  : Mean distance between clusters (UPGMA, balanced approach)
#   "single"   : Minimum distance (nearest-neighbor, can produce chaining)
#   "mcquitty" : Weighted average linkage (WPGMA)
#   "median"   : Median linkage (WPGMC)
#   "centroid" : Centroid linkage (UPGMC)

cluster_method <- "ward.D2"

# ==============================================================================
# Section 1c: Significance Markers
# ==============================================================================

# USER: Enable or disable significance markers ("*" for adj.P.Val < 0.05)
# Set to FALSE to hide stars, or TRUE to show them.
show_significance_stars <- FALSE

# ==============================================================================
# Section 1d: Visual Layout Parameters
# ==============================================================================

# USER: Customize track width and gap degree.
# track_width_factor: Default is "auto" (will scale based on number of datasets).
# Alternatively, set a numeric value (e.g., 0.15) to force a specific track width.
user_track_width <- "auto"

# gap_degree: Default is 55. Determines the size of the opening gap.
# user_gap_degree <- 25
user_gap_degree <- 40

# ==============================================================================
# Section 2: Define Gene Sets to Plot
# ==============================================================================

# USER: Modify this list to add/remove gene sets for heatmaps.
# Each entry: list(category, subcategory, gene_set_name)
# category and subcategory correspond to msigdbr arguments.
#
# Common categories:
#   "H"  (Hallmark)                   -> subcategory = NULL
#   "C2" (Curated gene sets)          -> subcategory = "CP:KEGG_LEGACY", "CP:REACTOME", "CGP", etc.
#   "C5" (Ontology gene sets)         -> subcategory = "GO:BP", "GO:MF", "GO:CC"
#   "C6" (Oncogenic signatures)       -> subcategory = NULL
#   "C7" (Immunologic signatures)     -> subcategory = NULL
#
# Example:
#   gene_sets_to_plot <- list(
#       list(category = "C2", subcategory = "CP:KEGG_LEGACY",
#            gene_set = "KEGG_OXIDATIVE_PHOSPHORYLATION"),
#       list(category = "H",  subcategory = NULL,
#            gene_set = "HALLMARK_E2F_TARGETS"),
#       list(category = "C2", subcategory = "CP:REACTOME",
#            gene_set = "REACTOME_CELL_CYCLE")
#   )

gene_sets_to_plot <- list(
    list(category = "C2", subcategory = "CP:KEGG_LEGACY",
         gene_set = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES")
)

# ==============================================================================
# Section 3: Define Color Palettes
# ==============================================================================
# Four color schemes suitable for diverging logFC heatmaps.
# All selected schemes will be generated. USER can select preferred one afterward.

color_schemes <- list(
    # Palette 1: Red-Blue (RdBu reversed) - Nature, Cell common
    "RdBu" = list(
        name   = "RdBu",
        label  = "Red-Blue (Nature/Cell style)",
        low    = "#2166AC",   # Blue (negative logFC)
        mid    = "#F7F7F7",   # White (zero)
        high   = "#B2182B"    # Red (positive logFC)
    ),
    # Palette 2: Red-Yellow-Blue (RdYlBu reversed) - Science common
    "RdYlBu" = list(
        name   = "RdYlBu",
        label  = "Red-Yellow-Blue (Science style)",
        low    = "#313695",   # Dark blue
        mid    = "#FFFFBF",   # Yellow
        high   = "#A50026"    # Dark red
    ),
    # Palette 3: Purple-Orange (PuOr) - colorblind-friendly
    "PuOr" = list(
        name   = "PuOr",
        label  = "Purple-Orange (colorblind-friendly)",
        low    = "#542788",   # Purple (negative)
        mid    = "#F7F7F7",   # White
        high   = "#E08214"    # Orange (positive)
    ),
    # Palette 4: Teal-Red (custom) - high-impact journal style
    "TealRed" = list(
        name   = "TealRed",
        label  = "Teal-Red (Lancet/NEJM style)",
        low    = "#008080",   # Teal (negative)
        mid    = "#F5F5F5",   # Near-white
        high   = "#CD2626"    # Red (positive)
    )
)

# USER: Specify which color schemes to generate.
# Select one or more from: "RdBu", "RdYlBu", "PuOr", "TealRed"
# To generate all schemes, keep it as: names(color_schemes)

selected_color_schemes <- names(color_schemes)

cat("====================================================================\n")
cat("Gene Set Protein Polar Heatmap (Circular Heatmap)\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 4: Retrieve Gene Set Definitions from MSigDB
# ==============================================================================

cat("Retrieving gene set definitions from MSigDB ...\n")

# Collect all unique gene sets to query
all_gene_set_info <- list()
for (gs_entry in gene_sets_to_plot) {
    key <- paste0(gs_entry$category, "|", ifelse(is.null(gs_entry$subcategory), "NULL", gs_entry$subcategory))
    if (is.null(all_gene_set_info[[key]])) {
        all_gene_set_info[[key]] <- list(
            category    = gs_entry$category,
            subcategory = gs_entry$subcategory,
            gene_sets   = c()
        )
    }
    all_gene_set_info[[key]]$gene_sets <- c(all_gene_set_info[[key]]$gene_sets, gs_entry$gene_set)
}

# Query msigdbr and build gene set -> gene symbol mapping
gene_set_genes <- list()  # gene_set_name -> character vector of gene symbols

for (key in names(all_gene_set_info)) {
    info <- all_gene_set_info[[key]]
    cat("  Querying MSigDB: category =", info$category)
    if (!is.null(info$subcategory)) cat(", subcategory =", info$subcategory)
    cat("\n")

    msig_df <- tryCatch({
        if (is.null(info$subcategory)) {
            msigdbr(species = "Homo sapiens", category = info$category)
        } else {
            msigdbr(species = "Homo sapiens", category = info$category,
                    subcategory = info$subcategory)
        }
    }, error = function(e) {
        cat("    [ERROR] msigdbr query failed:", e$message, "\n")
        NULL
    })

    if (is.null(msig_df)) next

    for (gs_name in info$gene_sets) {
        genes <- msig_df %>%
            filter(gs_name == !!gs_name) %>%
            pull(gene_symbol) %>%
            unique() %>%
            sort()

        if (length(genes) == 0) {
            cat("    [WARNING] Gene set not found:", gs_name, "\n")
        } else {
            gene_set_genes[[gs_name]] <- genes
            cat("    Gene set:", gs_name, "->", length(genes), "genes\n")
        }
    }
}

if (length(gene_set_genes) == 0) {
    stop("No valid gene sets found. Please check gene_sets_to_plot definitions.")
}

cat("\n")

# ==============================================================================
# Section 5: Helper Function for Drawing Circular Heatmap
# ==============================================================================

#' Draw a polar (circular) heatmap with dendrograms and significance markers
#'
#' @param logfc_mat Numeric matrix (genes x datasets) of logFC values
#' @param sig_mat   Logical matrix (genes x datasets) TRUE = adj.P.Val < 0.05
#' @param col_fun   colorRamp2 color mapping function
#' @param na_col    Color for NA cells
#' @param title_main Main title (comparison name)
#' @param title_sub  Subtitle (gene set name)
#' @param gene_dend  Dendrogram object for genes (NULL if no clustering)
#' @param show_dendro Whether to display the dendrogram
draw_polar_heatmap <- function(logfc_mat, sig_mat, col_fun, na_col = "grey90",
                                title_main = "", title_sub = "",
                                gene_dend = NULL, show_dendro = TRUE,
                                show_stars = FALSE,
                                gap_deg = 55, track_w = "auto") {

    n_genes    <- nrow(logfc_mat)
    n_datasets <- ncol(logfc_mat)

    # --- Adaptive sizing parameters ---
    gene_label_cex <- if (n_genes > 100) 0.50
                      else if (n_genes > 60) 0.65
                      else if (n_genes > 30) 0.80
                      else 1.00

    sig_cex <- if (n_genes > 100) 0.35
               else if (n_genes > 60) 0.45
               else 0.55

    # --- Dynamic height budget (total must stay under 1.0) ---
    # Reserve space for track margins: each track gets 0.01 (top + bottom)
    # Total tracks = 1 (labels) + n_datasets (heatmap) + 1 (dendro) = n_datasets + 2
    n_tracks     <- n_datasets + 2
    margin_total <- n_tracks * 0.01
    budget       <- 0.92 - margin_total   # usable radius after margins

    # Gene label track height (needs space for rotated text)
    label_track_h <- if (n_genes > 100) 0.12
                     else if (n_genes > 60) 0.15
                     else 0.20

    # Dendrogram track height
    dendro_track_h <- 0.10

    # Remaining budget for dataset tracks
    ds_budget <- budget - label_track_h - dendro_track_h
    track_h   <- if (track_w == "auto") {
        max(0.02, ds_budget / n_datasets)
    } else {
        as.numeric(track_w)
    }

    # --- Initialize circos ---
    circos.clear()
    circos.par(
        start.degree    = 135,
        gap.degree      = gap_deg,
        cell.padding    = c(0, 0, 0, 0),
        track.margin    = c(0.005, 0.005),
        canvas.xlim     = c(-1.25, 1.25),
        canvas.ylim     = c(-1.25, 1.25)
    )

    circos.initialize(sectors = "heatmap", xlim = c(0.5, n_genes + 0.5))

    # --- Track 1 (outermost): Gene labels ---
    circos.track(
        sectors      = "heatmap",
        ylim         = c(0, 1),
        track.height = label_track_h,
        bg.border    = NA,
        panel.fun    = function(x, y) {
            for (i in 1:n_genes) {
                circos.text(i, 0.5, rownames(logfc_mat)[i],
                    facing    = "clockwise",
                    niceFacing = TRUE,
                    adj       = c(0, 0.5),
                    cex       = gene_label_cex,
                    font      = 2,    # bold
                    col       = "black")
            }
        }
    )

    # --- Heatmap tracks (one concentric ring per dataset, outer to inner) ---
    R_ds <- numeric(n_datasets)
    ds_names <- colnames(logfc_mat)
    for (j in 1:n_datasets) {
        local({
            jj <- j
            circos.track(
                sectors      = "heatmap",
                ylim         = c(0, 1),
                track.height = track_h,
                bg.border    = NA,
                panel.fun    = function(x, y) {
                    for (i in 1:n_genes) {
                        val    <- logfc_mat[i, jj]
                        is_sig <- sig_mat[i, jj]

                        cell_col <- if (is.na(val)) na_col else col_fun(val)

                        circos.rect(i - 0.5, 0, i + 0.5, 1,
                            col = cell_col, border = NA)

                        if (show_stars && !is.na(is_sig) && is_sig) {
                            circos.text(i, 0.5, "*",
                                cex    = sig_cex,
                                col    = "black",
                                facing = "inside")
                        }
                    }

                    # Add dataset track name natively at the start of the plot (135 degrees)
                    circos.text(
                        x = 0.5,             # Left edge of the track
                        y = 0.5,             # Middle of the track radially
                        labels = paste0(ds_names[jj], "  "), # Add trailing spaces to avoid touching colors
                        facing = "bending.inside", # Follow the curvature of the track
                        adj = c(1, 0.5),     # Right-justify (pushes text into the gap)
                        cex = if (n_datasets > 8) 0.6 else 0.8,
                        col = "black"
                    )
                }
            )
            r_top <- get.cell.meta.data("cell.top.radius", sector.index = "heatmap", track.index = jj + 1)
            r_bot <- get.cell.meta.data("cell.bottom.radius", sector.index = "heatmap", track.index = jj + 1)
            R_ds[jj] <<- (r_top + r_bot) / 2
        })
    }

    # --- Dendrogram track (innermost) ---
    if (show_dendro && !is.null(gene_dend) && n_genes >= 2) {
        dend_height <- attr(gene_dend, "height")
        circos.track(
            sectors      = "heatmap",
            ylim         = c(0, dend_height),
            track.height = dendro_track_h,
            bg.border    = NA,
            panel.fun    = function(x, y) {
                circos.dendrogram(gene_dend, facing = "outside",
                                  max_height = dend_height)
            }
        )
    }

    circos.clear()

    # --- Annotations using grid (overlaid on the circos canvas) ---

    # Title and subtitle (top center)
    grid.text(title_main,
        x  = unit(0.5, "npc"), y = unit(0.97, "npc"),
        gp = gpar(fontsize = 12, fontface = "bold", col = "black"))
    grid.text(title_sub,
        x  = unit(0.5, "npc"), y = unit(0.94, "npc"),
        gp = gpar(fontsize = 9, col = "grey30"))

    # --- Legends (Bottom) ---
    lgd_logfc <- Legend(
        col_fun       = col_fun,
        title         = "logFC",
        direction     = "horizontal",
        title_gp      = gpar(fontsize = 9, fontface = "bold"),
        labels_gp     = gpar(fontsize = 8),
        legend_width  = unit(4, "cm"),
        grid_height   = unit(3, "mm"),
        border        = FALSE
    )
    
    lgd_list <- list(lgd_logfc)
    
    if (show_stars) {
        lgd_sig <- Legend(
            labels    = "adj.P.Val < 0.05",
            title     = "Significance",
            type      = "points",
            pch       = 8,
            size      = unit(3, "mm"),
            legend_gp = gpar(col = "black"),
            title_gp  = gpar(fontsize = 9, fontface = "bold"),
            labels_gp = gpar(fontsize = 8),
            border    = FALSE,
            direction = "horizontal"
        )
        lgd_list <- c(lgd_list, list(lgd_sig))
    }
    
    lgd_na <- Legend(
        labels    = "not measured",
        title     = "NA",
        legend_gp = gpar(fill = na_col, col = "black"),
        title_gp  = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        border    = TRUE,
        direction = "horizontal"
    )
    lgd_list <- c(lgd_list, list(lgd_na))
    
    pd <- packLegend(list = lgd_list, direction = "horizontal", gap = unit(8, "mm"))
    draw(pd, x = unit(0.5, "npc"), y = unit(0.02, "npc"), just = c("center", "bottom"))
}

# ==============================================================================
# Section 6: Load DEG Data and Generate Polar Heatmaps
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

    # --- Load DEG data from all selected datasets for this comparison ---
    deg_all <- list()

    for (ds in selected_datasets) {
        deg_file <- file.path(input_dir, ds, paste0("DEG_", comp_name, ".csv"))
        if (!file.exists(deg_file)) {
            cat("  [SKIP] DEG file not found for", ds, ":", basename(deg_file), "\n")
            next
        }

        deg_df <- tryCatch(
            fread(deg_file),
            error = function(e) {
                cat("  [ERROR] Failed to read", ds, ":", e$message, "\n")
                NULL
            }
        )

        if (is.null(deg_df) || nrow(deg_df) == 0) next

        # Ensure required columns exist
        required_cols <- c("gene_symbol", "logFC", "adj.P.Val")
        if (!all(required_cols %in% colnames(deg_df))) {
            cat("  [SKIP]", ds, ": missing required columns\n")
            next
        }

        # Add dataset column
        deg_df$dataset <- ds

        deg_all[[ds]] <- deg_df
    }

    if (length(deg_all) == 0) {
        cat("  [SKIP] No DEG data available for any dataset\n\n")
        next
    }

    # Combine all datasets
    deg_combined <- rbindlist(deg_all, fill = TRUE)

    # --- Generate heatmap for each gene set ---
    for (gs_name in names(gene_set_genes)) {
        gs_genes <- gene_set_genes[[gs_name]]

        cat("\n  Gene set:", gs_name, "(", length(gs_genes), "genes )\n")

        # Filter DEG data to genes in this gene set
        plot_data <- deg_combined[gene_symbol %in% gs_genes]

        if (nrow(plot_data) == 0) {
            cat("    [SKIP] No matching genes found in DEG data\n")
            next
        }

        # Keep only necessary columns
        plot_data <- plot_data[, .(gene_symbol, dataset, logFC, adj.P.Val)]

        # Handle duplicate gene_symbol per dataset (keep first occurrence)
        plot_data <- plot_data[!duplicated(plot_data[, .(gene_symbol, dataset)])]

        # Report coverage
        genes_found <- sort(unique(plot_data$gene_symbol))
        datasets_found <- sort(unique(plot_data$dataset))
        cat("    Genes found:", length(genes_found), "/", length(gs_genes), "\n")
        cat("    Datasets found:", length(datasets_found), "\n")

        # Create wide matrices for logFC and adj.P.Val
        logfc_wide <- dcast(plot_data, gene_symbol ~ dataset, value.var = "logFC")
        pval_wide  <- dcast(plot_data, gene_symbol ~ dataset, value.var = "adj.P.Val")

        # --- Build numeric matrices ---
        logfc_mat <- as.data.frame(logfc_wide)
        rownames(logfc_mat) <- logfc_mat$gene_symbol
        logfc_mat$gene_symbol <- NULL
        logfc_mat <- as.matrix(logfc_mat)

        pval_mat <- as.data.frame(pval_wide)
        rownames(pval_mat) <- pval_mat$gene_symbol
        pval_mat$gene_symbol <- NULL
        pval_mat <- as.matrix(pval_mat)

        # Significance matrix (TRUE if adj.P.Val < 0.05)
        sig_mat <- !is.na(pval_mat) & pval_mat < 0.05

        # Impute NA with 0 for distance calculation (NA = gene not measured)
        logfc_mat_imputed <- logfc_mat
        logfc_mat_imputed[is.na(logfc_mat_imputed)] <- 0

        # --- Gene clustering (determines order around the circle) ---
        gene_dend <- NULL
        if (cluster_genes && nrow(logfc_mat_imputed) >= 2) {
            gene_dist   <- dist(logfc_mat_imputed, method = cluster_distance)
            gene_hclust <- hclust(gene_dist, method = cluster_method)
            gene_order  <- rownames(logfc_mat_imputed)[gene_hclust$order]
            gene_dend   <- as.dendrogram(gene_hclust)
            cat("    Gene clustering: distance =", cluster_distance,
                ", method =", cluster_method, "\n")
        } else {
            gene_order <- sort(rownames(logfc_mat))
            if (cluster_genes && nrow(logfc_mat_imputed) < 2) {
                cat("    [NOTE] Too few genes for clustering, using alphabetical order\n")
            }
        }

        # --- Dataset clustering (determines track order, outer to inner) ---
        ds_hclust <- NULL
        if (cluster_datasets && ncol(logfc_mat_imputed) >= 2) {
            ds_dist   <- dist(t(logfc_mat_imputed), method = cluster_distance)
            ds_hclust <- hclust(ds_dist, method = cluster_method)
            ds_order  <- colnames(logfc_mat_imputed)[ds_hclust$order]
            cat("    Dataset clustering: distance =", cluster_distance,
                ", method =", cluster_method, "\n")
        } else {
            ds_order <- sort(colnames(logfc_mat))
            if (cluster_datasets && ncol(logfc_mat_imputed) < 2) {
                cat("    [NOTE] Too few datasets for clustering, using alphabetical order\n")
            }
        }

        # --- Reorder matrices by clustering results ---
        logfc_mat_plot <- logfc_mat[gene_order, ds_order, drop = FALSE]
        sig_mat_plot   <- sig_mat[gene_order, ds_order, drop = FALSE]

        # Determine symmetric color scale limits
        logfc_max <- max(abs(logfc_mat_plot), na.rm = TRUE)
        logfc_lim <- ceiling(logfc_max * 10) / 10  # round up to 1 decimal

        # --- Save data CSV and XLSX ---
        comp_data_dir <- file.path(data_out_dir, comp_name)
        dir.create(comp_data_dir, recursive = TRUE, showWarnings = FALSE)

        # Prepare export data: wide format with logFC and significance
        export_data <- merge(logfc_wide, pval_wide,
                              by = "gene_symbol", suffixes = c("_logFC", "_adjPVal"))
        export_data <- export_data[order(export_data$gene_symbol), ]

        # Save CSV
        data_csv <- file.path(comp_data_dir, paste0(gs_name, ".csv"))
        fwrite(export_data, data_csv)

        # Save XLSX
        data_xlsx <- file.path(comp_data_dir, paste0(gs_name, ".xlsx"))
        wb_data <- createWorkbook()
        addWorksheet(wb_data, "Heatmap_Data")
        writeData(wb_data, 1, export_data)

        header_style <- createStyle(
            fontName = "Arial", fontSize = 7,
            textDecoration = "bold", halign = "center", valign = "center",
            border = "bottom", fgFill = "#4472C4", fontColour = "white"
        )
        body_style <- createStyle(
            fontName = "Arial", fontSize = 7,
            halign = "left", valign = "center"
        )
        addStyle(wb_data, 1, header_style,
                 rows = 1, cols = 1:ncol(export_data), gridExpand = TRUE)
        if (nrow(export_data) > 0) {
            addStyle(wb_data, 1, body_style,
                     rows = 2:(nrow(export_data) + 1), cols = 1:ncol(export_data),
                     gridExpand = TRUE, stack = FALSE)
        }
        setColWidths(wb_data, 1, cols = 1:ncol(export_data), widths = "auto")
        saveWorkbook(wb_data, data_xlsx, overwrite = TRUE)

        cat("    Data saved:", basename(data_csv), "&", basename(data_xlsx), "\n")

        # --- Generate polar heatmaps for each selected color scheme ---
        for (cs_name in selected_color_schemes) {
            cs <- color_schemes[[cs_name]]

            # Color mapping function (circlize)
            col_fun <- colorRamp2(
                c(-logfc_lim, 0, logfc_lim),
                c(cs$low, cs$mid, cs$high)
            )

            # Plot dimensions (square for circular plots)
            n_genes_plot    <- nrow(logfc_mat_plot)
            n_datasets_plot <- ncol(logfc_mat_plot)
            plot_dim <- max(8, 6 + n_genes_plot * 0.015 + n_datasets_plot * 0.1)
            plot_dim <- min(plot_dim, 14)  # cap at 14 inches

            # Create comparison subdirectory within color scheme directory
            comp_cs_dir <- file.path(output_dir, cs_name, comp_name)
            dir.create(comp_cs_dir, recursive = TRUE, showWarnings = FALSE)

            # --- Save TIFF ---
            tiff_file <- file.path(comp_cs_dir, paste0(gs_name, ".tiff"))
            tiff(tiff_file, width = plot_dim, height = plot_dim,
                 units = "in", res = 300, compression = "lzw")

            draw_polar_heatmap(
                logfc_mat  = logfc_mat_plot,
                sig_mat    = sig_mat_plot,
                col_fun    = col_fun,
                na_col     = "grey90",
                title_main = gsub("_", " ", comp_name),
                title_sub  = paste0(gsub("_", " ", gs_name), " (", cs$label, ")"),
                gene_dend  = gene_dend,
                show_dendro = cluster_genes,
                show_stars  = show_significance_stars,
                gap_deg     = user_gap_degree,
                track_w     = user_track_width
            )

            dev.off()

            # --- Save PDF ---
            pdf_file <- file.path(comp_cs_dir, paste0(gs_name, ".pdf"))
            pdf(pdf_file, width = plot_dim, height = plot_dim)

            draw_polar_heatmap(
                logfc_mat  = logfc_mat_plot,
                sig_mat    = sig_mat_plot,
                col_fun    = col_fun,
                na_col     = "grey90",
                title_main = gsub("_", " ", comp_name),
                title_sub  = paste0(gsub("_", " ", gs_name), " (", cs$label, ")"),
                gene_dend  = gene_dend,
                show_dendro = cluster_genes,
                show_stars  = show_significance_stars,
                gap_deg     = user_gap_degree,
                track_w     = user_track_width
            )

            dev.off()
        }

        cat("    Polar heatmaps saved for selected color schemes\n")
    }

    cat("\n")
}

cat("====================================================================\n")
cat("Gene set protein polar heatmaps completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

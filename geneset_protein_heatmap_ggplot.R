################################################################################
# Gene Set Protein Heatmap
#
# Purpose:
#   Visualize protein-level DEG results (logFC) for genes within specified
#   MSigDB gene sets across multiple CPTAC datasets using heatmaps.
#   Color encodes logFC (fold change direction/magnitude).
#   Cells with adj.P.Val < 0.05 are marked with a black "*".
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
  library(ggdendro)
  library(patchwork)
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

selected_datasets <- all_datasets

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
# Options: "TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
#          "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
#          "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"
# To generate all comparisons, keep it as: comparisons$name

selected_comparisons <- comparisons$name

# ==============================================================================
# Section 1b: Clustering Settings
# ==============================================================================

# USER: Enable or disable hierarchical clustering for x-axis (datasets) and
#       y-axis (genes). When disabled, axes fall back to alphabetical order.

cluster_genes    <- TRUE   # Cluster y-axis (genes)?
cluster_datasets <- TRUE   # Cluster x-axis (datasets)?

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
       gene_set = "KEGG_OXIDATIVE_PHOSPHORYLATION")
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
cat("Gene Set Protein Heatmap\n")
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
# Section 5: Load DEG Data and Generate Heatmaps
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

for (j in seq_len(nrow(comparisons))) {
  comp_name  <- comparisons$name[j]
  comp_label <- comparisons$label[j]
  
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
    
    # Create wide matrices for logFC and significance
    logfc_wide <- dcast(plot_data, gene_symbol ~ dataset, value.var = "logFC")
    pval_wide  <- dcast(plot_data, gene_symbol ~ dataset, value.var = "adj.P.Val")
    
    # Create significance label ("*" if adj.P.Val < 0.05)
    plot_data$sig_label <- ifelse(!is.na(plot_data$adj.P.Val) & plot_data$adj.P.Val < 0.05,
                                  "*", "")
    
    # --- Determine axis ordering via hierarchical clustering or alphabetical ---
    # Build numeric matrix for clustering (genes as rows, datasets as columns)
    logfc_mat <- as.data.frame(logfc_wide)
    rownames(logfc_mat) <- logfc_mat$gene_symbol
    logfc_mat$gene_symbol <- NULL
    logfc_mat <- as.matrix(logfc_mat)
    
    # Impute NA with 0 for distance calculation (NA = gene not measured)
    logfc_mat_imputed <- logfc_mat
    logfc_mat_imputed[is.na(logfc_mat_imputed)] <- 0
    
    # Y-axis: gene ordering
    if (cluster_genes && nrow(logfc_mat_imputed) >= 2) {
      gene_dist <- dist(logfc_mat_imputed, method = cluster_distance)
      gene_hclust <- hclust(gene_dist, method = cluster_method)
      gene_order <- rownames(logfc_mat_imputed)[gene_hclust$order]
      cat("    Gene clustering: distance =", cluster_distance,
          ", method =", cluster_method, "\n")
    } else {
      gene_order <- sort(unique(plot_data$gene_symbol), decreasing = TRUE)
      if (cluster_genes && nrow(logfc_mat_imputed) < 2) {
        cat("    [NOTE] Too few genes for clustering, using alphabetical order\n")
      }
    }
    
    # X-axis: dataset ordering
    if (cluster_datasets && ncol(logfc_mat_imputed) >= 2) {
      ds_dist <- dist(t(logfc_mat_imputed), method = cluster_distance)
      ds_hclust <- hclust(ds_dist, method = cluster_method)
      ds_order <- colnames(logfc_mat_imputed)[ds_hclust$order]
      cat("    Dataset clustering: distance =", cluster_distance,
          ", method =", cluster_method, "\n")
    } else {
      ds_order <- sort(unique(as.character(plot_data$dataset)))
      if (cluster_datasets && ncol(logfc_mat_imputed) < 2) {
        cat("    [NOTE] Too few datasets for clustering, using alphabetical order\n")
      }
    }
    
    # Apply ordering
    plot_data$gene_symbol <- factor(plot_data$gene_symbol, levels = gene_order)
    plot_data$dataset <- factor(plot_data$dataset, levels = ds_order)
    
    # Determine symmetric color scale limits
    logfc_max <- max(abs(plot_data$logFC), na.rm = TRUE)
    logfc_lim <- ceiling(logfc_max * 10) / 10  # round up to 1 decimal
    
    # --- Save data CSV and XLSX ---
    comp_data_dir <- file.path(data_out_dir, comp_name)
    dir.create(comp_data_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Prepare export data: wide format with logFC and significance
    export_data <- merge(logfc_wide, pval_wide,
                         by = "gene_symbol", suffixes = c("_logFC", "_adjPVal"))
    # Sort by gene_symbol
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
    
    # --- Generate heatmaps for each selected color scheme ---
    for (cs_name in selected_color_schemes) {
      cs <- color_schemes[[cs_name]]
      
      # Build the main heatmap
      p_heat <- ggplot(plot_data, aes(x = dataset, y = gene_symbol)) +
        # Heatmap tiles
        geom_tile(aes(fill = logFC), color = NA) +
        # Significance markers
        geom_text(aes(label = sig_label),
                  color = "black", size = 4, vjust = 0.75, fontface = "bold") +
        # Color scale
        scale_fill_gradient2(
          low = cs$low, mid = cs$mid, high = cs$high,
          midpoint = 0,
          limits = c(-logfc_lim, logfc_lim),
          name = "logFC",
          na.value = "grey90"
        ) +
        labs(
          x = NULL,
          y = NULL,
          title = paste0(gsub("_", " ", comp_name)),
          subtitle = paste0(gsub("_", " ", gs_name))
        ) +
        theme(
          # White background, no outer border
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border     = element_blank(),
          panel.grid       = element_blank(),
          # Axes
          axis.line        = element_blank(),
          axis.ticks       = element_line(color = "black", linewidth = 0.3),
          axis.text.x      = element_text(color = "black", size = 8,
                                          angle = 45, hjust = 1, vjust = 1),
          axis.text.y      = element_text(color = "black", size = 7,
                                          face = "italic"),
          # Title
          plot.title       = element_text(hjust = 0.5, face = "bold",
                                          size = 10, color = "black"),
          plot.subtitle    = element_text(hjust = 0.5, face = "plain",
                                          size = 8, color = "grey30"),
          # Legend
          legend.background = element_rect(fill = "white", color = NA),
          legend.key        = element_rect(fill = "white", color = NA),
          legend.text       = element_text(size = 7),
          legend.title      = element_text(size = 8, face = "bold"),
          legend.position   = "right",
          # Margin
          plot.margin      = margin(2, 2, 2, 2)
        )
      
      # --- Build dendrograms ---
      has_top_dendro  <- cluster_datasets && exists("ds_hclust")
      has_left_dendro <- cluster_genes && exists("gene_hclust")
      
      if (has_top_dendro) {
        # Top dendrogram (datasets): horizontal
        ddata_top <- dendro_data(as.dendrogram(ds_hclust), type = "rectangle")
        p_dendro_top <- ggplot(segment(ddata_top)) +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                       color = "black", linewidth = 0.4) +
          scale_x_continuous(expand = c(0.5 / ncol(logfc_mat_imputed), 0.5 / ncol(logfc_mat_imputed))) +
          scale_y_reverse(expand = expansion(mult = c(0.05, 0))) +
          theme_void() +
          theme(
            plot.background = element_rect(fill = "white", color = NA),
            plot.margin = margin(5, 2, 0, 2)
          )
      }
      
      if (has_left_dendro) {
        # Left dendrogram (genes): vertical, rotated
        ddata_left <- dendro_data(as.dendrogram(gene_hclust), type = "rectangle")
        p_dendro_left <- ggplot(segment(ddata_left)) +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                       color = "black", linewidth = 0.4) +
          scale_x_continuous(expand = c(0.5 / nrow(logfc_mat_imputed), 0.5 / nrow(logfc_mat_imputed))) +
          scale_y_reverse(expand = expansion(mult = c(0.05, 0))) +
          coord_flip() +
          theme_void() +
          theme(
            plot.background = element_rect(fill = "white", color = NA),
            plot.margin = margin(2, 0, 2, 5)
          )
      }
      
      # --- Compose layout with patchwork ---
      # Dendrogram width/height proportions
      dendro_ratio <- 0.12  # fraction of total for dendrogram panels
      
      n_datasets <- length(unique(plot_data$dataset))
      n_genes    <- length(unique(plot_data$gene_symbol))
      
      if (has_top_dendro && has_left_dendro) {
        # Both dendrograms: 2x2 grid
        # Top-left: empty spacer, Top-right: top dendrogram
        # Bottom-left: left dendrogram, Bottom-right: heatmap
        p_spacer <- plot_spacer() +
          theme(plot.background = element_rect(fill = "white", color = NA))
        
        p_combined <- (p_spacer + p_dendro_top) /
          (p_dendro_left + p_heat) +
          plot_layout(
            widths  = c(dendro_ratio, 1 - dendro_ratio),
            heights = c(dendro_ratio, 1 - dendro_ratio)
          )
        
        plot_width  <- max(5.5, 2.0 + n_datasets * 0.55 + 1.5)
        plot_height <- max(5.0, 1.5 + n_genes * 0.18 + 1.5)
        
      } else if (has_top_dendro) {
        # Top dendrogram only
        p_combined <- p_dendro_top / p_heat +
          plot_layout(heights = c(dendro_ratio, 1 - dendro_ratio))
        
        plot_width  <- max(4.5, 1.5 + n_datasets * 0.55 + 1.5)
        plot_height <- max(4.5, 1.5 + n_genes * 0.18 + 1.2)
        
      } else if (has_left_dendro) {
        # Left dendrogram only
        p_combined <- p_dendro_left + p_heat +
          plot_layout(widths = c(dendro_ratio, 1 - dendro_ratio))
        
        plot_width  <- max(5.0, 1.8 + n_datasets * 0.55 + 1.5)
        plot_height <- max(4.0, 1.2 + n_genes * 0.18 + 1.0)
        
      } else {
        # No dendrograms
        p_combined <- p_heat
        
        plot_width  <- max(4.5, 1.5 + n_datasets * 0.55 + 1.5)
        plot_height <- max(4.0, 1.2 + n_genes * 0.18 + 1.0)
      }
      
      # Create comparison subdirectory within color scheme directory
      comp_cs_dir <- file.path(output_dir, cs_name, comp_name)
      dir.create(comp_cs_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Save TIFF
      tiff_file <- file.path(comp_cs_dir, paste0(gs_name, ".tiff"))
      tiff(tiff_file, width = plot_width, height = plot_height,
           units = "in", res = 300, compression = "lzw")
      print(p_combined)
      dev.off()
      
      # Save PDF
      pdf_file <- file.path(comp_cs_dir, paste0(gs_name, ".pdf"))
      pdf(pdf_file, width = plot_width, height = plot_height)
      print(p_combined)
      dev.off()
    }
    
    cat("    Heatmaps saved for selected color schemes\n")
  }
  
  cat("\n")
}

cat("====================================================================\n")
cat("Gene set protein heatmaps completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

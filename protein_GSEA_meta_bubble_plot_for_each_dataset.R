################################################################################
# Protein GSEA Meta Bubble Plot for Each Dataset
#
# Purpose:
#   Visualize protein-level GSEA results across datasets for specific pathways
#   using bubble plots. All specified pathways (across collections) are combined
#   into a single plot per comparison.
#   Color encodes NES (enrichment direction/magnitude).
#   Size encodes -log10(p-value) (statistical significance).
#   Bubbles with pval < 0.05 have a thick black border.
#
# Input:
#   - mRNA_GSEA_meta vs_protein_GSEA_meta/{collection}/META_GSEA_{comparison}.csv
#
# Output:
#   - mRNA_GSEA_meta vs_protein_GSEA_meta/bubble_plot/
#       {color_scheme}/{comparison}.tiff and .pdf
#       data/{comparison}.csv and .xlsx
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Input directory (integration output)
input_dir <- file.path(base_path, "mRNA_GSEA_meta vs_protein_GSEA_meta")

# Output directory
output_dir <- file.path(input_dir, "bubble_plot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

# USER: Specify which comparisons to generate bubble plots for.
# Options: "TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
#          "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
#          "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"
# To generate all comparisons, keep it as: comparisons$name

# selected_comparisons <- comparisons$name
selected_comparisons <- c("MUT_GOF_vs_MUT_LOF")


# ==============================================================================
# Section 2: Define Pathways to Plot
# ==============================================================================
# USER: modify this list to add/remove pathways for bubble plots.
# Each entry: collection folder name -> vector of pathway names.
# All pathways across collections will be combined into ONE bubble plot.

pathways_to_plot <- list(
    "Hallmark" = c(
        "HALLMARK_E2F_TARGETS",
        #"HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_G2M_CHECKPOINT",
        #"HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_MTORC1_SIGNALING"#,
        #"HALLMARK_ESTROGEN_RESPONSE_LATE",
        #"HALLMARK_ESTROGEN_RESPONSE_EARLY",
        #"HALLMARK_MYOGENESIS"
    ),
    "C2_CP_KEGG_LEGACY" = c(
        "KEGG_DNA_REPLICATION",
        "KEGG_CELL_CYCLE",
        #"KEGG_MISMATCH_REPAIR",
        "KEGG_OXIDATIVE_PHOSPHORYLATION",
        "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"
    ),
    "C2_CGP" = c(
        #"FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS"
    )
)

# ==============================================================================
# Section 3: Define Color Palettes
# ==============================================================================
# All 4 color schemes will be generated. USER can select preferred one afterward.

color_schemes <- list(
    # Palette 1: Red-Blue (RdBu reversed) - Nature, Cell common
    "RdBu" = list(
        name   = "RdBu",
        label  = "Red-Blue (Nature/Cell style)",
        low    = "#2166AC",   # Blue (negative NES)
        mid    = "#F7F7F7",   # White (zero)
        high   = "#B2182B"    # Red (positive NES)
    ),
    # Palette 2: Red-Yellow-Blue (RdYlBu reversed) - Science common
    "RdYlBu" = list(
        name   = "RdYlBu",
        label  = "Red-Yellow-Blue (Science style)",
        low    = "#313695",   # Dark blue
        mid    = "#FFFFBF",   # Yellow
        high   = "#A50026"    # Dark red
    ),
    # Palette 3: Purple-Green (PRGn reversed) - colorblind-friendly
    "PRGn" = list(
        name   = "PRGn",
        label  = "Purple-Green (colorblind-friendly)",
        low    = "#762A83",   # Purple (negative)
        mid    = "#F7F7F7",   # White
        high   = "#1B7837"    # Green (positive)
    ),
    # Palette 4: Red-Yellow-Green (RdYlGn reversed)
    "RdYlGn" = list(
        name   = "RdYlGn",
        label  = "Red-Yellow-Green",
        low    = "#D73027",   # Red (negative)
        mid    = "#FFFFBF",   # Yellow
        high   = "#1A9850"    # Green (positive)
    )
)
# USER: Specify which color schemes to generate.
# Select one or more from: "RdBu", "RdYlBu", "PRGn", "RdYlGn"
# To generate all schemes, keep it as: names(color_schemes)

# selected_color_schemes <- names(color_schemes)
selected_color_schemes <- c("RdBu", "PRGn")

cat("====================================================================\n")
cat("Protein GSEA Meta Bubble Plot for Each Dataset\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 4: Data Extraction and Bubble Plot Generation
# ==============================================================================

# Filter and validate selected comparisons
comparisons <- comparisons[comparisons$name %in% selected_comparisons, , drop = FALSE]
if (nrow(comparisons) == 0) {
    stop("No valid comparisons selected. Please choose from: ", paste(c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF", "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt", "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"), collapse = ", "))
}

# Validate selected color schemes
selected_color_schemes <- intersect(selected_color_schemes, names(color_schemes))
if (length(selected_color_schemes) == 0) {
    stop("No valid color schemes selected. Please choose from: ", paste(names(color_schemes), collapse = ", "))
}

# Create output subdirectories for each selected color scheme and data
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

    # --- Load and combine data from all specified collections ---
    all_long_list <- list()

    for (coll_name in names(pathways_to_plot)) {
        target_pathways <- pathways_to_plot[[coll_name]]

        int_file <- file.path(input_dir, coll_name, paste0("META_GSEA_", comp_name, ".csv"))
        if (!file.exists(int_file)) {
            cat("  [SKIP] File not found for collection:", coll_name, "\n")
            next
        }

        full_data <- fread(int_file)

        # Filter to target pathways
        plot_data <- full_data[pathway %in% target_pathways]
        if (nrow(plot_data) == 0) {
            cat("  [SKIP] No target pathways in:", coll_name, "\n")
            next
        }

        # Retrieve Z_meta_protein for y-axis ordering
        z_meta_prot <- plot_data[, .(pathway, Z_meta_protein)]

        # Extract protein NES and pval per dataset
        nes_prot_cols <- grep("^NES_.*_protein$", colnames(plot_data), value = TRUE)
        cancer_types  <- sub("^NES_(.*)_protein$", "\\1", nes_prot_cols)

        for (ct in cancer_types) {
            nes_col  <- paste0("NES_", ct, "_protein")
            pval_col <- paste0("pval_", ct, "_protein")

            if (nes_col %in% colnames(plot_data) && pval_col %in% colnames(plot_data)) {
                tmp <- data.frame(
                    pathway        = plot_data$pathway,
                    collection     = coll_name,
                    dataset        = ct,
                    NES            = plot_data[[nes_col]],
                    pval           = plot_data[[pval_col]],
                    Z_meta_protein = z_meta_prot$Z_meta_protein,
                    stringsAsFactors = FALSE
                )
                all_long_list[[paste(coll_name, ct, sep = "_")]] <- tmp
            }
        }
    }

    if (length(all_long_list) == 0) {
        cat("  [SKIP] No data available for any collection\n\n")
        next
    }

    long_df <- bind_rows(all_long_list)

    # Compute -log10(pval) for bubble size
    long_df$neg_log10_pval <- -log10(pmax(long_df$pval, 1e-300))

    # Remove rows with NA NES
    long_df <- long_df[!is.na(long_df$NES), ]

    if (nrow(long_df) == 0) {
        cat("  [SKIP] No valid data for plotting\n\n")
        next
    }

    # Mark significance (for thick border)
    long_df$is_sig <- long_df$pval < 0.05

    # Format pathway labels (remove collection prefixes for display)
    long_df$pathway_label <- long_df$pathway
    long_df$pathway_label <- gsub("^HALLMARK_", "", long_df$pathway_label)
    long_df$pathway_label <- gsub("^KEGG_", "", long_df$pathway_label)
    # Replace underscores with spaces for readability
    long_df$pathway_label <- gsub("_", " ", long_df$pathway_label)

    # Order pathways by Z_meta_protein (descending: largest at top)
    pathway_order <- long_df %>%
        group_by(pathway_label) %>%
        summarise(Z_meta = first(Z_meta_protein), .groups = "drop") %>%
        arrange(Z_meta) %>%    # ascending so that largest is at the top in ggplot
        pull(pathway_label)
    long_df$pathway_label <- factor(long_df$pathway_label, levels = pathway_order)

    # Order datasets alphabetically (left to right)
    long_df$dataset <- factor(long_df$dataset, levels = sort(unique(long_df$dataset)))

    # Determine symmetric color scale limits
    nes_max <- max(abs(long_df$NES), na.rm = TRUE)
    nes_lim <- ceiling(nes_max * 10) / 10  # round up to 1 decimal

    # Split data by significance for layered plotting
    df_nonsig <- long_df[long_df$is_sig == FALSE, ]
    df_sig    <- long_df[long_df$is_sig == TRUE, ]

    # --- Save data CSV and XLSX ---
    file_prefix <- comp_name

    # Save CSV
    data_csv <- file.path(data_out_dir, paste0(file_prefix, ".csv"))
    fwrite(long_df, data_csv)

    # Save XLSX
    data_xlsx <- file.path(data_out_dir, paste0(file_prefix, ".xlsx"))
    wb_data <- createWorkbook()
    addWorksheet(wb_data, comp_label)
    writeData(wb_data, 1, long_df)

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
             rows = 1, cols = 1:ncol(long_df), gridExpand = TRUE)
    if (nrow(long_df) > 0) {
        addStyle(wb_data, 1, body_style,
                 rows = 2:(nrow(long_df) + 1), cols = 1:ncol(long_df),
                 gridExpand = TRUE, stack = FALSE)
    }
    setColWidths(wb_data, 1, cols = 1:ncol(long_df), widths = "auto")
    saveWorkbook(wb_data, data_xlsx, overwrite = TRUE)

    cat("  Data saved:", basename(data_csv), "&", basename(data_xlsx), "\n")

    # --- Generate bubble plots for each selected color scheme ---
    for (cs_name in selected_color_schemes) {
        cs <- color_schemes[[cs_name]]

        p <- ggplot(long_df, aes(x = dataset, y = pathway_label)) +
            # Layer 1: non-significant bubbles (thin border)
            geom_point(data = df_nonsig,
                       aes(size = neg_log10_pval, fill = NES),
                       shape = 21, color = "black", stroke = 0.3) +
            # Layer 2: significant bubbles (thick black border)
            geom_point(data = df_sig,
                       aes(size = neg_log10_pval, fill = NES),
                       shape = 21, color = "black", stroke = 1.0) +
            scale_x_discrete(limits = sort(unique(long_df$dataset))) +
            scale_fill_gradient2(
                low = cs$low, mid = cs$mid, high = cs$high,
                midpoint = 0,
                limits = c(-nes_lim, nes_lim),
                name = "NES"
            ) +
            scale_size_continuous(
                range = c(1, 8),
                name = expression(-log[10](pval))
            ) +
            labs(
                x = NULL,
                y = NULL,
                title = paste0(gsub("_", " ", comp_name), " (", cs$label, ")")
            ) +
            theme(
                # White background, no gridlines
                plot.background  = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA),
                panel.border     = element_blank(),
                panel.grid       = element_blank(),
                # Axes: draw x and y axis lines
                axis.line.x      = element_line(color = "black", linewidth = 0.4),
                axis.line.y      = element_line(color = "black", linewidth = 0.4),
                axis.ticks       = element_line(color = "black", linewidth = 0.3),
                axis.text.x      = element_text(color = "black", size = 8,
                                                angle = 45, hjust = 1, vjust = 1),
                axis.text.y      = element_text(color = "black", size = 8),
                # Title
                plot.title       = element_text(hjust = 0.5, face = "bold",
                                                size = 10, color = "black"),
                # Legend
                legend.background = element_rect(fill = "white", color = NA),
                legend.key        = element_rect(fill = "white", color = NA),
                legend.text       = element_text(size = 7),
                legend.title      = element_text(size = 8, face = "bold"),
                # Margin
                plot.margin      = margin(10, 15, 10, 10)
            )

        # Calculate dynamic plot dimensions
        n_datasets <- length(unique(long_df$dataset))
        n_pathways <- length(unique(long_df$pathway_label))
        plot_width  <- max(4.5, 2.5 + n_datasets * 0.4 + 1.5)
        plot_height <- max(3.5, 0.8 + n_pathways * 0.3 + 0.8)

        # Save TIFF
        tiff_file <- file.path(output_dir, cs_name,
                               paste0(file_prefix, ".tiff"))
        tiff(tiff_file, width = plot_width, height = plot_height,
             units = "in", res = 300, compression = "lzw")
        print(p)
        dev.off()

        # Save PDF
        pdf_file <- file.path(output_dir, cs_name,
                              paste0(file_prefix, ".pdf"))
        pdf(pdf_file, width = plot_width, height = plot_height)
        print(p)
        dev.off()
    }

    cat("  Bubble plots saved for selected color schemes\n\n")
}

cat("====================================================================\n")
cat("Bubble plots completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

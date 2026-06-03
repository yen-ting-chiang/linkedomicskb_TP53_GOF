################################################################################
# Phosphoprotein without PN GSEA Meta Bubble Plot for Each Dataset
#
# Purpose:
#   Visualize phosphoprotein-level GSEA results across datasets for specific
#   PTMsigDB pathways using bubble plots. Two separate sets of bubble plots
#   are generated:
#     (1) phospho_GSEA_without_PN (without protein normalization)
#     (2) phospho_GSEA_with_PN (with protein normalization)
#   Color encodes NES (enrichment direction/magnitude).
#   Size encodes -log10(p-value) (statistical significance).
#   Bubbles with pval < 0.05 have a thick black border.
#
# Input:
#   - phospho_GSEA_without_PN_meta_vs_ phospho_GSEA_with_PN_meta/
#       {collection}/META_GSEA_{comparison}.csv
#     (from phospho_GSEA_without_PN_meta_vs_ phospho_GSEA_with_PN_meta.R)
#
# Output:
#   - phospho_GSEA_without_PN_meta_vs_phospho_GSEA_with_PN_meta/bubble_plot/
#       phospho_GSEA_without_PN/{color_scheme}/{comparison}.tiff and .pdf
#       phospho_GSEA_with_PN/{color_scheme}/{comparison}.tiff and .pdf
#       data/phospho_GSEA_without_PN/{comparison}.csv and .xlsx
#       data/phospho_GSEA_with_PN/{comparison}.csv and .xlsx
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

# Input directory (integration output; note the space in folder name)
input_dir <- file.path(base_path, "phospho_GSEA_without_PN_meta_vs_phospho_GSEA_with_PN_meta")

# Output directory (without the extra space)
output_dir <- file.path(base_path,
                         "phospho_GSEA_without_PN_meta_vs_phospho_GSEA_with_PN_meta",
                         "bubble_plot")
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

selected_comparisons <- comparisons$name
#selected_comparisons <- c("MUT_GOF_vs_MUT_LOF")


# ==============================================================================
# Section 2: Define Pathways to Plot
# ==============================================================================
# USER: modify this list to add/remove pathways for bubble plots.
# Each entry: collection folder name -> vector of pathway names.
# All pathways across collections will be combined into ONE bubble plot.

pathways_to_plot <- list(
    "PTMsigDB" = c(
        "KINASE-PSP_CDK1",
        "KINASE-PSP_CDK2",
        "KINASE-iKiP_CSNK2A2.CK2A2",
        "KINASE-PSP_CK2A1/CSNK2A1",
        "KINASE-PSP_ATM"
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

# Two analysis levels to plot separately
analysis_levels <- c("phospho_GSEA_without_PN", "phospho_GSEA_with_PN")

cat("====================================================================\n")
cat("Phosphoprotein without PN GSEA Meta Bubble Plot for Each Dataset\n")
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

# Create output subdirectories for each analysis level, color scheme, and data
for (level in analysis_levels) {
    for (cs_name in selected_color_schemes) {
        dir.create(file.path(output_dir, level, cs_name),
                   recursive = TRUE, showWarnings = FALSE)
    }
    dir.create(file.path(output_dir, "data", level),
               recursive = TRUE, showWarnings = FALSE)
}

for (j in seq_len(nrow(comparisons))) {
    comp_name  <- comparisons$name[j]
    comp_label <- comparisons$label[j]

    cat("==================================================================\n")
    cat("Comparison:", comp_name, "\n")
    cat("==================================================================\n")

    # --- Loop over each analysis level (without_PN and with_PN) ---
    for (level in analysis_levels) {

        cat("\n  Level:", level, "\n")

        # --- Load and combine data from all specified collections ---
        all_long_list <- list()

        for (coll_name in names(pathways_to_plot)) {
            target_pathways <- pathways_to_plot[[coll_name]]

            int_file <- file.path(input_dir, coll_name,
                                  paste0("META_GSEA_", comp_name, ".csv"))
            if (!file.exists(int_file)) {
                cat("    [SKIP] File not found for collection:", coll_name, "\n")
                next
            }

            full_data <- fread(int_file)

            # Filter to target pathways
            plot_data <- full_data[pathway %in% target_pathways]
            if (nrow(plot_data) == 0) {
                cat("    [SKIP] No target pathways in:", coll_name, "\n")
                next
            }

            # Retrieve Z_meta for y-axis ordering (use the current level's Z_meta)
            z_meta_col <- paste0("Z_meta_", level)
            if (!(z_meta_col %in% colnames(plot_data))) {
                cat("    [SKIP] Column", z_meta_col, "not found in:", coll_name, "\n")
                next
            }
            z_meta_vals <- plot_data[[z_meta_col]]

            # Extract NES and pval per dataset for the current level
            nes_pattern  <- paste0("^NES_.*_", level, "$")
            nes_cols     <- grep(nes_pattern, colnames(plot_data), value = TRUE)
            cancer_types <- sub(paste0("^NES_(.*)_", level, "$"), "\\1", nes_cols)

            for (ct in cancer_types) {
                nes_col  <- paste0("NES_", ct, "_", level)
                pval_col <- paste0("pval_", ct, "_", level)

                if (nes_col %in% colnames(plot_data) && pval_col %in% colnames(plot_data)) {
                    tmp <- data.frame(
                        pathway    = plot_data$pathway,
                        collection = coll_name,
                        dataset    = ct,
                        NES        = plot_data[[nes_col]],
                        pval       = plot_data[[pval_col]],
                        Z_meta     = z_meta_vals,
                        stringsAsFactors = FALSE
                    )
                    all_long_list[[paste(coll_name, ct, sep = "_")]] <- tmp
                }
            }
        }

        if (length(all_long_list) == 0) {
            cat("    [SKIP] No data available for any collection\n")
            next
        }

        long_df <- bind_rows(all_long_list)

        # Compute -log10(pval) for bubble size
        long_df$neg_log10_pval <- -log10(pmax(long_df$pval, 1e-300))

        # Remove rows with NA NES
        long_df <- long_df[!is.na(long_df$NES), ]

        if (nrow(long_df) == 0) {
            cat("    [SKIP] No valid data for plotting\n")
            next
        }

        # Mark significance (for thick border)
        long_df$is_sig <- long_df$pval < 0.05

        # Format pathway labels for display
        # PTMsigDB pathways: keep as-is (do not strip prefixes)
        long_df$pathway_label <- long_df$pathway

        # Order pathways by Z_meta (descending: largest at top)
        pathway_order <- long_df %>%
            group_by(pathway_label) %>%
            summarise(Z_meta_order = first(Z_meta), .groups = "drop") %>%
            arrange(Z_meta_order) %>%    # ascending so that largest is at the top in ggplot
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
        data_out_dir <- file.path(output_dir, "data", level)

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

        cat("    Data saved:", basename(data_csv), "&", basename(data_xlsx), "\n")

        # --- Generate bubble plots for each selected color scheme ---
        # Pretty level label for plot title
        level_label <- ifelse(level == "phospho_GSEA_without_PN",
                              "Phospho GSEA without PN",
                              "Phospho GSEA with PN")

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
                    title = paste0(gsub("_", " ", comp_name),
                                   " (", level_label, ", ", cs$label, ")")
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
            tiff_file <- file.path(output_dir, level, cs_name,
                                   paste0(file_prefix, ".tiff"))
            tiff(tiff_file, width = plot_width, height = plot_height,
                 units = "in", res = 300, compression = "lzw")
            print(p)
            dev.off()

            # Save PDF
            pdf_file <- file.path(output_dir, level, cs_name,
                                  paste0(file_prefix, ".pdf"))
            pdf(pdf_file, width = plot_width, height = plot_height)
            print(p)
            dev.off()
        }

        cat("    Bubble plots saved for", level, "\n")
    }

    cat("\n")
}

cat("====================================================================\n")
cat("Bubble plots completed!\n")
cat("Output directory:", output_dir, "\n")
cat("  Analysis levels:", paste(analysis_levels, collapse = ", "), "\n")
cat("  Color scheme subdirectories:", paste(selected_color_schemes, collapse = ", "), "\n")
cat("  Data subdirectory: data/\n")
cat("====================================================================\n")

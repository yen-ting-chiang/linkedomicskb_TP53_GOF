################################################################################
# mRNA Meta-GSEA vs Protein Meta-GSEA Integration
#
# Purpose:
#   Integrate mRNA-level and protein-level cross-cancer meta-GSEA results
#   by merging matching comparison files on the "pathway" column for each collection.
#
# Input:
#   - mRNA_differential_analysis_subgroup_adjusted/mRNA_DEG_meta_analysis/mRNA_meta_GSEA/
#       {collection}/META_GSEA_{comparison}.csv
#   - protein_differential_analysis_subgroup_adjusted/protein_DEG_meta_analysis/protein_meta_GSEA/
#       {collection}/META_GSEA_{comparison}.csv
#
# Output:
#   - mRNA_meta_GSEA_vs_protein_meta_GSEA/
#       {collection}/META_GSEA_{comparison}.csv and .xlsx
#       integration_summary.csv and .xlsx
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

# Input directories
mRNA_gsea_dir <- file.path(base_path,
    "mRNA_differential_analysis_subgroup_adjusted", "mRNA_DEG_meta_analysis", "mRNA_meta_GSEA")
protein_gsea_dir <- file.path(base_path,
    "protein_differential_analysis_subgroup_adjusted", "protein_DEG_meta_analysis", "protein_meta_GSEA")

# Output directory
output_dir <- file.path(base_path, "mRNA_meta_GSEA_vs_protein_meta_GSEA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 9 comparisons (must match the naming convention of META_GSEA files)
comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

# Identify collections by looking at subdirectories in mRNA_gsea_dir
collections <- list.dirs(mRNA_gsea_dir, full.names = FALSE, recursive = FALSE)
# Exclude any summary files or non-collection folders (if any)
collections <- collections[collections != ""]

cat("====================================================================\n")
cat("mRNA Meta-GSEA vs Protein Meta-GSEA Integration\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Integration Loop (per collection, per comparison)
# ==============================================================================

integration_summary <- list()

for (coll_name in collections) {
    cat("==================================================================\n")
    cat("Collection:", coll_name, "\n")
    cat("==================================================================\n")
    
    coll_out_dir <- file.path(output_dir, coll_name)
    dir.create(coll_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (j in seq_len(nrow(comparisons))) {
        comp_name <- comparisons$name[j]
        comp_label <- comparisons$label[j]
        
        cat("  Integration for:", comp_name, "\n")
        
        mRNA_file <- file.path(mRNA_gsea_dir, coll_name, paste0("META_GSEA_", comp_name, ".csv"))
        protein_file <- file.path(protein_gsea_dir, coll_name, paste0("META_GSEA_", comp_name, ".csv"))
        
        if (!file.exists(mRNA_file) || !file.exists(protein_file)) {
            cat("    [SKIP] Missing input files\n")
            integration_summary[[paste(coll_name, comp_name, sep = "_")]] <- data.frame(
                collection = coll_name,
                comparison = comp_name,
                n_mRNA_pathways = NA,
                n_protein_pathways = NA,
                n_merged_pathways = NA,
                n_both_sig_005 = NA,
                stringsAsFactors = FALSE
            )
            next
        }
        
        mRNA_data <- fread(mRNA_file)
        protein_data <- fread(protein_file)
        
        # Select core columns, prefix with mRNA_ and protein_
        if ("pathway" %in% colnames(mRNA_data) && "pathway" %in% colnames(protein_data)) {
            mRNA_core <- mRNA_data[, .(
                pathway,
                mRNA_pval = pval,
                mRNA_padj = padj,
                mRNA_log2err = log2err,
                mRNA_ES = ES,
                mRNA_NES = NES,
                mRNA_size = size,
                mRNA_leadingEdge = leadingEdge
            )]
            
            protein_core <- protein_data[, .(
                pathway,
                protein_pval = pval,
                protein_padj = padj,
                protein_log2err = log2err,
                protein_ES = ES,
                protein_NES = NES,
                protein_size = size,
                protein_leadingEdge = leadingEdge
            )]
            
            # Remove duplicate pathways just in case
            mRNA_core <- mRNA_core[!duplicated(pathway)]
            protein_core <- protein_core[!duplicated(pathway)]
            
            # Inner join on pathway
            merged <- merge(mRNA_core, protein_core, by = "pathway")
            
            # Compute a conservative padj threshold (max of both)
            merged$padj_max <- pmax(merged$mRNA_padj, merged$protein_padj, na.rm = TRUE)
            
            # Sort by padj_max
            merged <- merged[order(padj_max)]
            
            n_both_sig <- sum(merged$padj_max < 0.05, na.rm = TRUE)
            
            cat("    Merged pathways:", nrow(merged), "\n")
            cat("    Both significant (padj_max < 0.05):", n_both_sig, "\n")
            
            # Save CSV
            out_csv <- file.path(coll_out_dir, paste0("META_GSEA_", comp_name, ".csv"))
            fwrite(merged, out_csv)
            
            # Save XLSX
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
            
            # Highlight significant rows
            sig_rows <- which(merged$padj_max < 0.05) + 1
            if (length(sig_rows) > 0) {
                sig_style <- createStyle(fgFill = "#FFFFCC")
                addStyle(wb, 1, sig_style,
                         rows = sig_rows, cols = 1:ncol(merged), gridExpand = TRUE,
                         stack = TRUE)
            }
            
            out_xlsx <- file.path(coll_out_dir, paste0("META_GSEA_", comp_name, ".xlsx"))
            saveWorkbook(wb, out_xlsx, overwrite = TRUE)
            
            # Store summary
            integration_summary[[paste(coll_name, comp_name, sep = "_")]] <- data.frame(
                collection = coll_name,
                comparison = comp_name,
                n_mRNA_pathways = nrow(mRNA_core),
                n_protein_pathways = nrow(protein_core),
                n_merged_pathways = nrow(merged),
                n_both_sig_005 = n_both_sig,
                stringsAsFactors = FALSE
            )
        } else {
            cat("    [SKIP] Missing pathway column\n")
        }
    }
    cat("\n")
}

# ==============================================================================
# Section 3: Integration Summary
# ==============================================================================

cat("====================================================================\n")
cat("Generating Integration Summary\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(integration_summary)

if (nrow(summary_df) > 0) {
    fwrite(summary_df, file.path(output_dir, "integration_summary.csv"))
    
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
    
    cat("Integration Complete!\n")
    cat("Results saved to:", output_dir, "\n")
} else {
    cat("No data merged.\n")
}

# ==============================================================================
# Section 4: Scatter Plots for GSEA Integration
# ==============================================================================

cat("\n====================================================================\n")
cat("Generating Scatter Plots\n")
cat("====================================================================\n\n")

plot_dir <- file.path(output_dir, "scatter_plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# User can modify this vector to select which collections to plot
collections_to_plot <- c("Hallmark", "C2_CGP", "C2_CP_KEGG_LEGACY", "C3_TFT_TFT_LEGACY")

for (coll_name in collections_to_plot) {
    cat("Plotting collection:", coll_name, "\n")
    
    for (j in seq_len(nrow(comparisons))) {
        comp_name <- comparisons$name[j]
        
        int_file <- file.path(output_dir, coll_name, paste0("META_GSEA_", comp_name, ".csv"))
        if (!file.exists(int_file)) {
            next
        }
        
        df <- fread(int_file)
        
        # Keep only rows where both NES are finite
        df <- df[is.finite(mRNA_NES) & is.finite(protein_NES)]
        
        if (nrow(df) < 5) {
            next
        }
        
        # Assign colors based on rules
        df$fill_color <- "medium_grey" # Rule 3
        df$fill_color[df$mRNA_padj > 0.05 & df$protein_padj > 0.05] <- "light_grey" # Rule 2
        
        # Identify highlighted pathways (Rule 1)
        # Red: mRNA_NES > 0 & protein_NES > 0 & both padj < 0.05
        path_red <- df$pathway[df$mRNA_NES > 0 & df$protein_NES > 0 & df$mRNA_padj < 0.05 & df$protein_padj < 0.05]
        
        # Blue: mRNA_NES < 0 & protein_NES < 0 & both padj < 0.05
        path_blue <- df$pathway[df$mRNA_NES < 0 & df$protein_NES < 0 & df$mRNA_padj < 0.05 & df$protein_padj < 0.05]
        
        # Orange: protein_NES > 0 & protein_padj < 0.05, excluding (mRNA_NES > 0 & mRNA_padj < 0.05)
        path_orange <- df$pathway[!(df$mRNA_NES > 0 & df$mRNA_padj < 0.05) & df$protein_NES > 0 & df$protein_padj < 0.05]
        
        # Green: protein_NES < 0 & protein_padj < 0.05, excluding (mRNA_NES < 0 & mRNA_padj < 0.05)
        path_green <- df$pathway[!(df$mRNA_NES < 0 & df$mRNA_padj < 0.05) & df$protein_NES < 0 & df$protein_padj < 0.05]
        
        # Apply Rule 1 colors (highest priority, applied last)
        df$fill_color[df$pathway %in% path_green] <- "green"
        df$fill_color[df$pathway %in% path_orange] <- "orange"
        df$fill_color[df$pathway %in% path_blue]  <- "blue"
        df$fill_color[df$pathway %in% path_red]  <- "red"
        
        # Prepare data for label filtering
        df$NES_sum <- df$mRNA_NES + df$protein_NES
        
        # Filter red labels (top 5 NES_sum if > 5)
        if (length(path_red) > 5) {
            red_df <- df[df$pathway %in% path_red, ]
            red_df <- red_df[order(-NES_sum)]
            label_red <- red_df$pathway[1:5]
        } else {
            label_red <- path_red
        }
        
        # Filter blue labels (bottom 5 NES_sum if > 5)
        if (length(path_blue) > 5) {
            blue_df <- df[df$pathway %in% path_blue, ]
            blue_df <- blue_df[order(NES_sum)]
            label_blue <- blue_df$pathway[1:5]
        } else {
            label_blue <- path_blue
        }
        
        # Rule 1 union for labeling
        rule1_label_paths <- unique(c(label_red, label_blue, path_orange, path_green))
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
        
        # Calculate average size for point scaling
        df$avg_size <- (df$mRNA_size + df$protein_size) / 2
        
        # Split data for layered plotting (grey in background, colored in foreground)
        df_grey <- df[fill_color %in% c("light_grey", "medium_grey")]
        df_colored <- df[fill_color %in% c("red", "blue", "orange", "green")]
        df_label <- df[show_label == TRUE]
        
        axis_max <- max(abs(c(df$protein_NES, df$mRNA_NES)), na.rm = TRUE) * 1.05
        
        p <- ggplot(df, aes(x = protein_NES, y = mRNA_NES)) +
            # y=x reference line
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
            # Background grey points
            geom_point(data = df_grey, aes(color = fill_color, size = avg_size), shape = 16, alpha = 0.5) +
            # Foreground colored points
            geom_point(data = df_colored, aes(color = fill_color, size = avg_size), shape = 16, alpha = 0.85) +
            # Pathway labels
            geom_text_repel(
                data = df_label,
                aes(label = label),
                size = 2.2,
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
            scale_size_continuous(range = c(0.5, 3.5), name = "Average Size") +
            coord_cartesian(xlim = c(-axis_max, axis_max), ylim = c(-axis_max, axis_max)) +
            labs(
                x = "Protein NES",
                y = "mRNA NES",
                title = paste(gsub("_", " ", comp_name), "-", coll_name)
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
                plot.margin      = margin(10, 10, 10, 10)
            )
            
        # Save TIFF
        tiff_file <- file.path(plot_dir, paste0("scatter_", coll_name, "_", comp_name, ".tiff"))
        tiff(tiff_file, width = 4.5, height = 3.5, units = "in", res = 300, compression = "lzw")
        print(p)
        dev.off()
        
        # Save PDF
        pdf_file <- file.path(plot_dir, paste0("scatter_", coll_name, "_", comp_name, ".pdf"))
        pdf(pdf_file, width = 4.5, height = 3.5)
        print(p)
        dev.off()
        
        cat("    Saved plots for", comp_name, "\n")
    }
}

cat("====================================================================\n")
cat("Scatter plots completed!\n")
cat("Scatter plots saved to:", plot_dir, "\n")
cat("====================================================================\n")

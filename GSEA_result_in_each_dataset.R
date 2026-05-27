################################################################################
# Integrate mRNA and Protein GSEA Results
#
# Merges GSEA results from mRNA_GSEA_subgroup_adjusted and 
# protein_GSEA_subgroup_adjusted based on the 'pathway' column.
#
# Outputs to GSEA_result_in_each_dataset folder.
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"
mrna_dir <- file.path(base_path, "mRNA_GSEA_subgroup_adjusted")
prot_dir <- file.path(base_path, "protein_GSEA_subgroup_adjusted")
out_dir <- file.path(base_path, "GSEA_result_in_each_dataset")

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

comparisons <- c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
                 "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
                 "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN")

collections <- c(
    "Hallmark", "C2_CGP", "C2_CP_BIOCARTA", "C2_CP_KEGG_LEGACY",
    "C2_CP_KEGG_MEDICUS", "C2_CP_PID", "C2_CP_REACTOME", "C2_CP_WIKIPATHWAYS",
    "C3_MIR_MIRDB", "C3_MIR_MIR_LEGACY", "C3_TFT_GTRD", "C3_TFT_TFT_LEGACY",
    "C4_3CA", "C4_CGN", "C4_CM", "C5_GO_BP", "C5_GO_CC", "C5_GO_MF",
    "C6_Oncogenic", "C9_CellType"
)

cat("====================================================================\n")
cat("Integrating mRNA and Protein GSEA Results\n")
cat("====================================================================\n\n")

for (cancer in cancer_types) {
    cat("==================================================================\n")
    cat(sprintf("Processing dataset: %s\n", cancer))
    cat("==================================================================\n")
    
    for (coll in collections) {
        
        # Check if directories exist for either mRNA or protein before creating output dir
        mrna_coll_dir <- file.path(mrna_dir, cancer, coll)
        prot_coll_dir <- file.path(prot_dir, cancer, coll)
        
        if (!dir.exists(mrna_coll_dir) && !dir.exists(prot_coll_dir)) {
            next
        }
        
        cat(sprintf("  --- Collection: %s ---\n", coll))
        
        # Output directory for this dataset and collection
        ds_coll_out <- file.path(out_dir, cancer, coll)
        dir.create(ds_coll_out, recursive = TRUE, showWarnings = FALSE)
        
        for (comp in comparisons) {
            file_name <- paste0("GSEA_", comp, ".csv")
            mrna_file <- file.path(mrna_coll_dir, file_name)
            prot_file <- file.path(prot_coll_dir, file_name)
            
            # Check if either file exists
            if (!file.exists(mrna_file) && !file.exists(prot_file)) {
                next
            }
            
            mrna_data <- data.frame(pathway = character())
            prot_data <- data.frame(pathway = character())
            
            if (file.exists(mrna_file)) {
                tmp_mrna <- fread(mrna_file)
                if (nrow(tmp_mrna) > 0) {
                    # Rename columns to add _mRNA suffix except pathway
                    cols_to_rename <- colnames(tmp_mrna)[colnames(tmp_mrna) != "pathway"]
                    colnames(tmp_mrna)[colnames(tmp_mrna) != "pathway"] <- paste0(cols_to_rename, "_mRNA")
                    mrna_data <- tmp_mrna
                }
            }
            
            if (file.exists(prot_file)) {
                tmp_prot <- fread(prot_file)
                if (nrow(tmp_prot) > 0) {
                    # Rename columns to add _protein suffix except pathway
                    cols_to_rename <- colnames(tmp_prot)[colnames(tmp_prot) != "pathway"]
                    colnames(tmp_prot)[colnames(tmp_prot) != "pathway"] <- paste0(cols_to_rename, "_protein")
                    prot_data <- tmp_prot
                }
            }
            
            if (nrow(mrna_data) == 0 && nrow(prot_data) == 0) {
                next
            }
            
            # Merge by pathway
            merged_data <- full_join(mrna_data, prot_data, by = "pathway")
            
            # Save CSV
            out_csv <- file.path(ds_coll_out, file_name)
            fwrite(merged_data, out_csv)
            
            # Save XLSX
            out_xlsx <- file.path(ds_coll_out, paste0("GSEA_", comp, ".xlsx"))
            wb <- createWorkbook()
            # Use label from comparisons for sheet name if needed, or just comparison name
            # Sheet names in Excel have a 31 character limit
            sheet_name <- comp
            if (nchar(sheet_name) > 31) {
                sheet_name <- substr(sheet_name, 1, 31)
            }
            addWorksheet(wb, sheet_name)
            writeData(wb, 1, merged_data)
            saveWorkbook(wb, out_xlsx, overwrite = TRUE)
            
            cat(sprintf("    Merged: %s\n", comp))
        }
    }
    cat("\n")
}

cat("====================================================================\n")
cat("Integration Complete!\n")
cat(sprintf("Results saved to: %s\n", out_dir))
cat("====================================================================\n\n")

# ==============================================================================
# Section 6: Bubble Plot for Specific Pathways Across Datasets
# ==============================================================================

cat("====================================================================\n")
cat("Generating Cross-Dataset Bubble Plots\n")
cat("====================================================================\n\n")

# User Configuration for Bubble Plot
# ------------------------------------------------------------------------------
target_collection <- "Hallmark"
# You can change the target pathways here:
target_pathways <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

# Color palette options for NES Positive / Negative
# Option 1: Nature/Lancet style (Red/Blue) -> c("#E64B35", "#4DBBD5")
# Option 2: Science style (Orange/Purple)  -> c("#F39B7F", "#8491B4")
# Option 3: JAMA style (Red/Green)         -> c("#DC0000", "#00A087")
# Option 4: Classic (Red/Blue)             -> c("red3", "royalblue")

# Select your preferred color palette here:
selected_palette <- c("Positive" = "#E64B35", "Negative" = "#4DBBD5") # Default is Nature style
# ------------------------------------------------------------------------------

bubble_out_dir <- file.path(out_dir, "bubble")
dir.create(bubble_out_dir, recursive = TRUE, showWarnings = FALSE)

# Iterate through each specified pathway
for (pathway_name in target_pathways) {
    cat(sprintf("Creating bubble plots for pathway: %s\n", pathway_name))
    
    # Iterate through each comparison (e.g., TP53mt_vs_TP53wt)
    for (comp in comparisons) {
        
        # We will collect data for this specific pathway and comparison across all datasets
        plot_data_list <- list()
        
        for (cancer in cancer_types) {
            csv_path <- file.path(out_dir, cancer, target_collection, paste0("GSEA_", comp, ".csv"))
            
            if (file.exists(csv_path)) {
                df <- fread(csv_path)
                # Check if the target pathway exists in this dataframe
                df_pathway <- df %>% filter(pathway == pathway_name)
                
                if (nrow(df_pathway) > 0) {
                    # Extract mRNA data
                    if ("NES_mRNA" %in% colnames(df_pathway) && "padj_mRNA" %in% colnames(df_pathway)) {
                        plot_data_list[[paste0(cancer, "_mRNA")]] <- data.frame(
                            Dataset = cancer,
                            Omics = "mRNA",
                            NES = as.numeric(df_pathway$NES_mRNA),
                            padj = as.numeric(df_pathway$padj_mRNA),
                            stringsAsFactors = FALSE
                        )
                    }
                    
                    # Extract Protein data
                    if ("NES_protein" %in% colnames(df_pathway) && "padj_protein" %in% colnames(df_pathway)) {
                        plot_data_list[[paste0(cancer, "_protein")]] <- data.frame(
                            Dataset = cancer,
                            Omics = "Protein",
                            NES = as.numeric(df_pathway$NES_protein),
                            padj = as.numeric(df_pathway$padj_protein),
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
        
        if (length(plot_data_list) == 0) {
            next
        }
        
        # Combine data
        plot_df <- bind_rows(plot_data_list)
        
        # Remove NA values to avoid plotting errors
        plot_df <- plot_df[!is.na(plot_df$NES) & !is.na(plot_df$padj), ]
        
        if (nrow(plot_df) == 0) {
            next
        }
        
        # Data preparation for plotting
        # 1. Omics factor: mRNA on top, Protein on bottom
        plot_df$Omics <- factor(plot_df$Omics, levels = c("Protein", "mRNA"))
        
        # 2. Dataset factor to keep order
        plot_df$Dataset <- factor(plot_df$Dataset, levels = cancer_types)
        
        # 3. NES direction and Absolute NES for size
        plot_df$Direction <- ifelse(plot_df$NES > 0, "Positive", "Negative")
        plot_df$Abs_NES <- abs(plot_df$NES)
        
        # 4. -log10(padj) for color intensity
        # Replace 0 padj with a very small number to avoid infinite -log10
        min_padj <- min(plot_df$padj[plot_df$padj > 0], na.rm = TRUE)
        if (is.infinite(min_padj) || is.na(min_padj)) min_padj <- 1e-10
        plot_df$padj_safe <- ifelse(plot_df$padj == 0, min_padj / 10, plot_df$padj)
        plot_df$NegLogPadj <- -log10(plot_df$padj_safe)
        
        # Clean pathway name for display (remove "HALLMARK_")
        display_pathway <- sub("^HALLMARK_", "", pathway_name)
        
        # Save plotting data
        safe_pathway_name <- gsub("[^A-Za-z0-9_]", "_", display_pathway)
        data_prefix <- file.path(bubble_out_dir, paste0("BubbleData_", safe_pathway_name, "_", comp))
        fwrite(plot_df, paste0(data_prefix, ".csv"))
        
        wb_data <- createWorkbook()
        addWorksheet(wb_data, "Data")
        writeData(wb_data, "Data", plot_df)
        saveWorkbook(wb_data, paste0(data_prefix, ".xlsx"), overwrite = TRUE)
        
        # Create Plot
        p <- ggplot(plot_df, aes(x = Dataset, y = Omics)) +
            geom_point(aes(size = Abs_NES, color = Direction, alpha = NegLogPadj)) +
            scale_color_manual(values = selected_palette) +
            scale_alpha_continuous(range = c(0.2, 1), name = "-log10(padj)") +
            scale_size_continuous(range = c(3, 10), name = "|NES|") +
            labs(
                title = display_pathway,
                x = NULL,
                y = NULL
            ) +
            theme_minimal(base_size = 14) +
            theme(
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black", linewidth = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
                axis.text.y = element_text(color = "black", face = "bold"),
                plot.title = element_text(face = "bold", hjust = 0.5, color = "black")
            )
            
        # Save plots
        plot_prefix <- file.path(bubble_out_dir, paste0("BubblePlot_", safe_pathway_name, "_", comp))
        
        suppressWarnings({
            ggsave(paste0(plot_prefix, ".pdf"), plot = p, width = 8, height = 4, bg = "white")
            ggsave(paste0(plot_prefix, ".tiff"), plot = p, width = 8, height = 4, dpi = 300, compression = "lzw", bg = "white")
        })
        
        cat(sprintf("    Generated plots for: %s\n", comp))
    }
}

cat("====================================================================\n")
cat("Bubble Plots Generation Complete!\n")
cat(sprintf("Results saved to: %s\n", bubble_out_dir))
cat("====================================================================\n")

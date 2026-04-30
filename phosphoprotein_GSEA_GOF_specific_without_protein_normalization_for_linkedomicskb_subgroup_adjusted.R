################################################################################
# TP53 Phosphoprotein GSEA - GOF Specific Pathway Extraction (LinkedOmicsKB)
#
# This script filters the PTMsigDB GSEA summary results (without protein normalization)
# to identify pathways associated with GOF (Gain of Function) activity.
#
# Filtering Logic:
# 1. "GOF_vs_LOF" or "Hot_vs_LOF" must appear in "Pos_Significant_In" or "Neg_Significant_In".
# 2. Reversal check:
#    - If "GOF_vs_LOF" or "Hot_vs_LOF" is in "Pos_Significant_In", then "LOF_vs_wt" 
#      must NOT be in "Neg_Significant_In".
#    - If "GOF_vs_LOF" or "Hot_vs_LOF" is in "Neg_Significant_In", then "LOF_vs_wt" 
#      must NOT be in "Pos_Significant_In".
#
# Outputs: 
#   - Per-dataset: GOF_Specific_Pathway_Summary_[Dataset]_PTMsigDB.csv/xlsx
#   - Cross-dataset: GOF_Specific_Cross_Dataset_Summary_PTMsigDB.csv/xlsx
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"
gsea_dir <- file.path(base_path, "phosphoprotein_GSEA_without_protein_normalization_subgroup_adjusted")

datasets <- data.frame(
    folder = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
               "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    cancer_type = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
                    "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    stringsAsFactors = FALSE
)

collections <- c("PTMsigDB")

# Helper function to check if string exists in comma-separated list
contains_label <- function(text, label) {
    # Use regex to find exact match in comma-separated list
    grepl(paste0("(^|, )", label, "(, |$)"), text)
}

cat("====================================================================\n")
cat("TP53 Phosphoprotein GSEA - GOF Specific Extraction (without protein normalization)\n")
cat("====================================================================\n\n")

if (!dir.exists(gsea_dir)) {
    stop("GSEA directory not found. Please run the GSEA summary script first.")
}

# List to collect results for cross-dataset summary
all_gof_results <- list()

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]

    cat("Processing GOF-specific summary for Dataset:", ds_folder, "(", ds_cancer, ")\n")

    ds_out <- file.path(gsea_dir, ds_folder)
    if (!dir.exists(ds_out)) {
        cat("  [SKIP] Dataset GSEA directory not found\n")
        next
    }

    for (coll in collections) {
        coll_dir <- file.path(ds_out, coll)
        if (!dir.exists(coll_dir)) {
            next
        }

        # Read the existing detailed summary
        summary_file <- file.path(coll_dir, paste0("Detailed_Pathway_Summary_", ds_folder, "_", coll, ".csv"))
        
        if (!file.exists(summary_file)) {
            cat("  [SKIP] Detailed summary file not found:", basename(summary_file), "\n")
            next
        }

        summary_df <- fread(summary_file)
        
        if (nrow(summary_df) == 0) {
            next
        }

        # Apply GOF-specific filtering logic
        gof_specific_df <- summary_df %>%
            filter(
                ( (contains_label(Pos_Significant_In, "GOF_vs_LOF") | contains_label(Pos_Significant_In, "Hot_vs_LOF")) & 
                  !contains_label(Neg_Significant_In, "LOF_vs_wt") ) |
                ( (contains_label(Neg_Significant_In, "GOF_vs_LOF") | contains_label(Neg_Significant_In, "Hot_vs_LOF")) & 
                  !contains_label(Pos_Significant_In, "LOF_vs_wt") )
            )

        if (nrow(gof_specific_df) > 0) {
            cat("  Found", nrow(gof_specific_df), "GOF-specific pathways\n")
            
            # Add dataset info for later cross-dataset summary
            gof_specific_df <- gof_specific_df %>%
                mutate(dataset = ds_folder)
            
            all_gof_results[[paste0(ds_folder, "_", coll)]] <- gof_specific_df

            # Save per-dataset outputs
            out_prefix <- file.path(coll_dir, paste0("GOF_Specific_Pathway_Summary_", ds_folder, "_", coll))
            fwrite(gof_specific_df, paste0(out_prefix, ".csv"))
            
            wb <- createWorkbook()
            addWorksheet(wb, "GOF_Specific_Pathways")
            writeData(wb, "GOF_Specific_Pathways", gof_specific_df)
            header_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom", fgFill = "#C6E0B4")
            addStyle(wb, "GOF_Specific_Pathways", header_style, rows = 1, cols = 1:ncol(gof_specific_df), gridExpand = TRUE)
            setColWidths(wb, "GOF_Specific_Pathways", cols = 1:ncol(gof_specific_df), widths = "auto")
            saveWorkbook(wb, paste0(out_prefix, ".xlsx"), overwrite = TRUE)
        }
    }
    cat("\n")
}

# ==============================================================================
# Cross-Dataset Summary Logic
# ==============================================================================
cat("Generating Cross-Dataset Summary...\n")

if (length(all_gof_results) > 0) {
    combined_gof <- bind_rows(all_gof_results)
    
    # Tag directions for each entry
    combined_gof <- combined_gof %>%
        mutate(
            Is_Pos = (contains_label(Pos_Significant_In, "GOF_vs_LOF") | contains_label(Pos_Significant_In, "Hot_vs_LOF")),
            Is_Neg = (contains_label(Neg_Significant_In, "GOF_vs_LOF") | contains_label(Neg_Significant_In, "Hot_vs_LOF"))
        )
    
    # Aggregate by pathway
    cross_summary <- combined_gof %>%
        group_by(pathway) %>%
        summarize(
            Frequency_Pos = sum(Is_Pos),
            Datasets_Pos = paste(dataset[Is_Pos], collapse = ", "),
            Frequency_Neg = sum(Is_Neg),
            Datasets_Neg = paste(dataset[Is_Neg], collapse = ", "),
            Total_Datasets = n_distinct(dataset)
        ) %>%
        arrange(desc(Total_Datasets), desc(Frequency_Pos + Frequency_Neg))

    # Save Cross-Dataset outputs
    cross_csv <- file.path(gsea_dir, "GOF_Specific_Cross_Dataset_Summary_PTMsigDB.csv")
    fwrite(cross_summary, cross_csv)
    
    cross_xlsx <- file.path(gsea_dir, "GOF_Specific_Cross_Dataset_Summary_PTMsigDB.xlsx")
    wb_cross <- createWorkbook()
    addWorksheet(wb_cross, "GOF_Cross_Summary")
    writeData(wb_cross, "GOF_Cross_Summary", cross_summary)
    header_style_cross <- createStyle(textDecoration = "bold", halign = "center", border = "bottom", fgFill = "#A9D08E")
    addStyle(wb_cross, "GOF_Cross_Summary", header_style_cross, rows = 1, cols = 1:ncol(cross_summary), gridExpand = TRUE)
    setColWidths(wb_cross, "GOF_Cross_Summary", cols = 1:ncol(cross_summary), widths = "auto")
    saveWorkbook(wb_cross, cross_xlsx, overwrite = TRUE)
    
    cat("  Cross-dataset summary saved to:", basename(cross_xlsx), "\n")
} else {
    cat("  No GOF-specific pathways found across any dataset.\n")
}

cat("\n====================================================================\n")
cat("GOF-specific Extraction & Summary Complete!\n")
cat("====================================================================\n")

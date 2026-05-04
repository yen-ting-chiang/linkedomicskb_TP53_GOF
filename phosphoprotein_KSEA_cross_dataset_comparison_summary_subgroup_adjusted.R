################################################################################
# TP53 Phosphoprotein KSEA - Cross-Dataset Comparison Summary (LinkedOmicsKB)
#
# This script aggregates the KSEA results (PhosphoSitePlus) across all 10 datasets 
# for each of the 9 TP53 mutation comparisons.
#
# For each comparison and each kinase, it calculates:
#   - Frequency_Pos: Number of datasets where FDR < 0.05 and z.score > 0
#   - Datasets_Pos: List of datasets where positive significance was found
#   - Frequency_Neg: Number of datasets where FDR < 0.05 and z.score < 0
#   - Datasets_Neg: List of datasets where negative significance was found
#
# Outputs: Cross_Dataset_Comparison_Summary_Kinases.csv/xlsx in phosphoprotein_KSEA_subgroup_adjusted/
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"
ksea_dir <- file.path(base_path, "phosphoprotein_KSEA_subgroup_adjusted")

datasets <- data.frame(
    folder = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
               "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    cancer_type = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
                    "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    stringsAsFactors = FALSE
)

comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

cat("====================================================================\n")
cat("TP53 Phosphoprotein KSEA - Cross-Dataset Comparison Summary\n")
cat("====================================================================\n\n")

if (!dir.exists(ksea_dir)) {
    stop("KSEA directory not found. Please run the KSEA script first.")
}

comp_list <- list()

for (j in seq_len(nrow(comparisons))) {
    comp_name <- comparisons$name[j]
    comp_label <- comparisons$label[j]
    
    cat("Summarizing comparison:", comp_label, "...\n")
    
    all_ds_results <- list()

    for (i in seq_len(nrow(datasets))) {
        ds_folder <- datasets$folder[i]
        
        # KSEA files are directly in the dataset folder
        ksea_file <- file.path(ksea_dir, ds_folder, paste0("KSEA_", comp_name, ".csv"))
        
        if (file.exists(ksea_file)) {
            res <- tryCatch(fread(ksea_file), error = function(e) NULL)
            
            if (!is.null(res) && nrow(res) > 0 && "Kinase.Gene" %in% colnames(res)) {
                # Filter for significance and add dataset tag
                sig_res <- res %>%
                    filter(FDR < 0.05) %>%
                    mutate(dataset = ds_folder) %>%
                    select(Kinase.Gene, z.score, dataset)
                
                if (nrow(sig_res) > 0) {
                    all_ds_results[[ds_folder]] <- sig_res
                }
            }
        }
    }

    if (length(all_ds_results) > 0) {
        combined_ds <- bind_rows(all_ds_results)
        
        # Aggregate by kinase for THIS comparison
        comp_summary <- combined_ds %>%
            group_by(Kinase.Gene) %>%
            summarize(
                Comparison = comp_label,
                Frequency_Pos = sum(z.score > 0),
                Datasets_Pos = paste(dataset[z.score > 0], collapse = ", "),
                Frequency_Neg = sum(z.score < 0),
                Datasets_Neg = paste(dataset[z.score < 0], collapse = ", "),
                Total_Datasets = n_distinct(dataset)
            ) %>%
            mutate(Total_Frequency = Frequency_Pos + Frequency_Neg) %>%
            arrange(desc(Total_Frequency))
        
        comp_list[[comp_label]] <- comp_summary
    }
}

if (length(comp_list) > 0) {
    master_summary <- bind_rows(comp_list) %>%
        arrange(Comparison, desc(Total_Frequency))

    # Save outputs
    out_csv <- file.path(ksea_dir, "Cross_Dataset_Comparison_Summary_Kinases.csv")
    fwrite(master_summary, out_csv)
    
    out_xlsx <- file.path(ksea_dir, "Cross_Dataset_Comparison_Summary_Kinases.xlsx")
    wb <- createWorkbook()
    
    header_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom", fgFill = "#BDD7EE")

    # Add a sheet for each comparison
    for (label in names(comp_list)) {
        sheet_name <- substr(label, 1, 31) 
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, comp_list[[label]])
        addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(comp_list[[label]]), gridExpand = TRUE)
        setColWidths(wb, sheet_name, cols = 1:ncol(comp_list[[label]]), widths = "auto")
    }
    
    # Also add a Master sheet
    addWorksheet(wb, "Master_Summary")
    writeData(wb, "Master_Summary", master_summary)
    addStyle(wb, "Master_Summary", header_style, rows = 1, cols = 1:ncol(master_summary), gridExpand = TRUE)
    setColWidths(wb, "Master_Summary", cols = 1:ncol(master_summary), widths = "auto")
    
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    
    cat("\nCross-dataset comparison summary complete!\n")
    cat("  CSV:", basename(out_csv), "\n")
    cat("  XLSX:", basename(out_xlsx), "\n")
} else {
    cat("\nNo significant KSEA results found across any comparisons/datasets.\n")
}

cat("\n====================================================================\n")
cat("Aggregation Complete!\n")
cat("====================================================================\n")

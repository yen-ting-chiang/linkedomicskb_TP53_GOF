################################################################################
# TP53 Phosphoprotein Gene-Level Preranked GSEA (with protein normalization)
#
# Input: Differential phosphosite (DPS) results from
#   phosphoprotein_differential_analysis_subgroup_adjusted/
#
# Two gene-level aggregation methods:
#   1. Max Abs: gene score = t-statistic of the phosphosite with the largest
#      |t| per gene (direction preserved)
#      Reference: Krug et al. Cell 2020 (CPTAC); Mertins et al. Nature 2016
#   2. Median: gene score = median t-statistic across all phosphosites per gene
#      Reference: Hernandez-Armenta et al. Mol Syst Biol 2017;
#                 Ochoa et al. PNAS 2016
#
# Gene sets: Hallmark, C2 (CGP, BIOCARTA, KEGG_LEGACY, KEGG_MEDICUS, PID,
#            REACTOME, WIKIPATHWAYS), C3 (MIR:MIRDB, MIR:MIR_LEGACY,
#            TFT:GTRD, TFT:TFT_LEGACY), C4 (3CA, CGN, CM),
#            C5 (GO:BP, GO:CC, GO:MF), C6 (Oncogenic)
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspots vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspots vs TP53wt,
#   DN vs TP53wt, Non-DN vs TP53wt, DN vs non-DN
#
# Methodology: fgsea preranked GSEA with msigdbr gene sets
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(msigdbr)
  library(fgsea)
  library(openxlsx)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")
datasets <- data.frame(
  folder = cancer_types,
  cancer_type = cancer_types,
  stringsAsFactors = FALSE
)

# Input: phosphosite DPS results (with protein normalization)
dps_dir <- file.path(base_path, "phosphoprotein_differential_analysis_subgroup_adjusted")

# Output: gene-level GSEA results
output_dir <- file.path(base_path, "phosphoprotein_with_PN_gene_level_GSEA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 9 comparisons with their file prefixes
comparisons <- data.frame(
  name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
           "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
           "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
  label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
            "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
            "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
  stringsAsFactors = FALSE
)

# Aggregation method selection
# Options: "both" (default), "max_abs_t", or "median_t"
methods_to_run <- "both"
aggregation_methods <- if (methods_to_run == "both") {
  c("max_abs_t", "median_t")
} else {
  methods_to_run
}

# fgsea parameters
minSize <- 15
maxSize <- 500

cat("====================================================================\n")
cat("TP53 Phosphoprotein Gene-Level Preranked GSEA\n")
cat("(with protein normalization)\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Load Gene Sets from MSigDB
# ==============================================================================

cat("Loading MSigDB gene sets...\n")

# Build a global mapping from ensembl_gene to gene_symbol
cat("Building Ensembl to Gene Symbol mapping...\n")
all_genes_df <- suppressMessages(msigdbr(species = "Homo sapiens"))
ens2sym_df <- all_genes_df[!is.na(all_genes_df$ensembl_gene) & all_genes_df$ensembl_gene != "", ]
ens2sym_df <- ens2sym_df[!duplicated(ens2sym_df$ensembl_gene), ]
ens2sym_dict <- setNames(ens2sym_df$gene_symbol, ens2sym_df$ensembl_gene)
rm(all_genes_df, ens2sym_df) # free memory

# Helper: convert msigdbr output to named list for fgsea
msigdbr_to_list <- function(df) {
  df <- df[!is.na(df$ensembl_gene) & df$ensembl_gene != "", ]
  split(df$ensembl_gene, df$gs_name)
}

collections_to_fetch <- list(
  list(name="Hallmark",               cat="H",  sub=NULL),
  list(name="C2_CGP",                 cat="C2", sub="CGP"),
  list(name="C2_CP_BIOCARTA",         cat="C2", sub="CP:BIOCARTA"),
  list(name="C2_CP_KEGG_LEGACY",      cat="C2", sub="CP:KEGG_LEGACY"),
  list(name="C2_CP_KEGG_MEDICUS",     cat="C2", sub="CP:KEGG_MEDICUS"),
  list(name="C2_CP_PID",              cat="C2", sub="CP:PID"),
  list(name="C2_CP_REACTOME",         cat="C2", sub="CP:REACTOME"),
  list(name="C2_CP_WIKIPATHWAYS",     cat="C2", sub="CP:WIKIPATHWAYS"),
  list(name="C3_MIR_MIRDB",           cat="C3", sub="MIR:MIRDB"),
  list(name="C3_MIR_MIR_LEGACY",      cat="C3", sub="MIR:MIR_LEGACY"),
  list(name="C3_TFT_GTRD",            cat="C3", sub="TFT:GTRD"),
  list(name="C3_TFT_TFT_LEGACY",      cat="C3", sub="TFT:TFT_LEGACY"),
  list(name="C4_3CA",                 cat="C4", sub="3CA"),
  list(name="C4_CGN",                 cat="C4", sub="CGN"),
  list(name="C4_CM",                  cat="C4", sub="CM"),
  list(name="C5_GO_BP",               cat="C5", sub="GO:BP"),
  list(name="C5_GO_CC",               cat="C5", sub="GO:CC"),
  list(name="C5_GO_MF",               cat="C5", sub="GO:MF"),
  list(name="C6_Oncogenic",           cat="C6", sub=NULL)
)

collections <- list()

for (info in collections_to_fetch) {
  if (is.null(info$sub)) {
    df <- suppressMessages(msigdbr(species = "Homo sapiens", collection = info$cat))
  } else {
    df <- suppressMessages(msigdbr(species = "Homo sapiens", collection = info$cat, subcollection = info$sub))
  }
  df <- df %>% dplyr::select(gs_name, ensembl_gene)
  collections[[info$name]] <- msigdbr_to_list(df)
  cat(sprintf("  %-22s : %d gene sets\n", info$name, length(collections[[info$name]])))
}
cat("\n")

# ==============================================================================
# Section 3: Gene-Level Aggregation Functions
# ==============================================================================

#' Aggregate phosphosite-level DPS results to gene-level scores
#'
#' @param dps_file Path to DPS CSV file (from differential phosphosite analysis)
#' @param method Aggregation method: "max_abs_t" or "median_t"
#' @return Named numeric vector of gene-level scores (Ensembl gene IDs, version stripped)
#'         sorted in decreasing order, or NULL if insufficient data
aggregate_phospho_to_gene <- function(dps_file, method = "max_abs_t") {
  if (!file.exists(dps_file)) return(NULL)
  
  dps <- fread(dps_file)
  if (!("phosphosite" %in% colnames(dps)) || !("t" %in% colnames(dps))) {
    cat("    [WARN] Missing phosphosite/t columns in:", basename(dps_file), "\n")
    return(NULL)
  }
  
  # Extract Ensembl gene ID from phosphosite column
  # Format: ENSG00000183765.22|ENSP00000385747.1|S260|KRKFAIGSAREADPA|1
  # Gene ID is the first field before "|"
  dps$ensembl_gene <- sub("\\|.*$", "", dps$phosphosite)
  
  # Strip Ensembl version suffix (e.g., ENSG00000183765.22 -> ENSG00000183765)
  dps$ensembl_gene <- sub("\\.[0-9]+$", "", dps$ensembl_gene)
  
  # Remove rows with missing t-statistics or gene IDs
  dps <- dps[is.finite(dps$t) & !is.na(dps$ensembl_gene) & dps$ensembl_gene != "", ]
  
  if (nrow(dps) == 0) {
    cat("    [WARN] No valid phosphosites after filtering\n")
    return(NULL)
  }
  
  # Aggregate to gene level
  if (method == "max_abs_t") {
    # For each gene, select the phosphosite with the largest |t|,
    # then use its original (signed) t-statistic as the gene score
    gene_scores <- dps %>%
      group_by(ensembl_gene) %>%
      slice_max(abs(t), n = 1, with_ties = FALSE) %>%
      ungroup()
    stats_vec <- setNames(gene_scores$t, gene_scores$ensembl_gene)
  } else if (method == "median_t") {
    # For each gene, compute the median t-statistic across all phosphosites
    gene_scores <- dps %>%
      group_by(ensembl_gene) %>%
      summarise(t_median = median(t, na.rm = TRUE), .groups = "drop")
    stats_vec <- setNames(gene_scores$t_median, gene_scores$ensembl_gene)
  } else {
    stop("Unknown aggregation method: ", method)
  }
  
  # Remove duplicates (should not exist after aggregation, but safety check)
  stats_vec <- stats_vec[!duplicated(names(stats_vec))]
  
  # Remove non-finite values
  stats_vec <- stats_vec[is.finite(stats_vec)]
  
  if (length(stats_vec) < 100) {
    cat("    [WARN] Too few genes after aggregation:", length(stats_vec), "\n")
    return(NULL)
  }
  
  # Sort descending (required by fgsea)
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  stats_vec
}

# ==============================================================================
# Section 4: Run Preranked GSEA
# ==============================================================================

#' Run fgsea preranked GSEA from a gene-level score vector
#' @param stats_vec Named numeric vector of gene scores (sorted descending)
#' @param pathways Named list of gene sets
#' @param mapping_dict Named vector mapping Ensembl IDs to Gene Symbols
#' @return fgsea result data.frame or NULL
run_preranked_gsea <- function(stats_vec, pathways, mapping_dict = ens2sym_dict) {
  if (is.null(stats_vec) || length(stats_vec) < 100) return(NULL)
  
  # Run fgsea
  res <- tryCatch({
    suppressWarnings(fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats_vec,
      minSize = minSize,
      maxSize = maxSize,
      eps = 0
    ))
  }, error = function(e) {
    cat("    [ERROR] fgsea failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # Convert leadingEdge list to semicolon-separated string for CSV/XLSX output
  # And map Ensembl IDs to Gene Symbols
  res$leadingEdge <- vapply(res$leadingEdge, function(x) {
    symbols <- mapping_dict[x]
    symbols[is.na(symbols) | symbols == ""] <- x[is.na(symbols) | symbols == ""]
    paste(symbols, collapse = ";")
  }, character(1))
  
  # Sort by NES descending
  res <- res[order(res$NES, decreasing = TRUE), ]
  as.data.frame(res)
}

# ==============================================================================
# Section 5: Main Processing Loop
# ==============================================================================

for (agg_method in aggregation_methods) {
  
  cat("####################################################################\n")
  cat("Aggregation Method:", agg_method, "\n")
  if (agg_method == "max_abs_t") {
    cat("  (gene score = t-statistic of phosphosite with largest |t| per gene)\n")
  } else {
    cat("  (gene score = median t-statistic across all phosphosites per gene)\n")
  }
  cat("####################################################################\n\n")
  
  # Method-specific output subdirectory
  method_output_dir <- file.path(output_dir, agg_method)
  dir.create(method_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Storage for summary statistics (one per collection)
  all_summaries <- list()
  for (coll_name in names(collections)) {
    all_summaries[[coll_name]] <- list()
  }
  
  for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]
    
    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ") |", agg_method, "\n")
    cat("==================================================================\n")
    
    ds_dps_dir <- file.path(dps_dir, ds_folder)
    if (!dir.exists(ds_dps_dir)) {
      cat("  [SKIP] DPS directory not found\n\n")
      next
    }
    
    ds_out <- file.path(method_output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)
    
    for (coll_name in names(collections)) {
      pathways <- collections[[coll_name]]
      cat("\n  --- Collection:", coll_name, "---\n")
      
      # Create collection output subdirectory
      coll_out <- file.path(ds_out, coll_name)
      dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)
      
      # Summary row for this dataset + collection
      summary_row <- list(dataset = ds_folder, cancer_type = ds_cancer)
      
      for (j in seq_len(nrow(comparisons))) {
        comp_name <- comparisons$name[j]
        comp_label <- comparisons$label[j]
        
        dps_file <- file.path(ds_dps_dir, paste0("DPS_", comp_name, ".csv"))
        cat("    ", comp_label, ":", sep = "")
        
        # Aggregate phosphosites to gene-level scores
        stats_vec <- aggregate_phospho_to_gene(dps_file, method = agg_method)
        
        if (!is.null(stats_vec)) {
          cat(" ", length(stats_vec), "genes |", sep = "")
        }
        
        # Run GSEA
        res <- run_preranked_gsea(stats_vec, pathways)
        
        if (!is.null(res)) {
          n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
          n_up <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
          n_down <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
          cat(" ", nrow(res), "sets tested,", n_sig, "significant\n")
          
          # Save CSV (also sorted by NES descending)
          fwrite(res, file.path(coll_out, paste0("GSEA_", comp_name, ".csv")))
          
          # Save XLSX with padj < 0.05 highlighting
          wb <- createWorkbook()
          addWorksheet(wb, comp_label)
          writeData(wb, 1, res)
          
          # Header style
          hs <- createStyle(
            textDecoration = "bold", halign = "center",
            border = "bottom", fgFill = "#4472C4", fontColour = "white"
          )
          addStyle(wb, 1, hs, rows = 1, cols = 1:ncol(res), gridExpand = TRUE)
          
          # Highlight rows where padj < 0.05
          padj_col <- which(colnames(res) == "padj")
          if (length(padj_col) == 1) {
            sig_rows <- which(res$padj < 0.05) + 1  # +1 for header row
            if (length(sig_rows) > 0) {
              sig_style <- createStyle(fgFill = "#C6EFCE")
              addStyle(wb, 1, sig_style,
                       rows = sig_rows, cols = 1:ncol(res),
                       gridExpand = TRUE, stack = TRUE)
            }
          }
          
          setColWidths(wb, 1, cols = 1:ncol(res), widths = "auto")
          saveWorkbook(wb, file.path(coll_out, paste0("GSEA_", comp_name, ".xlsx")),
                       overwrite = TRUE)
          
          summary_row[[paste0(comp_name, "_total")]] <- nrow(res)
          summary_row[[paste0(comp_name, "_sig")]] <- n_sig
          summary_row[[paste0(comp_name, "_up")]] <- n_up
          summary_row[[paste0(comp_name, "_down")]] <- n_down
        } else {
          cat(" no DPS file or insufficient data\n")
          summary_row[[paste0(comp_name, "_total")]] <- NA
          summary_row[[paste0(comp_name, "_sig")]] <- NA
          summary_row[[paste0(comp_name, "_up")]] <- NA
          summary_row[[paste0(comp_name, "_down")]] <- NA
        }
      }
      
      all_summaries[[coll_name]][[ds_folder]] <-
        as.data.frame(summary_row, stringsAsFactors = FALSE)
    }
    cat("\n")
  }
  
  # ==========================================================================
  # Section 6: Summary Statistics (per collection, per method)
  # ==========================================================================
  
  cat("====================================================================\n")
  cat("Generating Summary Statistics for:", agg_method, "\n")
  cat("====================================================================\n\n")
  
  for (coll_name in names(all_summaries)) {
    cat("--- Summary for", coll_name, "(", agg_method, ") ---\n")
    
    summary_df <- bind_rows(all_summaries[[coll_name]])
    
    # CSV
    fwrite(summary_df, file.path(method_output_dir,
                                 paste0("GSEA_summary_", coll_name, ".csv")))
    
    # XLSX
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
    saveWorkbook(wb_sum, file.path(method_output_dir,
                                   paste0("GSEA_summary_", coll_name, ".xlsx")), overwrite = TRUE)
    
    # Print significant counts only
    sig_cols <- grep("_sig$", colnames(summary_df), value = TRUE)
    print_df <- summary_df[, c("dataset", "cancer_type", sig_cols)]
    print(as.data.frame(print_df), row.names = FALSE)
    cat("\n")
  }
}

cat("====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("  max_abs_t/ : Gene score = strongest phosphosite per gene\n")
cat("  median_t/  : Gene score = median t across all phosphosites\n")
cat("====================================================================\n")

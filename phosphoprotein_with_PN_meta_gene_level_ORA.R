################################################################################
# TP53 Phosphoprotein Gene-Level ORA from Meta-Analysis DPS Results
# (with protein normalization)
#
# Input: META_DPS_*.csv files from
#   phosphoprotein_differential_analysis_subgroup_adjusted/
#   phosphoprotein_DPS_meta_analysis/
#
# For each comparison, extract parent genes from significant phosphosites
# and perform three types of Over-Representation Analysis (ORA):
#   1. positive_ORA: genes with at least one phosphosite "up" & padj < 0.05
#   2. negative_ORA: genes with at least one phosphosite "down" & padj < 0.05
#   3. total_ORA:    genes with at least one phosphosite padj < 0.05 (any direction)
#
# Background gene universe: all parent genes from all phosphosites in that
# comparison (deduplicated).
#
# Parent gene is extracted from gene_symbol_phosphosite column:
#   e.g., "MCM6_S762" -> parent gene = "MCM6"
#
# MSigDB collections (same as protein_GSEA_for_linkedomicskb_subgroup_adjusted.R):
#   Hallmark, C2 (CGP, BIOCARTA, KEGG_LEGACY, KEGG_MEDICUS, PID, REACTOME,
#   WIKIPATHWAYS), C3 (MIR:MIRDB, MIR:MIR_LEGACY, TFT:GTRD, TFT:TFT_LEGACY),
#   C4 (3CA, CGN, CM), C5 (GO:BP, GO:CC, GO:MF), C6 (Oncogenic)
#
# ORA method: clusterProfiler::enricher with custom TERM2GENE and universe
#
# 9 Comparisons:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspot vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspot vs TP53wt,
#   DN vs TP53wt, NonDN vs TP53wt, DN vs NonDN
#
# Output: phosphoprotein_differential_analysis_subgroup_adjusted/
#         phosphoprotein_DPS_meta_analysis/
#         phosphoprotein_meta_gene_level_ORA/
#
# Methodology references:
#   - Boyle EI, et al. GO::TermFinder. Bioinformatics. 2004;20(18):3710-5.
#     PMID: 15297299
#   - Yu G, et al. clusterProfiler: an R package for comparing biological themes
#     among gene clusters. OMICS. 2012;16(5):284-7. PMID: 22455463
#   - Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark
#     gene set collection. Cell Syst. 2015;1(6):417-25. PMID: 26771021
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(msigdbr)
  library(clusterProfiler)
  library(openxlsx)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Input directory: meta-analysis DPS results
meta_dir <- file.path(base_path,
                      "phosphoprotein_differential_analysis_subgroup_adjusted",
                      "phosphoprotein_DPS_meta_analysis")

# Output directory
output_dir <- file.path(meta_dir, "phosphoprotein_meta_gene_level_ORA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use Windows short path names to avoid MAX_PATH (260 char) limitation
# The full path exceeds 260 characters, causing silent dir.create failures
if (.Platform$OS.type == "windows") {
  meta_dir <- shortPathName(meta_dir)
  output_dir <- shortPathName(output_dir)
  cat("  Using short paths to avoid Windows MAX_PATH limitation\n")
  cat("  meta_dir:   ", meta_dir, "\n")
  cat("  output_dir: ", output_dir, "\n\n")
}

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

# Three ORA types
ora_types <- c("positive_ORA", "negative_ORA", "total_ORA")

cat("====================================================================\n")
cat("TP53 Phosphoprotein Gene-Level ORA from Meta-Analysis DPS Results\n")
cat("(with protein normalization)\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Load Gene Sets from MSigDB
# ==============================================================================

cat("Loading MSigDB gene sets (gene_symbol format for ORA)...\n")

# Helper: convert msigdbr output to TERM2GENE data.frame for clusterProfiler
msigdbr_to_term2gene <- function(df) {
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", ]
  df <- df[, c("gs_name", "gene_symbol")]
  colnames(df) <- c("term", "gene")
  df <- df[!duplicated(paste0(df$term, "_", df$gene)), ]
  as.data.frame(df, stringsAsFactors = FALSE)
}

# Same MSigDB collections as protein_GSEA_for_linkedomicskb_subgroup_adjusted.R
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
  df <- df %>% dplyr::select(gs_name, gene_symbol)
  collections[[info$name]] <- msigdbr_to_term2gene(df)
  cat(sprintf("  %-22s : %d gene sets, %d term-gene pairs\n",
              info$name,
              length(unique(collections[[info$name]]$term)),
              nrow(collections[[info$name]])))
}
cat("\n")

# ==============================================================================
# Section 3: Parent Gene Extraction Function
# ==============================================================================

#' Extract parent gene from gene_symbol_phosphosite
#'
#' The gene_symbol_phosphosite column has format "GENE_SITE" (e.g., "MCM6_S762").
#' Parent gene = everything before the last underscore followed by a residue+position
#' pattern (e.g., S/T/Y followed by digits).
#'
#' For multi-underscore names like "TP53BP1_S1758", we use a regex to strip the
#' trailing _[STY][0-9]+ portion.
#'
#' @param gene_symbol_phosphosite Character vector
#' @return Character vector of parent gene symbols
extract_parent_gene <- function(gene_symbol_phosphosite) {
  # Remove trailing _[STY][0-9]+ (phosphosite residue+position)
  parent <- sub("_[STY][0-9]+$", "", gene_symbol_phosphosite)
  parent
}

# ==============================================================================
# Section 4: ORA Function
# ==============================================================================

#' Run ORA using clusterProfiler::enricher
#'
#' @param gene_list Character vector of query genes
#' @param universe Character vector of background genes
#' @param term2gene TERM2GENE data.frame (columns: term, gene)
#' @param pvalueCutoff P-value cutoff for enricher (default 1, to return all results)
#' @param qvalueCutoff Q-value cutoff for enricher (default 1, to return all results)
#' @param minGSSize Minimum gene set size (default 5)
#' @param maxGSSize Maximum gene set size (default 500)
#' @return data.frame of ORA results or NULL
run_ora <- function(gene_list, universe, term2gene,
                    pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = 5, maxGSSize = 500) {
  
  # Remove empty/NA entries
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  universe <- universe[!is.na(universe) & universe != ""]
  
  if (length(gene_list) == 0) {
    cat("      [SKIP] No query genes\n")
    return(NULL)
  }
  
  # Ensure query genes are a subset of universe
  gene_list <- intersect(gene_list, universe)
  
  if (length(gene_list) == 0) {
    cat("      [SKIP] No query genes overlap with universe\n")
    return(NULL)
  }
  
  ora_res <- tryCatch({
    enricher(
      gene         = gene_list,
      universe     = universe,
      TERM2GENE    = term2gene,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize
    )
  }, error = function(e) {
    cat("      [ERROR] enricher failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) {
    return(NULL)
  }
  
  res_df <- as.data.frame(ora_res)
  
  # Sort by p.adjust
  res_df <- res_df[order(res_df$p.adjust), ]
  
  res_df
}

# ==============================================================================
# Section 5: Main Processing Loop
# ==============================================================================

# Storage for summary statistics
all_summaries <- list()
for (coll_name in names(collections)) {
  all_summaries[[coll_name]] <- list()
}

for (j in seq_len(nrow(comparisons))) {
  comp_name <- comparisons$name[j]
  comp_label <- comparisons$label[j]
  
  cat("==================================================================\n")
  cat("Gene-Level ORA for:", comp_name, "\n")
  cat("==================================================================\n")
  
  # Read meta-analysis DPS file
  meta_file <- file.path(meta_dir, paste0("META_DPS_", comp_name, ".csv"))
  if (!file.exists(meta_file)) {
    cat("  [SKIP] META_DPS file not found:", basename(meta_file), "\n\n")
    for (coll_name in names(collections)) {
      all_summaries[[coll_name]][[comp_name]] <- data.frame(
        comparison = comp_name,
        stringsAsFactors = FALSE
      )
    }
    next
  }
  
  meta <- fread(meta_file)
  
  # Verify required columns
  if (!all(c("gene_symbol_phosphosite", "direction", "padj") %in% colnames(meta))) {
    cat("  [SKIP] Missing required columns (gene_symbol_phosphosite, direction, padj)\n\n")
    next
  }
  
  # Filter out rows with missing gene_symbol_phosphosite
  meta <- meta[!is.na(meta$gene_symbol_phosphosite) & meta$gene_symbol_phosphosite != "", ]
  
  # Extract parent gene for all phosphosites
  meta$parent_gene <- extract_parent_gene(meta$gene_symbol_phosphosite)
  
  # Background universe: all unique parent genes in this comparison
  universe_genes <- unique(meta$parent_gene)
  universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
  
  cat("  Total phosphosites:", nrow(meta), "\n")
  cat("  Background universe (unique parent genes):", length(universe_genes), "\n")
  
  # Define the three ORA gene lists
  # 1. positive_ORA: direction == "up" & padj < 0.05
  pos_genes <- unique(meta$parent_gene[meta$direction == "up" & meta$padj < 0.05])
  pos_genes <- pos_genes[!is.na(pos_genes) & pos_genes != ""]
  
  # 2. negative_ORA: direction == "down" & padj < 0.05
  neg_genes <- unique(meta$parent_gene[meta$direction == "down" & meta$padj < 0.05])
  neg_genes <- neg_genes[!is.na(neg_genes) & neg_genes != ""]
  
  # 3. total_ORA: padj < 0.05 (any direction)
  total_genes <- unique(meta$parent_gene[meta$padj < 0.05])
  total_genes <- total_genes[!is.na(total_genes) & total_genes != ""]
  
  cat("  Positive (up, padj < 0.05) parent genes:", length(pos_genes), "\n")
  cat("  Negative (down, padj < 0.05) parent genes:", length(neg_genes), "\n")
  cat("  Total (padj < 0.05) parent genes:", length(total_genes), "\n\n")
  
  # Build a named list of query gene sets for the three ORA types
  ora_gene_lists <- list(
    positive_ORA = pos_genes,
    negative_ORA = neg_genes,
    total_ORA    = total_genes
  )
  
  # Create comparison-level output directory
  comp_out <- file.path(output_dir, comp_name)
  dir.create(comp_out, recursive = TRUE, showWarnings = FALSE)
  
  for (coll_name in names(collections)) {
    term2gene <- collections[[coll_name]]
    cat("  --- Collection:", coll_name, "---\n")
    
    # Create collection output subdirectory
    coll_out <- file.path(comp_out, coll_name)
    dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)
    
    # Summary row for this comparison + collection
    summary_row <- list(comparison = comp_name)
    
    for (ora_type in ora_types) {
      query_genes <- ora_gene_lists[[ora_type]]
      cat("    ", ora_type, ":", sep = "")
      
      if (length(query_genes) == 0) {
        cat(" no query genes, skipping\n")
        summary_row[[paste0(ora_type, "_n_query")]] <- 0
        summary_row[[paste0(ora_type, "_total")]] <- NA
        summary_row[[paste0(ora_type, "_sig")]] <- NA
        next
      }
      
      res <- run_ora(query_genes, universe_genes, term2gene)
      
      if (!is.null(res)) {
        n_sig <- sum(res$p.adjust < 0.05, na.rm = TRUE)
        cat(" ", nrow(res), "sets tested,", n_sig, "significant (padj < 0.05)\n")
        
        # Save CSV
        fwrite(res, file.path(coll_out, paste0("ORA_", ora_type, "_", comp_name, ".csv")))
        
        # Save XLSX with formatting
        wb <- createWorkbook()
        addWorksheet(wb, paste0(comp_label, "_", ora_type))
        writeData(wb, 1, res)
        
        # Header style
        hs <- createStyle(
          textDecoration = "bold", halign = "center",
          border = "bottom", fgFill = "#4472C4", fontColour = "white"
        )
        addStyle(wb, 1, hs, rows = 1, cols = 1:ncol(res), gridExpand = TRUE)
        
        # Highlight rows where p.adjust < 0.05
        padj_col <- which(colnames(res) == "p.adjust")
        if (length(padj_col) == 1) {
          sig_rows <- which(res$p.adjust < 0.05) + 1  # +1 for header row
          if (length(sig_rows) > 0) {
            sig_style <- createStyle(fgFill = "#C6EFCE")
            addStyle(wb, 1, sig_style,
                     rows = sig_rows, cols = 1:ncol(res),
                     gridExpand = TRUE, stack = TRUE)
          }
        }
        
        setColWidths(wb, 1, cols = 1:ncol(res), widths = "auto")
        saveWorkbook(wb, file.path(coll_out, paste0("ORA_", ora_type, "_", comp_name, ".xlsx")),
                     overwrite = TRUE)
        
        summary_row[[paste0(ora_type, "_n_query")]] <- length(query_genes)
        summary_row[[paste0(ora_type, "_total")]] <- nrow(res)
        summary_row[[paste0(ora_type, "_sig")]] <- n_sig
      } else {
        cat(" no enriched terms found\n")
        summary_row[[paste0(ora_type, "_n_query")]] <- length(query_genes)
        summary_row[[paste0(ora_type, "_total")]] <- 0
        summary_row[[paste0(ora_type, "_sig")]] <- 0
      }
    }
    
    all_summaries[[coll_name]][[comp_name]] <-
      as.data.frame(summary_row, stringsAsFactors = FALSE)
  }
  cat("\n")
}

# ==============================================================================
# Section 6: Summary Statistics (per collection)
# ==============================================================================

cat("====================================================================\n")
cat("Generating ORA Summary Statistics\n")
cat("====================================================================\n\n")

for (coll_name in names(all_summaries)) {
  cat("--- Summary for", coll_name, "---\n")
  
  summary_df <- bind_rows(all_summaries[[coll_name]])
  
  if (nrow(summary_df) == 0) {
    cat("  No results\n\n")
    next
  }
  
  # CSV
  fwrite(summary_df, file.path(output_dir,
                               paste0("ORA_summary_", coll_name, ".csv")))
  
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
  saveWorkbook(wb_sum, file.path(output_dir,
                                 paste0("ORA_summary_", coll_name, ".xlsx")), overwrite = TRUE)
  
  # Print summary
  sig_cols <- grep("_sig$", colnames(summary_df), value = TRUE)
  if (length(sig_cols) > 0) {
    print_df <- summary_df[, c("comparison", sig_cols), drop = FALSE]
    print(as.data.frame(print_df), row.names = FALSE)
  }
  cat("\n")
}

cat("====================================================================\n")
cat("Gene-Level ORA Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

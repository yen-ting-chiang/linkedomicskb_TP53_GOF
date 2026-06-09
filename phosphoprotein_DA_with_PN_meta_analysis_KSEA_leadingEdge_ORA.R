################################################################################
# KSEA Leading Edge Substrates Over-Representation Analysis (ORA)
# (Meta-analysis, with protein normalization)
#
# Purpose:
#   For each selected kinase SEPARATELY from the meta-KSEA results, extract
#   its leading edge substrates, convert them to unique parent genes, and
#   perform ORA against MSigDB gene set collections.
#
# Input:
#   1. META_KSEA_*.csv files from
#        phosphoprotein_differential_analysis_subgroup_adjusted/
#        phosphoprotein_DPS_meta_analysis/meta_KSEA/
#      These files contain per-kinase enrichment results with a
#      "leadingEdge_substrates" column (semicolon-separated substrate IDs in
#      GENE_SITE format, e.g., "BUB1B_S670;CDCA5_S33;...")
#
#   2. META_DPS_*.csv files from
#        phosphoprotein_differential_analysis_subgroup_adjusted/
#        phosphoprotein_DPS_meta_analysis/
#      Used to define the background gene universe: all unique parent genes from
#      the "gene_symbol_phosphosite" column.
#
# Kinase selection:
#   Default: CDK1, CDK2, CSNK2A1, AKT1, MAPK1
#   Each kinase is analyzed SEPARATELY (independent ORA per kinase).
#   (user-configurable in Section 1)
#
# Background gene universe:
#   All unique parent genes extracted from gene_symbol_phosphosite in
#   META_DPS_<comparison>.csv (per comparison).
#
# ORA gene list (per kinase):
#   Unique parent genes from the leadingEdge_substrates of that single kinase.
#
# MSigDB collections:
#   Same as phosphoprotein_with_PN_meta_gene_level_ORA.R:
#     Hallmark, C2 (CGP, BIOCARTA, KEGG_LEGACY, KEGG_MEDICUS, PID, REACTOME,
#     WIKIPATHWAYS), C3 (MIR:MIRDB, MIR:MIR_LEGACY, TFT:GTRD, TFT:TFT_LEGACY),
#     C4 (3CA, CGN, CM), C5 (GO:BP, GO:CC, GO:MF), C6 (Oncogenic)
#
# Output:
#   phosphoprotein_differential_analysis_subgroup_adjusted/
#   phosphoprotein_DPS_meta_analysis/
#   meta_KSEA/leadingEdge_ORA/<comparison>/<kinase>/<collection>/
#
# ORA method:
#   clusterProfiler::enricher with custom TERM2GENE and background universe
#
# Methodology references:
#   - Boyle EI, et al. GO::TermFinder. Bioinformatics. 2004;20(18):3710-5.
#     PMID: 15297299
#   - Yu G, et al. clusterProfiler: an R package for comparing biological themes
#     among gene clusters. OMICS. 2012;16(5):284-7. PMID: 22455463
#   - Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark
#     gene set collection. Cell Syst. 2015;1(6):417-25. PMID: 26771021
#   - Casado P, et al. Kinase-substrate enrichment analysis provides insights
#     into the heterogeneity of signaling pathway activation in leukemia cells.
#     Sci Signal. 2013;6(268):rs6. PMID: 23532336
################################################################################

# ==============================================================================
# Section 1: Setup and User Configuration
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(msigdbr)
  library(clusterProfiler)
  library(openxlsx)
  library(ggplot2)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

# Meta-analysis DPS directory (for background gene universe)
meta_dps_dir <- file.path(base_path,
                          "phosphoprotein_differential_analysis_subgroup_adjusted",
                          "phosphoprotein_DPS_meta_analysis")

# Meta-KSEA input directory
meta_ksea_dir <- file.path(meta_dps_dir, "meta_KSEA")

# Output directory
output_dir <- file.path(meta_ksea_dir, "leadingEdge_ORA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use Windows short path names to avoid MAX_PATH (260 char) limitation
if (.Platform$OS.type == "windows") {
  meta_dps_dir <- shortPathName(meta_dps_dir)
  meta_ksea_dir <- shortPathName(meta_ksea_dir)
  output_dir <- shortPathName(output_dir)
  cat("  Using short paths to avoid Windows MAX_PATH limitation\n")
  cat("  meta_dps_dir:  ", meta_dps_dir, "\n")
  cat("  meta_ksea_dir: ", meta_ksea_dir, "\n")
  cat("  output_dir:    ", output_dir, "\n\n")
}

# ---- User-configurable: Kinases of interest ----
# Each kinase is analyzed SEPARATELY (independent ORA per kinase).
selected_kinases <- c("CDK1", "CDK2", "CSNK2A1", "AKT1", "MAPK1")

# ---- User-configurable: Comparisons to run ----
# Default: all 9 comparisons.
comparisons <- data.frame(
  name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
           "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
           "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
  label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
            "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
            "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
  stringsAsFactors = FALSE
)

# To run only specific comparisons, uncomment and modify the line below:
# comparisons <- comparisons[comparisons$name %in% c("MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt"), ]

# ---- User-configurable: Collections to run ----
# Set collections_to_use to a character vector of collection names to limit,
# or leave as NULL to run all 19 collections.
# Available: "Hallmark", "C2_CGP", "C2_CP_BIOCARTA", "C2_CP_KEGG_LEGACY",
#            "C2_CP_KEGG_MEDICUS", "C2_CP_PID", "C2_CP_REACTOME",
#            "C2_CP_WIKIPATHWAYS", "C3_MIR_MIRDB", "C3_MIR_MIR_LEGACY",
#            "C3_TFT_GTRD", "C3_TFT_TFT_LEGACY", "C4_3CA", "C4_CGN",
#            "C4_CM", "C5_GO_BP", "C5_GO_CC", "C5_GO_MF", "C6_Oncogenic"
collections_to_use <- NULL  # NULL = use all collections

# To run only specific collections, uncomment and modify the line below:
# collections_to_use <- c("Hallmark", "C2_CP_REACTOME", "C5_GO_BP", "C6_Oncogenic")

# ---- User-configurable: Bubble Plot Settings ----
# Which comparison(s), kinase(s), and collection(s) to plot bubble plots for.
# Default: TP53mt_vs_TP53wt comparison, CDK1 and CDK2 kinases, Hallmark collection.
bubble_comparisons <- c("TP53mt_vs_TP53wt")
bubble_kinases     <- c("CDK1", "CDK2")
bubble_collections <- c("Hallmark")

# To plot all combinations, uncomment the following lines:
# bubble_comparisons <- comparisons$name
# bubble_kinases     <- selected_kinases
# bubble_collections <- NULL  # NULL = all collections used in ORA

# Bubble plot color palette
# Select a color gradient for bubble plot p.adjust mapping.
# The gradient runs from low p.adjust (most significant) to high p.adjust.
#
# Available palettes:
#   1 = NPG Red to Blue         : low="#E64B35", high="#4DBBD5" (default)
#   2 = Lancet Red to Teal      : low="#ED0000", high="#0099B4"
#   3 = NEJM Red to Blue        : low="#BC3C29", high="#0072B5"
#   4 = JAMA Warm to Cool       : low="#DF8F44", high="#374E55"
#   5 = JCO Gold to Blue        : low="#EFC000", high="#0073C2"
#   6 = Viridis-inspired        : low="#FDE725", high="#440154"
#   7 = Magma-inspired          : low="#FCFDBF", high="#000004"
#   8 = Nature Red to Green     : low="#E64B35", high="#00A087"
palette_bubble_choice <- 4  # Change this number to select bubble palette (1-8)

bubble_palettes <- list(
  c("#E64B35", "#4DBBD5"),  # 1: NPG Red to Blue
  c("#ED0000", "#0099B4"),  # 2: Lancet Red to Teal
  c("#BC3C29", "#0072B5"),  # 3: NEJM Red to Blue
  c("#DF8F44", "#374E55"),  # 4: JAMA Warm to Cool
  c("#EFC000", "#0073C2"),  # 5: JCO Gold to Blue
  c("#FDE725", "#440154"),  # 6: Viridis-inspired
  c("#FCFDBF", "#000004"),  # 7: Magma-inspired
  c("#E64B35", "#00A087")   # 8: Nature Red to Green
)

bubble_color_low  <- bubble_palettes[[palette_bubble_choice]][1]
bubble_color_high <- bubble_palettes[[palette_bubble_choice]][2]

# Bubble plot output directory
bubble_output_dir <- file.path(base_path,
                               "phosphoprotein_differential_analysis_subgroup_adjusted",
                               "phosphoprotein_DPS_meta_analysis",
                               "meta_KSEA", "leadingEdge_ORA", "ORA_bubble_plot")
dir.create(bubble_output_dir, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows") {
  bubble_output_dir <- shortPathName(bubble_output_dir)
}

cat("====================================================================\n")
cat("KSEA Leading Edge Substrates ORA (per kinase, separately)\n")
cat("(Meta-analysis, with protein normalization)\n")
cat("====================================================================\n\n")

cat("Selected kinases (each analyzed separately):",
    paste(selected_kinases, collapse = ", "), "\n")
cat("Comparisons to process:", nrow(comparisons), "\n\n")

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

# Same MSigDB collections as phosphoprotein_with_PN_meta_gene_level_ORA.R
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

# Apply user collection filter
if (!is.null(collections_to_use)) {
  valid <- intersect(collections_to_use, names(collections))
  if (length(valid) == 0) {
    stop("None of the specified collections_to_use match available collections.")
  }
  collections <- collections[valid]
  cat("\n  Filtered to user-selected collections:", paste(valid, collapse = ", "), "\n")
}

cat("\n")

# ==============================================================================
# Section 3: Helper Functions
# ==============================================================================

#' Extract parent gene from substrate ID (GENE_SITE format)
#'
#' The substrate ID has format "GENE_SITE" (e.g., "MCM6_S762", "TP53BP1_S1678").
#' Parent gene = everything before the last underscore followed by a residue+position
#' pattern ([STY] followed by digits).
#'
#' @param substrate_id Character vector of substrate IDs
#' @return Character vector of parent gene symbols
extract_parent_gene <- function(substrate_id) {
  sub("_[STY][0-9]+$", "", substrate_id)
}

#' Parse leading edge substrates from KSEA result
#'
#' The leadingEdge_substrates column contains semicolon-separated substrate IDs.
#' This function splits them and returns a character vector.
#'
#' @param le_string A single character string (semicolon-separated)
#' @return Character vector of substrate IDs (or empty if input is NA/"_"/empty)
parse_leading_edge <- function(le_string) {
  if (is.na(le_string) || le_string == "" || le_string == "_") {
    return(character(0))
  }
  trimws(unlist(strsplit(le_string, ";")))
}

#' Run ORA using clusterProfiler::enricher
#'
#' @param gene_list Character vector of query genes
#' @param universe Character vector of background genes
#' @param term2gene TERM2GENE data.frame (columns: term, gene)
#' @param pvalueCutoff P-value cutoff for enricher (default 1, return all)
#' @param qvalueCutoff Q-value cutoff for enricher (default 1, return all)
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
  res_df <- res_df[order(res_df$p.adjust), ]
  res_df
}

# ==============================================================================
# Section 4: Main Processing Loop (per comparison x per kinase)
# ==============================================================================

# Storage for summary statistics: keyed by collection -> list of rows
all_summaries <- list()
for (coll_name in names(collections)) {
  all_summaries[[coll_name]] <- list()
}

for (j in seq_len(nrow(comparisons))) {
  comp_name <- comparisons$name[j]
  comp_label <- comparisons$label[j]
  
  cat("==================================================================\n")
  cat("Leading Edge ORA for:", comp_name, "\n")
  cat("==================================================================\n")
  
  # ------------------------------------------------------------------
  # Step 1: Read META_KSEA file
  # ------------------------------------------------------------------
  ksea_file <- file.path(meta_ksea_dir, paste0("META_KSEA_", comp_name, ".csv"))
  if (!file.exists(ksea_file)) {
    cat("  [SKIP] META_KSEA file not found:", basename(ksea_file), "\n\n")
    next
  }
  
  ksea <- fread(ksea_file)
  
  # Verify required columns
  if (!all(c("Kinase.Gene", "leadingEdge_substrates") %in% colnames(ksea))) {
    cat("  [SKIP] Missing required columns (Kinase.Gene, leadingEdge_substrates)\n\n")
    next
  }
  
  # ------------------------------------------------------------------
  # Step 2: Build background gene universe from META_DPS file
  # ------------------------------------------------------------------
  dps_files <- list.files(meta_dps_dir,
                          pattern = glob2rx(paste0("META_DPS_", comp_name, "*.csv")),
                          full.names = TRUE)
  
  # If glob finds nothing, try exact match
  if (length(dps_files) == 0) {
    dps_exact <- file.path(meta_dps_dir, paste0("META_DPS_", comp_name, ".csv"))
    if (file.exists(dps_exact)) {
      dps_files <- dps_exact
    }
  }
  
  if (length(dps_files) == 0) {
    cat("  [SKIP] META_DPS file not found for:", comp_name, "\n\n")
    next
  }
  
  # Read all matching META_DPS files and collect gene_symbol_phosphosite
  all_bg_phosphosites <- character(0)
  for (dps_file in dps_files) {
    dps <- fread(dps_file, select = "gene_symbol_phosphosite")
    all_bg_phosphosites <- c(all_bg_phosphosites,
                             dps$gene_symbol_phosphosite[!is.na(dps$gene_symbol_phosphosite) &
                                                           dps$gene_symbol_phosphosite != ""])
  }
  
  # Convert to unique parent genes for background
  universe_genes <- unique(extract_parent_gene(all_bg_phosphosites))
  universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
  
  cat("  Background universe (unique parent genes from META_DPS):", length(universe_genes), "\n\n")
  
  # ------------------------------------------------------------------
  # Step 3: Loop over each kinase SEPARATELY
  # ------------------------------------------------------------------
  for (kinase_name in selected_kinases) {
    
    cat("  *** Kinase:", kinase_name, "***\n")
    
    # Check if this kinase exists in KSEA results
    ksea_row <- ksea[ksea$Kinase.Gene == kinase_name, ]
    
    if (nrow(ksea_row) == 0) {
      cat("    [SKIP] Kinase not found in KSEA results\n\n")
      # Record empty summary for each collection
      for (coll_name in names(collections)) {
        all_summaries[[coll_name]][[paste0(comp_name, "__", kinase_name)]] <-
          data.frame(comparison = comp_name, kinase = kinase_name,
                     n_query_genes = 0, total_sets = NA, sig_sets = NA,
                     stringsAsFactors = FALSE)
      }
      next
    }
    
    # Parse leading edge substrates for this kinase
    le_str <- as.character(ksea_row$leadingEdge_substrates[1])
    le_subs <- parse_leading_edge(le_str)
    cat("    Leading edge substrates:", length(le_subs), "\n")
    
    if (length(le_subs) == 0) {
      cat("    [SKIP] No leading edge substrates\n\n")
      for (coll_name in names(collections)) {
        all_summaries[[coll_name]][[paste0(comp_name, "__", kinase_name)]] <-
          data.frame(comparison = comp_name, kinase = kinase_name,
                     n_query_genes = 0, total_sets = NA, sig_sets = NA,
                     stringsAsFactors = FALSE)
      }
      next
    }
    
    # Convert to unique parent genes
    le_parent_genes <- unique(extract_parent_gene(le_subs))
    le_parent_genes <- le_parent_genes[!is.na(le_parent_genes) & le_parent_genes != ""]
    cat("    Unique leading edge parent genes:", length(le_parent_genes), "\n")
    
    # Check overlap with universe
    overlap_n <- length(intersect(le_parent_genes, universe_genes))
    cat("    Query genes overlapping with universe:", overlap_n, "/",
        length(le_parent_genes), "\n")
    
    if (overlap_n == 0) {
      cat("    [SKIP] No query genes overlap with background universe\n\n")
      for (coll_name in names(collections)) {
        all_summaries[[coll_name]][[paste0(comp_name, "__", kinase_name)]] <-
          data.frame(comparison = comp_name, kinase = kinase_name,
                     n_query_genes = length(le_parent_genes),
                     total_sets = 0, sig_sets = 0,
                     stringsAsFactors = FALSE)
      }
      next
    }
    
    # Create kinase-level output directory: <comparison>/<kinase>/
    kinase_out <- file.path(output_dir, comp_name, kinase_name)
    dir.create(kinase_out, recursive = TRUE, showWarnings = FALSE)
    
    # Save the leading edge gene list for this kinase
    le_info <- data.frame(
      substrate = le_subs,
      parent_gene = extract_parent_gene(le_subs),
      stringsAsFactors = FALSE
    )
    fwrite(le_info, file.path(kinase_out,
                              paste0("leadingEdge_substrates_", kinase_name, "_", comp_name, ".csv")))
    
    # Run ORA for each collection
    for (coll_name in names(collections)) {
      term2gene <- collections[[coll_name]]
      cat("    --- Collection:", coll_name, "---\n")
      
      # Create collection output subdirectory
      coll_out <- file.path(kinase_out, coll_name)
      dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)
      
      # Summary row
      summary_row <- list(comparison = comp_name,
                          kinase = kinase_name,
                          n_query_genes = length(le_parent_genes))
      
      cat("      ORA:", sep = "")
      
      res <- run_ora(le_parent_genes, universe_genes, term2gene)
      
      if (!is.null(res)) {
        n_sig <- sum(res$p.adjust < 0.05, na.rm = TRUE)
        cat(" ", nrow(res), "sets tested,", n_sig,
            "significant (padj < 0.05)\n")
        
        # Save CSV
        fwrite(res, file.path(coll_out,
                              paste0("ORA_", kinase_name, "_", comp_name, ".csv")))
        
        # Save XLSX with formatting
        wb <- createWorkbook()
        sht_name <- paste0(kinase_name, "_", comp_label)
        # Truncate sheet name to 31 characters (Excel limit)
        if (nchar(sht_name) > 31) sht_name <- substr(sht_name, 1, 31)
        addWorksheet(wb, sht_name)
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
        saveWorkbook(wb, file.path(coll_out,
                                   paste0("ORA_", kinase_name, "_", comp_name, ".xlsx")),
                     overwrite = TRUE)
        
        summary_row[["total_sets"]] <- nrow(res)
        summary_row[["sig_sets"]] <- n_sig
      } else {
        cat(" no enriched terms found\n")
        summary_row[["total_sets"]] <- 0
        summary_row[["sig_sets"]] <- 0
      }
      
      all_summaries[[coll_name]][[paste0(comp_name, "__", kinase_name)]] <-
        as.data.frame(summary_row, stringsAsFactors = FALSE)
    }
    
    cat("\n")
  }
  cat("\n")
}

# ==============================================================================
# Section 5: Summary Statistics (per collection)
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
                               paste0("ORA_leadingEdge_summary_", coll_name, ".csv")))
  
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
                                 paste0("ORA_leadingEdge_summary_", coll_name, ".xlsx")), overwrite = TRUE)
  
  # Print summary table
  print(as.data.frame(summary_df), row.names = FALSE)
  cat("\n")
}

cat("====================================================================\n")
cat("KSEA Leading Edge ORA Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 6: Bubble Plots for Selected ORA Results
# ==============================================================================

cat("====================================================================\n")
cat("Generating ORA Bubble Plots\n")
cat("====================================================================\n\n")

cat("Bubble plot comparisons:", paste(bubble_comparisons, collapse = ", "), "\n")
cat("Bubble plot kinases:", paste(bubble_kinases, collapse = ", "), "\n")
cat("Bubble plot palette:", palette_bubble_choice, "->",
    bubble_color_low, "(low) to", bubble_color_high, "(high)\n")
cat("Bubble plot output:", bubble_output_dir, "\n\n")

# Determine which collections to plot
if (is.null(bubble_collections)) {
  bubble_collections_final <- names(collections)
} else {
  bubble_collections_final <- intersect(bubble_collections, names(collections))
}
cat("Bubble plot collections:", paste(bubble_collections_final, collapse = ", "), "\n\n")

# Prefix-to-strip mapping for cleaner pathway names in bubble plots
collection_prefix_map <- list(
  Hallmark           = "(?i)^hallmark[_ ]",
  C2_CGP             = "(?i)^CGP[_ ]",
  C2_CP_BIOCARTA     = "(?i)^biocarta[_ ]",
  C2_CP_KEGG_LEGACY  = "(?i)^kegg[_ ]",
  C2_CP_KEGG_MEDICUS = "(?i)^kegg[_ ]",
  C2_CP_PID          = "(?i)^pid[_ ]",
  C2_CP_REACTOME     = "(?i)^reactome[_ ]",
  C2_CP_WIKIPATHWAYS = "(?i)^wp[_ ]",
  C3_MIR_MIRDB       = "(?i)^mir[_ ]",
  C3_MIR_MIR_LEGACY  = "(?i)^mir[_ ]",
  C3_TFT_GTRD        = "(?i)^gtrd[_ ]",
  C3_TFT_TFT_LEGACY  = "(?i)^tft[_ ]",
  C4_3CA             = "(?i)^3ca[_ ]",
  C4_CGN             = "(?i)^cgn[_ ]",
  C4_CM              = "(?i)^cm[_ ]",
  C5_GO_BP           = "(?i)^gobp[_ ]",
  C5_GO_CC           = "(?i)^gocc[_ ]",
  C5_GO_MF           = "(?i)^gomf[_ ]",
  C6_Oncogenic       = "(?i)^oncogenic[_ ]"
)

bubble_count <- 0

for (bp_comp in bubble_comparisons) {
  for (bp_kinase in bubble_kinases) {
    for (bp_coll in bubble_collections_final) {
      
      cat("  Bubble plot:", bp_comp, "/", bp_kinase, "/", bp_coll, "\n")
      
      # Locate the ORA CSV result file
      ora_csv <- file.path(output_dir, bp_comp, bp_kinase, bp_coll,
                           paste0("ORA_", bp_kinase, "_", bp_comp, ".csv"))
      
      if (!file.exists(ora_csv)) {
        cat("    [SKIP] ORA result file not found:", basename(ora_csv), "\n")
        next
      }
      
      res_df <- fread(ora_csv)
      
      if (nrow(res_df) == 0) {
        cat("    [SKIP] ORA result file is empty\n")
        next
      }
      
      # Calculate numeric GeneRatio for plotting
      res_df$GeneRatioNum <- sapply(res_df$GeneRatio, function(x) {
        parts <- as.numeric(strsplit(as.character(x), "/")[[1]])
        if (length(parts) == 2 && parts[2] != 0) parts[1] / parts[2] else 0
      })
      
      # Sort by p.adjust
      res_df <- res_df[order(res_df$p.adjust), ]
      
      # Identify significant rows (padj < 0.05)
      sig_df <- res_df[res_df$p.adjust < 0.05, ]
      
      # If significant pathways < 5, pad to 5. Otherwise keep up to 15.
      if (nrow(sig_df) < 5) {
        plot_df <- head(res_df, 5)
      } else {
        plot_df <- head(sig_df, 15)
      }
      
      if (nrow(plot_df) == 0) {
        cat("    [SKIP] No data for bubble plot\n")
        next
      }
      
      # Clean pathway names: remove collection prefix for cleaner plotting
      prefix_pattern <- collection_prefix_map[[bp_coll]]
      if (!is.null(prefix_pattern)) {
        plot_df$Description <- gsub(prefix_pattern, "", plot_df$Description)
      }
      
      # Replace underscores with spaces for readability
      plot_df$Description <- gsub("_", " ", plot_df$Description)
      
      # Order Description factor by p.adjust (descending so smallest padj at top)
      plot_df$Description <- factor(plot_df$Description,
                                    levels = rev(plot_df$Description))
      
      # Control plot title visibility
      show_plot_title <- TRUE
      plot_title_text <- if (show_plot_title) {
        paste0(bp_coll, " ORA: ", bp_kinase, " (", bp_comp, ")")
      } else {
        NULL
      }
      
      # Bubble plot (custom ggplot, same style as reference script)
      p <- ggplot(plot_df, aes(x = GeneRatioNum, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = bubble_color_low, high = bubble_color_high,
                             name = "p.adjust") +
        labs(title = plot_title_text, x = "GeneRatio", y = "") +
        scale_x_continuous(expand = expansion(mult = 0.15)) +
        scale_y_discrete(expand = expansion(mult = 0.1)) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(color = "black"),
          axis.text.y = element_text(size = 7, face = "bold", color = "black"),
          axis.text.x = element_text(size = 7, face = "bold", color = "black"),
          axis.title = element_text(size = 7, face = "bold", color = "black"),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5,
                                    color = "black"),
          legend.position = "right",
          legend.box.margin = margin(0, 0, 0, 10),
          legend.title = element_text(size = 7, face = "bold", color = "black"),
          legend.text = element_text(size = 7, color = "black"),
          legend.key.size = unit(0.3, "cm")
        )
      
      if (!show_plot_title) {
        p <- p + theme(plot.title = element_blank())
      }
      
      # Create comparison-level subdirectory within bubble_output_dir
      bp_out_subdir <- file.path(bubble_output_dir, bp_comp)
      dir.create(bp_out_subdir, recursive = TRUE, showWarnings = FALSE)
      
      file_prefix <- paste0("ORA_bubble_", bp_kinase, "_", bp_coll, "_", bp_comp)
      
      pdf(file.path(bp_out_subdir, paste0(file_prefix, ".pdf")),
          width = 3.8, height = 2.5)
      print(p)
      dev.off()
      
      tiff(file.path(bp_out_subdir, paste0(file_prefix, ".tiff")),
           width = 3.8, height = 2.5, units = "in", res = 300,
           compression = "lzw")
      print(p)
      dev.off()
      
      cat("    Bubble plot saved: ", file_prefix, ".pdf/.tiff\n")
      bubble_count <- bubble_count + 1
    }
  }
}

cat("\nTotal bubble plots generated:", bubble_count, "\n")
cat("Bubble plots saved to:", bubble_output_dir, "\n")

cat("====================================================================\n")
cat("KSEA Leading Edge ORA + Bubble Plot Pipeline Complete!\n")
cat("ORA results saved to:", output_dir, "\n")
cat("Bubble plots saved to:", bubble_output_dir, "\n")
cat("====================================================================\n")

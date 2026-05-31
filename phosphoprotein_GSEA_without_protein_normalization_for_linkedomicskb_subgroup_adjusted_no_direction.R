################################################################################
# TP53 Phosphoprotein-Level Preranked GSEA (LinkedOmicsKB Version)
#
# Uses moderated-t statistics from phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted/ as ranking
#
# Gene Set Collection:
#   - PTMsigDB (phosphosite level)
#     data_PTMsigDB_all_sites_v2.0.0.xlsx
#
# ID Mapping Strategy:
#   DPS phosphosite ID format:
#     ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
#     (ENSG | ENSP | residue+position | 15-mer flanking sequence | multiplicity)
#
#   PTMsigDB site.annotation format:
#     PPP1R12A_T696:15226371;20801872
#     (GENE_SYMBOL_RESIDUE+POSITION:PMIDs)
#
#   Matching approach: Use the 15-mer flanking sequence (field 4 of DPS ID)
#   to look up the corresponding PTMsigDB gene_site ID (from site.annotation).
#   PTMsigDB also contains flanking sequences (site.flanking column) that can
#   be matched to the DPS flanking sequences. This enables cross-referencing
#   between Ensembl-based phosphosite IDs and HGNC gene symbol-based PTMsigDB
#   identifiers.
#
# 10 Datasets: BRCA, CCRCC, COAD, GBM, HNSCC, LSCC, LUAD, OV, PDAC, UCEC
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspots vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspots vs TP53wt,
#   DN vs TP53wt, Non-DN vs TP53wt, DN vs non-DN
#
# Methodology: fgsea preranked GSEA with PTMsigDB
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(fgsea)
    library(openxlsx)
    library(readxl)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

datasets <- data.frame(
    folder = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
               "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    cancer_type = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
                    "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    stringsAsFactors = FALSE
)

# Input: phosphoprotein DPS results (without protein normalization)
dps_dir <- file.path(base_path, "phosphoprotein_differential_analysis_without_protein_normalization_subgroup_adjusted")

# Output: GSEA results (without protein normalization)
output_dir <- file.path(base_path, "phosphoprotein_GSEA_without_protein_normalization_subgroup_adjusted")
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

# fgsea parameters
minSize <- 5
maxSize <- 500

cat("====================================================================\n")
cat("TP53 Phosphoprotein-Level Preranked GSEA (LinkedOmicsKB - without protein normalization)\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Load Gene Sets from PTMsigDB
# ==============================================================================

cat("Loading PTMsigDB gene sets...\n")

#' Read PTMsigDB and build pathway lists keyed by GENE_SITE identifier
#' Also build a flanking-sequence-to-gene_site lookup table for ID mapping
#'
#' @param fp Path to PTMsigDB xlsx file
#' @return A list with two elements:
#'   - pathways: named list of character vectors (signature -> gene_site IDs)
#'   - flank_lookup: named character vector (flanking_seq -> gene_site)
read_ptmsigdb <- function(fp) {
    if (!file.exists(fp)) {
        stop("[PTMsigDB] File does not exist: ", fp)
    }

    df <- readxl::read_xlsx(fp)
    req <- c("signature", "site.annotation", "site.flanking")
    if (!all(req %in% names(df))) {
        stop("[PTMsigDB] xlsx is missing necessary fields: ",
             paste(setdiff(req, names(df)), collapse = ", "))
    }

    # Extract gene_site from site.annotation (e.g., "PPP1R12A_T696:15226371" -> "PPP1R12A_T696")
    gene_site <- toupper(sub(":.*$", "", trimws(df$site.annotation)))

    # Extract flanking sequence (15-mer)
    flanking <- toupper(trimws(df$site.flanking))

    # Build pathway lists keyed by gene_site
    by_sig <- split(gene_site, df$signature)
    pathways <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
    pathways[lengths(pathways) == 0] <- NULL

    # Build flanking-to-gene_site lookup
    # Use unique flanking sequences; for duplicates, keep first occurrence
    valid <- nzchar(flanking) & nzchar(gene_site) & nchar(flanking) == 15
    flank_lookup <- setNames(gene_site[valid], flanking[valid])
    flank_lookup <- flank_lookup[!duplicated(names(flank_lookup))]

    list(pathways = pathways, flank_lookup = flank_lookup)
}

ptmsigdb_file <- file.path(base_path, "data_PTMsigDB_all_sites_v2.0.0.xlsx")
ptmsigdb_data <- read_ptmsigdb(ptmsigdb_file)
pathways_ptm <- ptmsigdb_data$pathways
flank_lookup <- ptmsigdb_data$flank_lookup

cat("  PTMsigDB:", length(pathways_ptm), "phosphosite sets\n")
cat("  Flanking lookup table:", length(flank_lookup), "unique 15-mer entries\n\n")

# Named list for iteration
collections <- list(
    PTMsigDB = pathways_ptm
)

# ==============================================================================
# Section 3: Run Preranked GSEA
# ==============================================================================

#' Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
#'
#' DPS ID format: ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
#' This function extracts the 15-mer flanking sequence (field 4) and uses the
#' prebuilt lookup table to find the corresponding PTMsigDB gene_site ID
#' (e.g., "PPP1R12A_T696").
#'
#' @param dps_ids Character vector of DPS phosphosite IDs
#' @param lookup Named character vector (flanking -> gene_site)
#' @return Character vector of PTMsigDB gene_site IDs (NA if no match)
map_dps_to_ptmsigdb <- function(dps_ids, lookup) {
    parts <- strsplit(dps_ids, "|", fixed = TRUE)

    # Extract flanking sequence (field 4) and convert to uppercase
    flanking <- vapply(parts, function(x) {
        if (length(x) >= 4) toupper(x[4]) else NA_character_
    }, character(1))

    # Look up gene_site via flanking sequence
    mapped <- lookup[flanking]

    # Report mapping statistics
    n_total <- length(dps_ids)
    n_mapped <- sum(!is.na(mapped))
    cat(sprintf("    ID mapping: %d of %d phosphosites mapped to PTMsigDB (%.1f%%)\n",
                n_mapped, n_total, 100 * n_mapped / n_total))

    mapped
}

#' Run fgsea preranked GSEA from a DPS table
#' @param dps_file Path to DPS CSV file
#' @param pathways Named list of phosphosite sets (gene_site IDs)
#' @param lookup Flanking-to-gene_site lookup table
#' @return fgsea result data.frame or NULL
run_preranked_gsea <- function(dps_file, pathways, lookup) {
    if (!file.exists(dps_file)) return(NULL)

    dps <- fread(dps_file)
    if (!("phosphosite" %in% colnames(dps)) || !("t" %in% colnames(dps))) {
        cat("    [WARN] Missing phosphosite/t columns in:", basename(dps_file), "\n")
        return(NULL)
    }

    # Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
    gene_site_ids <- map_dps_to_ptmsigdb(dps$phosphosite, lookup)

    # Keep only successfully mapped phosphosites
    mapped_mask <- !is.na(gene_site_ids)
    if (sum(mapped_mask) < 50) {
        cat("    [WARN] Too few mapped phosphosites:", sum(mapped_mask), "\n")
        return(NULL)
    }

    mapped_ids <- gene_site_ids[mapped_mask]
    mapped_t <- dps$t[mapped_mask]

    # Build named vector of moderated-t statistics
    stats_vec <- setNames(mapped_t, mapped_ids)

    # Remove duplicates (keep first = most significant from topTable)
    stats_vec <- stats_vec[!duplicated(names(stats_vec))]

    # Remove non-finite values
    stats_vec <- stats_vec[is.finite(stats_vec)]

    if (length(stats_vec) < 50) {
        cat("    [WARN] Too few unique mapped phosphosites after dedup:", length(stats_vec), "\n")
        return(NULL)
    }

    # Sort descending (required by fgsea)
    stats_vec <- sort(stats_vec, decreasing = TRUE)

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
    res$leadingEdge <- vapply(res$leadingEdge,
                              function(x) paste(x, collapse = ";"),
                              character(1))

    # Sort by padj
    res <- res[order(res$padj), ]
    as.data.frame(res)
}

# ==============================================================================
# Section 4: Main Processing Loop
# ==============================================================================

# Storage for summary statistics (one per collection)
all_summaries <- list()
for (coll_name in names(collections)) {
    all_summaries[[coll_name]] <- list()
}

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    ds_deg_dir <- file.path(dps_dir, ds_folder)
    if (!dir.exists(ds_deg_dir)) {
        cat("  [SKIP] DPS directory not found\n\n")
        next
    }

    ds_out <- file.path(output_dir, ds_folder)
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

            # DPS results saved as DPS_{comp_name}.csv
            dps_file <- file.path(ds_deg_dir, paste0("DPS_", comp_name, ".csv"))
            cat("    ", comp_label, ":", sep = "")

            res <- run_preranked_gsea(dps_file, pathways, flank_lookup)

            if (!is.null(res)) {
                n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
                n_up <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
                n_down <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
                cat(" ", nrow(res), "sets tested,", n_sig, "significant\n")

                # Save CSV
                fwrite(res, file.path(coll_out, paste0("GSEA_", comp_name, ".csv")))

                # Save XLSX
                wb <- createWorkbook()
                addWorksheet(wb, comp_label)
                writeData(wb, 1, res)
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

# ==============================================================================
# Section 5: Summary Statistics (per collection)
# ==============================================================================

cat("====================================================================\n")
cat("Generating Summary Statistics\n")
cat("====================================================================\n\n")

for (coll_name in names(all_summaries)) {
    cat("--- Summary for", coll_name, "---\n")

    summary_df <- bind_rows(all_summaries[[coll_name]])

    # CSV
    fwrite(summary_df, file.path(output_dir,
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
    saveWorkbook(wb_sum, file.path(output_dir,
        paste0("GSEA_summary_", coll_name, ".xlsx")), overwrite = TRUE)

    # Print significant counts only
    sig_cols <- grep("_sig$", colnames(summary_df), value = TRUE)
    print_df <- summary_df[, c("dataset", "cancer_type", sig_cols)]
    print(as.data.frame(print_df), row.names = FALSE)
    cat("\n")
}

cat("====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

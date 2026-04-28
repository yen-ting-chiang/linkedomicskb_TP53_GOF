################################################################################
# TP53 Phosphoprotein KSEA (Kinase-Substrate Enrichment Analysis)
#
# Uses moderated-t statistics from phosphoprotein_differential_analysis_without_protein_normalization/ as ranking
#
# KSEA Package: KSEAapp
#
# ID Mapping Strategy:
#   DPS phosphosite ID format:
#     ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
#
#   KSEAapp requires Gene Symbols (e.g., CHEK2) and Residue (e.g., S260).
#   We use PTMsigDB (data_PTMsigDB_all_sites_v2.0.0.xlsx) as a bridge:
#   1. Map DPS ID -> 15-mer flanking sequence.
#   2. Map flanking sequence -> PTMsigDB gene_site (e.g., CHEK2_S260) via lookup.
#   3. Extract Gene and Residue from gene_site for KSEA.
#
# 10 Datasets: BRCA, CCRCC, COAD, GBM, HNSCC, LSCC, LUAD, OV, PDAC, UCEC
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspots vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspots vs TP53wt,
#   DN vs TP53wt, Non-DN vs TP53wt, DN vs non-DN
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(KSEAapp)
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
dps_dir <- file.path(base_path, "phosphoprotein_differential_analysis_without_protein_normalization")

# Output: KSEA results (without protein normalization)
output_dir <- file.path(base_path, "phosphoprotein_KSEA_without_protein_normalization")
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

# KSEA Parameters
use_networkin <- FALSE
networkin_cutoff <- 5

cat("====================================================================\n")
cat("TP53 Phosphoprotein KSEA (KSEAapp - LinkedOmicsKB - without protein normalization)\n")
cat("====================================================================\n\n")

# Load KSData from KSEAapp namespace
data("KSData", package = "KSEAapp", envir = environment())
cat("  Loaded KSData from KSEAapp package\n")

# ==============================================================================
# Section 2: Load Mapping Data (PTMsigDB)
# ==============================================================================

cat("Loading PTMsigDB for ID mapping...\n")
ptmsigdb_file <- file.path(base_path, "data_PTMsigDB_all_sites_v2.0.0.xlsx")

load_flank_lookup <- function(fp) {
    if (!file.exists(fp)) stop("PTMsigDB file not found: ", fp)
    df <- readxl::read_xlsx(fp)
    # site.annotation: e.g. "PPP1R12A_T696:15226371"
    # site.flanking: 15-mer sequence
    gene_site <- toupper(sub(":.*$", "", trimws(df$site.annotation)))
    flanking <- toupper(trimws(df$site.flanking))
    
    # Build lookup
    valid <- nzchar(flanking) & nzchar(gene_site) & nchar(flanking) == 15
    lookup <- setNames(gene_site[valid], flanking[valid])
    lookup <- lookup[!duplicated(names(lookup))]
    lookup
}

flank_lookup <- load_flank_lookup(ptmsigdb_file)
cat("  Flanking lookup table loaded:", length(flank_lookup), "entries\n\n")

# ==============================================================================
# Section 3: Run KSEA
# ==============================================================================

#' Run KSEA from a DPS table
#' @param dps_file Path to DPS CSV file
#' @param lookup Flanking-to-Gene_Site lookup table
#' @return KSEA result data.frame or NULL
run_ksea <- function(dps_file, lookup) {
    if (!file.exists(dps_file)) return(NULL)

    dps <- fread(dps_file)
    if (!("phosphosite" %in% colnames(dps)) || !("t" %in% colnames(dps))) {
        cat("    [WARN] Missing phosphosite/t columns in:", basename(dps_file), "\n")
        return(NULL)
    }
    
    # Filter and extract flanking
    dps <- dps[is.finite(dps$t) & !is.na(dps$phosphosite), ]
    parts <- strsplit(dps$phosphosite, "|", fixed = TRUE)
    flanking_dps <- vapply(parts, function(x) if(length(x) >= 4) toupper(x[4]) else NA_character_, character(1))
    
    # Map to Gene_Site
    gene_sites <- lookup[flanking_dps]
    
    # Create PX dataframe
    # KSEAapp needs Gene, Residue.Both (e.g. S260), p, FC
    valid_mask <- !is.na(gene_sites)
    if (sum(valid_mask) < 50) {
        cat("    [WARN] Too few mapped phosphosites:", sum(valid_mask), "\n")
        return(NULL)
    }
    
    mapped_gs <- gene_sites[valid_mask]
    mapped_p <- dps$P.Value[valid_mask]
    mapped_t <- dps$t[valid_mask]
    
    # Split Gene_Site into Gene and Residue
    # Format is typically GENE_S260 or GENE_S260-P
    # We want Gene = GENE, Residue.Both = S260
    genes <- sub("_.*$", "", mapped_gs)
    residues <- sub("^[^_]+_", "", mapped_gs)
    residues <- sub("-P$", "", residues) # Remove -P suffix if present
    
    PX <- data.frame(
        Protein = genes,
        Gene = genes,
        Peptide = "",
        Residue.Both = residues,
        p = mapped_p,
        FC = 2^mapped_t, # log2(FC) = t
        stringsAsFactors = FALSE
    )
    
    # Dedup PX to ensure one entry per Gene_Residue (keep max absolute t)
    PX$abs_t <- abs(mapped_t)
    PX <- PX %>% 
        group_by(Gene, Residue.Both) %>% 
        slice_max(abs_t, n = 1, with_ties = FALSE) %>% 
        ungroup() %>% 
        select(-abs_t)
    
    if (nrow(PX) < 50) {
        cat("    [WARN] Too few unique mapped sites:", nrow(PX), "\n")
        return(NULL)
    }

    # Run KSEA.Scores 
    res <- tryCatch({
        KSEA.Scores(KSData, PX, NetworKIN = use_networkin, NetworKIN.cutoff = networkin_cutoff)
    }, error = function(e) {
        cat("    [ERROR] KSEA failed:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(res) || nrow(res) == 0) return(NULL)

    # Sort by FDR (padj)
    res <- res[order(res$FDR), ]
    as.data.frame(res)
}

# ==============================================================================
# Section 4: Main Processing Loop
# ==============================================================================

all_summaries <- list()

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

    summary_row <- list(dataset = ds_folder, cancer_type = ds_cancer)

    for (j in seq_len(nrow(comparisons))) {
        comp_name <- comparisons$name[j]
        comp_label <- comparisons$label[j]

        dps_file <- file.path(ds_deg_dir, paste0("DPS_", comp_name, ".csv"))
        cat("    ", comp_label, ":", sep = "")

        res <- run_ksea(dps_file, flank_lookup)

        if (!is.null(res)) {
            n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
            n_up <- sum(res$FDR < 0.05 & res$z.score > 0, na.rm = TRUE)
            n_down <- sum(res$FDR < 0.05 & res$z.score < 0, na.rm = TRUE)
            cat(" ", nrow(res), "kinases tested,", n_sig, "significant\n")

            # Save CSV
            fwrite(res, file.path(ds_out, paste0("KSEA_", comp_name, ".csv")))

            # Save XLSX
            wb <- createWorkbook()
            addWorksheet(wb, comp_label)
            writeData(wb, 1, res)
            saveWorkbook(wb, file.path(ds_out, paste0("KSEA_", comp_name, ".xlsx")),
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

    all_summaries[[ds_folder]] <- as.data.frame(summary_row, stringsAsFactors = FALSE)
    cat("\n")
}

# ==============================================================================
# Section 5: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating KSEA Summary Statistics\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(all_summaries)
fwrite(summary_df, file.path(output_dir, "KSEA_summary_PhosphoSitePlus.csv"))

wb_sum <- createWorkbook()
addWorksheet(wb_sum, "Summary")
writeData(wb_sum, "Summary", summary_df)
header_style <- createStyle(
    textDecoration = "bold", halign = "center",
    border = "bottom", fgFill = "#4472C4", fontColour = "white"
)
addStyle(wb_sum, "Summary", header_style, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
setColWidths(wb_sum, "Summary", cols = 1:ncol(summary_df), widths = "auto")
saveWorkbook(wb_sum, file.path(output_dir, "KSEA_summary_PhosphoSitePlus.xlsx"), overwrite = TRUE)

sig_cols <- grep("_sig$", colnames(summary_df), value = TRUE)
print_df <- summary_df[, c("dataset", "cancer_type", sig_cols)]
print(as.data.frame(print_df), row.names = FALSE)
cat("\n")

cat("====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

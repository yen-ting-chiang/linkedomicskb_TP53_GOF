################################################################################
# TP53 Phosphoprotein Differential Phosphosite Analysis (limma)
#
# Comparisons:
#   1. TP53mt vs TP53wt
#   2. MUT_GOF vs MUT_LOF
#   3. Hotspots vs MUT_LOF
#   4. MUT_GOF vs TP53wt
#   5. MUT_LOF vs TP53wt
#   6. Hotspots vs TP53wt
#   7. DN vs TP53wt
#   8. Non-DN vs TP53wt
#   9. DN vs non-DN
#
# Covariates: sex, age, tumor purity (WES_purity). TMT plex excluded.
# Protein normalization NOT applied
#
# Methodology adapted from phosphoprotein_DPS_analysis/ pipeline
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(glue)
    library(limma)
    library(vroom)
    library(openxlsx)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

datasets <- data.frame(
    folder = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    cancer_type = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    stringsAsFactors = FALSE
)

output_dir <- file.path(base_path, "phosphoprotein_differential_analysis_without_protein_normalization")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_frac_complete <- 0.50  # phosphosite data is sparser
min_per_group <- 8

cat("====================================================================\n")
cat("TP53 Phosphoprotein Differential Phosphosite Analysis\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Utility Functions
# ==============================================================================

#' Load site-level phosphosite matrix (protein normalization disabled)
load_phosphosite_matrix <- function(ds_id, ds_dir, protein_adjust = FALSE) {
    fp <- file.path(ds_dir, paste0(ds_id, "_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"))
    if (!file.exists(fp)) stop("Cannot find phosphoprotein file: ", fp)
    cat("  Reading phosphosite matrix:", basename(fp), "\n")

    df <- suppressMessages(vroom::vroom(fp, delim = "\t"))
    df <- as.data.frame(df, check.names = FALSE)
    
    stopifnot("idx" %in% colnames(df))
    
    site_id <- as.character(df$idx)
    keep <- !is.na(site_id) & nzchar(site_id)
    df_keep <- df[keep, , drop = FALSE]
    site_id <- site_id[keep]
    
    samp_cols <- setdiff(colnames(df_keep), "idx")
    M_site <- as.matrix(df_keep[, samp_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    
    if (anyDuplicated(site_id)) {
        split_idx <- split(seq_len(nrow(M_site)), site_id)
        M_site_agg <- do.call(rbind, lapply(split_idx, function(idx) {
            if (length(idx) == 1) return(M_site[idx, ])
            apply(M_site[idx, , drop = FALSE], 2, function(x) stats::median(x, na.rm = TRUE))
        }))
        rownames(M_site_agg) <- names(split_idx)
        site_id <- rownames(M_site_agg)
        M_site <- M_site_agg
    } else {
        rownames(M_site) <- site_id
    }
    colnames(M_site) <- samp_cols
    
    # Protein normalization
    if (isTRUE(protein_adjust)) {
        prot_fp <- file.path(ds_dir, paste0(ds_id, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))
        if (!file.exists(prot_fp)) {
            cat("  WARNING: protein file not found, skipping normalization\n")
        } else {
            cat("  Applying protein normalization\n")
            prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
            prot_df <- as.data.frame(prot_df, check.names = FALSE)
            
            if ("idx" %in% colnames(prot_df)) {
                prot_genes <- as.character(prot_df$idx)
                prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), "idx"), drop = FALSE])
                rownames(prot_mat) <- prot_genes
                storage.mode(prot_mat) <- "numeric"
                
                gene_of_site <- sub("\\|.*$", "", site_id)
                
                common <- intersect(colnames(M_site), colnames(prot_mat))
                if (length(common) >= 2) {
                    M_site <- M_site[, common, drop = FALSE]
                    prot_mat <- prot_mat[, common, drop = FALSE]
                    
                    idx <- match(gene_of_site, rownames(prot_mat))
                    hasp <- !is.na(idx)
                    if (any(hasp)) {
                        M_site[hasp, ] <- M_site[hasp, , drop = FALSE] - prot_mat[idx[hasp], , drop = FALSE]
                        cat("    Normalized", sum(hasp), "of", nrow(M_site), "phosphosites\n")
                    }
                }
            }
        }
    }

    cat("  Phosphosite matrix:", nrow(M_site), "sites x", ncol(M_site), "samples\n")
    M_site
}

#' Impute and filter phosphosite matrix (gene-median imputation)
impute_and_filter_phospho <- function(mat, min_frac = 0.50) {
    keep <- rowMeans(is.finite(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (anyNA(m)) {
        meds <- apply(m, 1, function(v) median(v[is.finite(v)], na.rm = TRUE))
        for (i in seq_len(nrow(m))) {
            vi <- m[i, ]
            vi[!is.finite(vi)] <- meds[i]
            m[i, ] <- vi
        }
    }
    m
}

#' Get purity covariate from phenotype file
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    purity_df <- data.frame(
        Sample = sample_ids,
        Purity = rep(NA_real_, length(sample_ids)),
        stringsAsFactors = FALSE
    )
    
    fp <- file.path(ds_dir, paste0(ds_id, "_phenotype.txt"))
    if (file.exists(fp)) {
        ph_df <- suppressMessages(vroom::vroom(fp, delim = "\t"))
        if ("idx" %in% colnames(ph_df) && "WES_purity" %in% colnames(ph_df)) {
            match_idx <- match(sample_ids, ph_df$idx)
            purity_df$Purity <- as.numeric(ph_df$WES_purity[match_idx])
            cat("  Found WES_purity in", basename(fp), "\n")
        }
    }
    
    missing_pct <- sum(is.na(purity_df$Purity)) / nrow(purity_df) * 100
    cat(sprintf("  Purity covariate: %.1f%% missing\n", missing_pct))
    
    purity_df
}

#' Get sex and age covariates from meta file
get_sex_age_covariates <- function(ds_id, ds_dir, sample_ids) {
    cov_df <- data.frame(
        Sample = sample_ids,
        Sex = rep(NA_character_, length(sample_ids)),
        Age = rep(NA_real_, length(sample_ids)),
        stringsAsFactors = FALSE
    )
    
    fp <- file.path(ds_dir, paste0(ds_id, "_meta.txt"))
    if (file.exists(fp)) {
        meta_df <- suppressMessages(vroom::vroom(fp, delim = "\t"))
        if ("case_id" %in% colnames(meta_df)) {
            match_idx <- match(sample_ids, meta_df$case_id)
            if ("Sex" %in% colnames(meta_df)) {
                cov_df$Sex <- as.character(meta_df$Sex[match_idx])
            }
            if ("Age" %in% colnames(meta_df)) {
                cov_df$Age <- as.numeric(meta_df$Age[match_idx])
            }
            cat("  Found Sex/Age in", basename(fp), "\n")
        }
    }
    
    cat(sprintf("  Sex covariate: %.1f%% missing\n", sum(is.na(cov_df$Sex)|cov_df$Sex=="Unknown") / nrow(cov_df) * 100))
    cat(sprintf("  Age covariate: %.1f%% missing\n", sum(is.na(cov_df$Age)) / nrow(cov_df) * 100))
    
    cov_df
}

# ==============================================================================
# Section 3: Differential Expression Worker (limma)
# ==============================================================================

run_limma_group_comparison <- function(M, group_vec, purity_vec, sa_df) {
    common <- intersect(colnames(M), names(group_vec))
    if (length(common) == 0) return(NULL)

    M <- M[, common, drop = FALSE]
    grp <- factor(group_vec[common])
    
    tab <- table(grp)
    if (length(tab) < 2 || any(tab < min_per_group)) {
        cat("    [SKIP] Not enough samples per group (min =", min_per_group, "):", paste(tab, collapse=" vs "), "\n")
        return(NULL)
    }

    DF <- data.frame(
        row.names = common,
        group = grp,
        stringsAsFactors = FALSE
    )
    
    if (!is.null(sa_df)) {
        m_idx <- match(common, sa_df$Sample)
        if ("Sex" %in% names(sa_df)) {
            sex_val <- sa_df$Sex[m_idx]
            if (length(unique(na.omit(sex_val))) > 1) {
                DF$sex <- factor(sex_val)
            }
        }
        if ("Age" %in% names(sa_df)) {
            age_val <- sa_df$Age[m_idx]
            if (sum(!is.na(age_val)) > 3) {
                DF$age <- age_val
            }
        }
    }
    if (!is.null(purity_vec)) {
        p_val <- purity_vec$Purity[match(common, purity_vec$Sample)]
        if (sum(!is.na(p_val)) > 3) {
            DF$purity <- p_val
        }
    }

    keep_samples <- complete.cases(DF)
    if (sum(keep_samples) < 10) {
        cat("    [SKIP] Too few complete cases after covariates:", sum(keep_samples), "\n")
        return(NULL)
    }

    DF <- DF[keep_samples, , drop = FALSE]
    M <- M[, rownames(DF), drop = FALSE]
    
    tab_clean <- table(DF$group)
    if (length(tab_clean) < 2 || any(tab_clean < 2)) {
        cat("    [SKIP] Group dropped after covariate filtering\n")
        return(NULL)
    }

    des <- model.matrix(~ 0 + ., data = DF)
    
    lvls <- levels(DF$group)
    group_cols <- paste0("group", lvls)
    if (!all(group_cols %in% colnames(des))) {
        cat("    [ERROR] Group columns not found in design matrix\n")
        return(NULL)
    }

    # Saturation protection: drop covariates if residual df <= 0
    rnk <- qr(des)$rank
    res_df <- nrow(des) - rnk
    if (res_df <= 0) {
        drop_order <- intersect(c("sex", "age", "purity"), colnames(DF))
        for (cv in drop_order) {
            DF[[cv]] <- NULL
            cat("    [saturation] dropping covariate:", cv, "\n")
            des <- model.matrix(~ 0 + ., data = DF)
            rnk <- qr(des)$rank
            res_df <- nrow(des) - rnk
            if (res_df > 0) break
        }
        if (res_df <= 0) {
            cat("    [SKIP] No residual df after removing covariates\n")
            return(NULL)
        }
        lvls <- levels(DF$group)
        group_cols <- paste0("group", lvls)
        M <- M[, rownames(DF), drop = FALSE]
    }

    contrast_str <- paste0(group_cols[2], " - ", group_cols[1])
    contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = des)

    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= (rnk + 1L)
    if (sum(keep_rows) == 0) {
        cat("    [SKIP] No phosphosites pass df filter\n")
        return(NULL)
    }
    M <- M[keep_rows, , drop = FALSE]

    fit <- limma::lmFit(M, design = des)
    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    eb <- limma::eBayes(fit2, trend = TRUE)

    tbl <- limma::topTable(eb, number = Inf, sort.by = "P")
    tbl$phosphosite <- rownames(tbl)
    tbl <- tbl[, c("phosphosite", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    cat("    Results:", nrow(tbl), "phosphosites tested,",
        sum(tbl$adj.P.Val < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")

    tbl
}

#' Add gene_symbol_phosphosite column by mapping ENSEMBL IDs to HGNC gene symbols
#'
#' Parses phosphosite ID format: "ENSG...|ENSP...|SITE|SEQ|NUM"
#' Extracts ENSEMBL gene ID (1st field) and phosphosite (3rd field)
#' Produces "SYMBOL_SITE" format (e.g., "PDZD4_T150")
#'
#' @param deg_df Data frame with a "phosphosite" column
#' @param mapping Named character vector (names = versioned ENSEMBL gene IDs, values = gene symbols)
#' @return Data frame with gene_symbol_phosphosite column inserted before phosphosite column
add_gene_symbol_phosphosite <- function(deg_df, mapping) {
    # Parse phosphosite IDs: extract ENSEMBL gene ID and site
    parts <- strsplit(deg_df$phosphosite, "\\|")
    ensembl_gene <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
    site <- vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))

    # Map ENSEMBL gene ID to gene symbol
    symbol <- mapping[ensembl_gene]

    # Build gene_symbol_phosphosite: SYMBOL_SITE
    deg_df$gene_symbol_phosphosite <- ifelse(
        !is.na(symbol) & !is.na(site),
        paste0(symbol, "_", site),
        ""
    )

    # Reorder: gene_symbol_phosphosite first, then phosphosite, then the rest
    deg_df <- deg_df[, c("gene_symbol_phosphosite", "phosphosite",
                         setdiff(names(deg_df), c("gene_symbol_phosphosite", "phosphosite")))]
    deg_df
}

# ==============================================================================
# Section 4: Helper to run one comparison and save results
# ==============================================================================

run_and_save_comparison <- function(comparison_name, file_prefix, sheet_name,
                                    group_vec, mat, purity_vec, sa_df,
                                    ds_out, ds_summary, ensembl2symbol) {
    cat("\n  ---", comparison_name, "---\n")
    deg <- run_limma_group_comparison(mat, group_vec, purity_vec, sa_df)
    prefix_clean <- gsub("[^A-Za-z0-9_]", "_", file_prefix)

    if (!is.null(deg)) {
        deg <- add_gene_symbol_phosphosite(deg, ensembl2symbol)
        fwrite(deg, file.path(ds_out, paste0("DPS_", prefix_clean, ".csv")))
        wb <- createWorkbook()
        addWorksheet(wb, sheet_name)
        writeData(wb, 1, deg)
        saveWorkbook(wb, file.path(ds_out, paste0("DPS_", prefix_clean, ".xlsx")), overwrite = TRUE)
        ds_summary[[paste0(prefix_clean, "_total")]] <- nrow(deg)
        ds_summary[[paste0(prefix_clean, "_sig")]] <- sum(deg$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary[[paste0(prefix_clean, "_up")]] <- sum(deg$adj.P.Val < 0.05 & deg$logFC > 0, na.rm = TRUE)
        ds_summary[[paste0(prefix_clean, "_down")]] <- sum(deg$adj.P.Val < 0.05 & deg$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary[[paste0(prefix_clean, "_total")]] <- NA
        ds_summary[[paste0(prefix_clean, "_sig")]] <- NA
        ds_summary[[paste0(prefix_clean, "_up")]] <- NA
        ds_summary[[paste0(prefix_clean, "_down")]] <- NA
    }
    ds_summary
}

# ==============================================================================
# Section 5: Main Processing Loop
# ==============================================================================

tp53_file <- file.path(base_path, "TP53_mutation_classification",
                       "all_CPTAC_TP53_classification.csv")
tp53_all <- read.csv(tp53_file, stringsAsFactors = FALSE)
cat("Loaded TP53 classification:", nrow(tp53_all), "samples\n\n")

all_deg_summary <- list()

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]
    ds_dir <- file.path(base_path, ds_folder)

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    # --- Load phosphosite matrix (protein normalization disabled) ---
    mat0 <- tryCatch(
        load_phosphosite_matrix(ds_folder, ds_dir, protein_adjust = FALSE),
        error = function(e) {
            cat("  ERROR loading phosphosite matrix:", conditionMessage(e), "\n")
            NULL
        }
    )
    if (is.null(mat0)) { cat("  Skipping dataset\n\n"); next }

    # Impute and filter
    mat <- impute_and_filter_phospho(mat0, min_frac = min_frac_complete)
    cat("  After filtering:", nrow(mat), "phosphosites x", ncol(mat), "samples\n")

    # --- Build ENSEMBL-to-gene-symbol lookup for this dataset ---
    ensembl_ids_versioned <- unique(sub("\\|.*$", "", rownames(mat)))
    ensembl_ids_bare <- sub("\\.\\d+$", "", ensembl_ids_versioned)
    symbol_map <- tryCatch(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids_bare,
                              column = "SYMBOL", keytype = "ENSEMBL",
                              multiVals = "first"),
        error = function(e) { cat("  [WARN] Gene symbol mapping failed:", e$message, "\n"); setNames(rep(NA_character_, length(ensembl_ids_bare)), ensembl_ids_bare) }
    )
    # Create mapping keyed by versioned ENSEMBL gene IDs
    ensembl2symbol <- setNames(as.character(symbol_map[ensembl_ids_bare]), ensembl_ids_versioned)
    cat("  Gene symbol mapping:", sum(!is.na(ensembl2symbol)), "of", length(ensembl2symbol), "genes mapped\n")

    # --- Get covariates ---
    purity_vec <- get_purity_covariate(ds_folder, ds_dir, colnames(mat))
    sa_df <- get_sex_age_covariates(ds_folder, ds_dir, colnames(mat))
    cat("  Covariates loaded (purity, sex, age)\n")

    # --- Get TP53 classification ---
    tp53_ds <- tp53_all %>% filter(cancer_type == ds_cancer)
    tp53_mut <- tp53_ds %>% filter(mt == 1)

    wt_samples <- tp53_ds$sample_id[tp53_ds$wt == 1]
    gof_samples <- tp53_mut$sample_id[tp53_mut$GOF == 1]
    lof_samples <- tp53_mut$sample_id[tp53_mut$LOF == 1]
    hotspot_samples <- tp53_mut$sample_id[tp53_mut$hotspot == 1]
    dn_samples <- tp53_mut$sample_id[tp53_mut$DN == 1]
    nondn_samples <- tp53_mut$sample_id[tp53_mut$non_DN == 1]

    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)
    ds_summary <- list(dataset = ds_folder, cancer_type = ds_cancer)

    # === 9 Comparisons ===

    # 1. TP53mt vs TP53wt
    g1 <- setNames(ifelse(tp53_ds$mt == 1, "TP53mt", "TP53wt"), tp53_ds$sample_id)
    g1 <- factor(g1[!is.na(g1)], levels = c("TP53wt", "TP53mt"))
    cat("    TP53wt:", sum(g1 == "TP53wt"), " TP53mt:", sum(g1 == "TP53mt"), "\n")
    ds_summary <- run_and_save_comparison("Comparison 1: TP53mt vs TP53wt", "TP53mt_vs_TP53wt", "mt_vs_wt",
        g1, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 2. MUT_GOF vs MUT_LOF
    g2 <- setNames(c(rep("MUT_GOF", length(gof_samples)), rep("MUT_LOF", length(lof_samples))), c(gof_samples, lof_samples))
    g2 <- factor(g2, levels = c("MUT_LOF", "MUT_GOF"))
    ds_summary <- run_and_save_comparison("Comparison 2: GOF vs LOF", "MUT_GOF_vs_MUT_LOF", "GOF_vs_LOF",
        g2, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 3. Hotspot vs MUT_LOF
    g3 <- setNames(c(rep("Hotspot", length(hotspot_samples)), rep("MUT_LOF", length(lof_samples))), c(hotspot_samples, lof_samples))
    g3 <- g3[!duplicated(names(g3))]; g3 <- factor(g3, levels = c("MUT_LOF", "Hotspot"))
    ds_summary <- run_and_save_comparison("Comparison 3: Hotspot vs LOF", "Hotspot_vs_MUT_LOF", "Hot_vs_LOF",
        g3, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 4. MUT_GOF vs TP53wt
    g4 <- setNames(c(rep("MUT_GOF", length(gof_samples)), rep("TP53wt", length(wt_samples))), c(gof_samples, wt_samples))
    g4 <- factor(g4, levels = c("TP53wt", "MUT_GOF"))
    ds_summary <- run_and_save_comparison("Comparison 4: GOF vs TP53wt", "MUT_GOF_vs_TP53wt", "GOF_vs_wt",
        g4, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 5. MUT_LOF vs TP53wt
    g5 <- setNames(c(rep("MUT_LOF", length(lof_samples)), rep("TP53wt", length(wt_samples))), c(lof_samples, wt_samples))
    g5 <- factor(g5, levels = c("TP53wt", "MUT_LOF"))
    ds_summary <- run_and_save_comparison("Comparison 5: LOF vs TP53wt", "MUT_LOF_vs_TP53wt", "LOF_vs_wt",
        g5, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 6. Hotspot vs TP53wt
    g6 <- setNames(c(rep("Hotspot", length(hotspot_samples)), rep("TP53wt", length(wt_samples))), c(hotspot_samples, wt_samples))
    g6 <- factor(g6, levels = c("TP53wt", "Hotspot"))
    ds_summary <- run_and_save_comparison("Comparison 6: Hotspot vs TP53wt", "Hotspot_vs_TP53wt", "Hot_vs_wt",
        g6, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 7. DN vs TP53wt
    g7 <- setNames(c(rep("DN", length(dn_samples)), rep("TP53wt", length(wt_samples))), c(dn_samples, wt_samples))
    g7 <- factor(g7, levels = c("TP53wt", "DN"))
    ds_summary <- run_and_save_comparison("Comparison 7: DN vs TP53wt", "DN_vs_TP53wt", "DN_vs_wt",
        g7, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 8. Non-DN vs TP53wt
    g8 <- setNames(c(rep("NonDN", length(nondn_samples)), rep("TP53wt", length(wt_samples))), c(nondn_samples, wt_samples))
    g8 <- factor(g8, levels = c("TP53wt", "NonDN"))
    ds_summary <- run_and_save_comparison("Comparison 8: NonDN vs TP53wt", "NonDN_vs_TP53wt", "NonDN_vs_wt",
        g8, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    # 9. DN vs non-DN
    g9 <- setNames(c(rep("DN", length(dn_samples)), rep("NonDN", length(nondn_samples))), c(dn_samples, nondn_samples))
    g9 <- g9[!duplicated(names(g9))]; g9 <- factor(g9, levels = c("NonDN", "DN"))
    ds_summary <- run_and_save_comparison("Comparison 9: DN vs NonDN", "DN_vs_NonDN", "DN_vs_NonDN",
        g9, mat, purity_vec, sa_df, ds_out, ds_summary, ensembl2symbol)

    all_deg_summary[[ds_folder]] <- as.data.frame(ds_summary, stringsAsFactors = FALSE)
    cat("\n")
}

# ==============================================================================
# Section 6: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Summary Statistics\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(all_deg_summary)
fwrite(summary_df, file.path(output_dir, "DPS_summary_statistics.csv"))

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
saveWorkbook(wb_sum, file.path(output_dir, "DPS_summary_statistics.xlsx"), overwrite = TRUE)

cat("\n=== Differential Phosphosite Summary (padj < 0.05) ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

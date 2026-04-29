################################################################################
# TP53 Protein-Level Differential Expression Analysis (limma)
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
# Covariates: TMT plex (batch), sex, age, tumor purity
#
# Reference: Shahbandi et al. (2023) Cell Death Discovery
# Methodology adapted from protein_level_DEG_analysis/ pipeline
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
    library(openxlsx)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")
datasets <- data.frame(
    folder = cancer_types,
    cancer_type = cancer_types,
    stringsAsFactors = FALSE
)

# Output directory
output_dir <- file.path(base_path, "protein_differential_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum fraction of non-NA values per gene
min_frac_complete <- 0.75
# Minimum samples per group
min_per_group <- 8

cat("====================================================================\n")
cat("TP53 Protein Differential Expression Analysis\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Utility Functions (adapted from protein_level_DEG_analysis)
# ==============================================================================

#' Load protein quantification matrix
load_protein_matrix <- function(ds_dir, ds_cancer) {
    fp <- file.path(ds_dir, paste0(ds_cancer, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))
    if (!file.exists(fp)) stop(paste("File not found:", fp))

    cat("  Reading protein matrix:", basename(fp), "\n")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))

    # Identify gene column
    gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol)
    dat$Gene <- sub("\\|.*$", "", dat$Gene)

    # Identify sample columns
    not_sample <- c("Gene")
    sample_cols <- setdiff(names(dat), not_sample)

    # Create matrix
    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")

    rn <- m$Gene
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    # Handle duplicate genes by averaging
    if (anyDuplicated(rownames(m))) {
        cat("  Averaging duplicate genes\n")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) /
            as.vector(table(rownames(m)))
    }

    cat("  Matrix dimensions:", nrow(m), "genes x", ncol(m), "samples\n")
    m
}

#' Filter low-coverage genes (imputation removed)
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    m
}

#' Get tumor purity covariate
get_purity_covariate <- function(ds_dir, ds_cancer, sample_ids) {
    fp <- file.path(ds_dir, paste0(ds_cancer, "_phenotype.txt"))
    if (!file.exists(fp)) return(NULL)

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(clin)) return(NULL)

    purity_col <- "WES_purity"
    if (!(purity_col %in% colnames(clin))) return(NULL)

    id_col <- colnames(clin)[1]

    clin_df <- data.frame(
        sample_id = as.character(clin[[id_col]]),
        purity = suppressWarnings(as.numeric(clin[[purity_col]])),
        stringsAsFactors = FALSE
    )

    purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        matches <- clin_df$sample_id[clin_df$sample_id == sid |
                                         sapply(clin_df$sample_id, function(pid) grepl(pid, sid, fixed = TRUE))]
        if (length(matches) > 0) {
            purity_vec[sid] <- clin_df$purity[clin_df$sample_id == matches[1]]
        }
    }
    purity_vec
}

#' Get sex and age covariates
get_sex_age_covariates <- function(ds_dir, ds_cancer, sample_ids) {
    fp <- file.path(ds_dir, paste0(ds_cancer, "_meta.txt"))
    empty_df <- data.frame(
        sex = rep(NA_character_, length(sample_ids)),
        age = rep(NA_real_, length(sample_ids)),
        row.names = sample_ids, stringsAsFactors = FALSE
    )
    if (!file.exists(fp)) return(empty_df)

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(clin)) return(empty_df)

    id_col <- colnames(clin)[1]

    sex_values <- if ("Sex" %in% colnames(clin)) as.character(clin[["Sex"]]) else rep(NA_character_, nrow(clin))
    age_values <- if ("Age" %in% colnames(clin)) suppressWarnings(as.numeric(clin[["Age"]])) else rep(NA_real_, nrow(clin))

    clin_aligned <- data.frame(
        SAMPLE_ID = as.character(clin[[id_col]]), 
        sex = sex_values, 
        age = age_values,
        stringsAsFactors = FALSE
    )
    
    sex_res <- rep(NA_character_, length(sample_ids))
    age_res <- rep(NA_real_, length(sample_ids))
    
    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        matches <- clin_aligned$SAMPLE_ID[clin_aligned$SAMPLE_ID == sid |
                                          sapply(clin_aligned$SAMPLE_ID, function(pid) grepl(pid, sid, fixed = TRUE))]
        if (length(matches) > 0) {
            sex_res[i] <- clin_aligned$sex[clin_aligned$SAMPLE_ID == matches[1]]
            age_res[i] <- clin_aligned$age[clin_aligned$SAMPLE_ID == matches[1]]
        }
    }
    
    data.frame(sex = sex_res, age = age_res, row.names = sample_ids, stringsAsFactors = FALSE)
}

#' Safely coerce covariates (remove single-level factors, etc.)
coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]
        if (is.factor(v) || is.character(v) || is.logical(v)) {
            v <- factor(v)
            lv <- levels(droplevels(v[!is.na(v)]))
            if (length(lv) <= 1) {
                keep[cn] <- FALSE
                cat("    [covars] drop single-level factor:", cn, "\n")
            } else {
                df[[cn]] <- v
            }
        } else {
            df[[cn]] <- suppressWarnings(as.numeric(v))
        }
    }
    df[, keep, drop = FALSE]
}

# ==============================================================================
# Section 3: Limma DEG Function for Group Comparison
# ==============================================================================

#' Run limma DEG analysis for a binary group comparison
#'
#' @param mat Expression matrix (genes x samples), already imputed
#' @param group_vec Named factor with two levels (test vs reference)
#' @param batch_fac Batch factor (optional)
#' @param purity_vec Purity numeric vector (optional)
#' @param sa_df Sex/age data frame (optional)
#' @return Data frame with DEG results, or NULL if insufficient data
run_limma_group_comparison <- function(mat, group_vec,
                                       purity_vec = NULL,
                                       sa_df = NULL) {
    # Align samples
    common <- intersect(colnames(mat), names(group_vec))
    group_vec <- group_vec[common]
    group_vec <- droplevels(group_vec)

    if (nlevels(group_vec) < 2) {
        cat("    [SKIP] Less than 2 group levels\n")
        return(NULL)
    }

    tab <- table(group_vec)
    if (any(tab < min_per_group)) {
        cat("    [SKIP] Group too small:", paste(names(tab), "=", tab, collapse = ", "), "\n")
        return(NULL)
    }

    M <- mat[, common, drop = FALSE]
    so <- common  # sample order

    # Build covariate data frame
    DF <- data.frame(group = group_vec, row.names = so)



    # Add purity
    if (!is.null(purity_vec)) {
        p <- suppressWarnings(as.numeric(purity_vec[so]))
        if (sum(is.finite(p)) > 0) {
            DF$purity <- p
        }
    }

    # Add sex and age
    if (!is.null(sa_df)) {
        if ("sex" %in% colnames(sa_df)) {
            s <- sa_df[so, "sex"]
            if (sum(!is.na(s)) > 0) {
                DF$sex <- factor(s)
            }
        }
        if ("age" %in% colnames(sa_df)) {
            a <- suppressWarnings(as.numeric(sa_df[so, "age"]))
            if (sum(is.finite(a)) > 0) {
                DF$age <- a
            }
        }
    }

    # Remove all-NA columns
    all_na <- vapply(DF, function(z) all(is.na(z)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]

    # Clean covariates
    DF <- coerce_covariates_safely(DF)

    # Ensure group is still present
    if (!("group" %in% colnames(DF))) {
        cat("    [ERROR] Group variable dropped during covariate cleaning\n")
        return(NULL)
    }

    # Keep only complete cases for all variables to ensure M and des match
    DF <- na.omit(DF)
    
    if (nlevels(droplevels(DF$group)) < 2) {
        cat("    [SKIP] Less than 2 group levels remaining after dropping missing covariates\n")
        return(NULL)
    }
    
    M <- M[, rownames(DF), drop = FALSE]

    # Build design matrix
    des <- model.matrix(~ 0 + ., data = DF)

    # Make contrast: test level - reference level (second level vs first)
    lvls <- levels(DF$group)
    # Column names in design matrix
    group_cols <- paste0("group", lvls)

    # Check both group columns exist
    if (!all(group_cols %in% colnames(des))) {
        cat("    [ERROR] Group columns not found in design matrix\n")
        return(NULL)
    }

    contrast_str <- paste0(group_cols[2], " - ", group_cols[1])
    contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = des)

    # Filter genes with insufficient df
    rnk <- qr(des)$rank
    need <- rnk + 1L
    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= need
    if (sum(keep_rows) == 0) {
        cat("    [SKIP] No genes pass df filter\n")
        return(NULL)
    }
    M <- M[keep_rows, , drop = FALSE]

    # Run limma
    fit <- limma::lmFit(M, design = des)
    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    eb <- limma::eBayes(fit2, trend = TRUE)

    tbl <- limma::topTable(eb, number = Inf, sort.by = "P")
    tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    cat("    Results:", nrow(tbl), "genes tested,",
        sum(tbl$adj.P.Val < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")

    tbl
}

#' Add gene_symbol column by mapping ENSEMBL IDs to HGNC gene symbols
#'
#' @param deg_df Data frame with a "gene" column containing ENSEMBL IDs
#' @param mapping Named character vector (names = ENSEMBL IDs with version, values = gene symbols)
#' @return Data frame with gene_symbol column inserted before gene column
add_gene_symbol <- function(deg_df, mapping) {
    deg_df$gene_symbol <- mapping[deg_df$gene]
    deg_df$gene_symbol[is.na(deg_df$gene_symbol)] <- ""
    # Reorder: gene_symbol first, then gene, then the rest
    deg_df <- deg_df[, c("gene_symbol", "gene", setdiff(names(deg_df), c("gene_symbol", "gene")))]
    deg_df
}

# ==============================================================================
# Section 4: Main Processing Loop
# ==============================================================================

# Load TP53 classification results
tp53_file <- file.path(base_path, "TP53_mutation_classification",
                       "all_CPTAC_TP53_classification.csv")
tp53_all <- read.csv(tp53_file, stringsAsFactors = FALSE)
cat("Loaded TP53 classification:", nrow(tp53_all), "samples\n\n")

# Storage for summary
all_deg_summary <- list()
imputation_report_all <- list()

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]
    ds_dir <- file.path(base_path, ds_folder)

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    # --- Load protein matrix ---
    mat0 <- load_protein_matrix(ds_dir, ds_cancer)

    # Log-transform if needed (values > 100 suggest raw counts/intensities)
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) {
        cat("  Log2-transforming data\n")
        mat0 <- log2(mat0 + 1)
    }

    # Filter
    mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
    cat("  After filtering:", nrow(mat), "genes x", ncol(mat), "samples\n")

    # --- Build ENSEMBL-to-gene-symbol lookup for this dataset ---
    ensembl_ids_versioned <- rownames(mat)
    ensembl_ids_bare <- sub("\\.\\d+$", "", ensembl_ids_versioned)
    symbol_map <- tryCatch(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids_bare,
                              column = "SYMBOL", keytype = "ENSEMBL",
                              multiVals = "first"),
        error = function(e) { cat("  [WARN] Gene symbol mapping failed:", e$message, "\n"); setNames(rep(NA_character_, length(ensembl_ids_bare)), ensembl_ids_bare) }
    )
    # Create mapping keyed by versioned ENSEMBL IDs (to match deg$gene)
    ensembl2symbol <- setNames(as.character(symbol_map[ensembl_ids_bare]), ensembl_ids_versioned)
    cat("  Gene symbol mapping:", sum(!is.na(ensembl2symbol)), "of", length(ensembl2symbol), "genes mapped\n")

    # --- Get covariates ---
    purity_vec <- get_purity_covariate(ds_dir, ds_cancer, colnames(mat))
    sa_df <- get_sex_age_covariates(ds_dir, ds_cancer, colnames(mat))
    
    # --- Covariate Imputation at Dataset Level ---
    n_samples <- ncol(mat)
    
    # Purity
    p_miss_cnt <- 0; p_imp_val <- NA
    if (!is.null(purity_vec)) {
        p <- suppressWarnings(as.numeric(purity_vec))
        p_miss_cnt <- sum(is.na(p))
        if (sum(is.finite(p)) >= n_samples * 0.6) {
            p_imp_val <- median(p, na.rm = TRUE)
            p[is.na(p)] <- p_imp_val
            purity_vec <- setNames(p, names(purity_vec))
        } else {
            purity_vec <- NULL # Too much missing data, discard covariate
        }
    }
    
    # Age
    a_miss_cnt <- 0; a_imp_val <- NA
    if (!is.null(sa_df) && "age" %in% colnames(sa_df)) {
        a <- suppressWarnings(as.numeric(sa_df$age))
        a_miss_cnt <- sum(is.na(a))
        if (sum(is.finite(a)) >= n_samples * 0.8) {
            a_imp_val <- median(a, na.rm = TRUE)
            a[is.na(a)] <- a_imp_val
            sa_df$age <- a
        } else {
            sa_df$age <- NULL
        }
    }
    
    # Sex
    s_miss_cnt <- 0; s_imp_val <- NA
    if (!is.null(sa_df) && "sex" %in% colnames(sa_df)) {
        s <- sa_df$sex
        s_miss_cnt <- sum(is.na(s))
        if (sum(!is.na(s)) >= n_samples * 0.8) {
            s_imp_val <- names(sort(table(s), decreasing = TRUE))[1]
            s[is.na(s)] <- s_imp_val
            sa_df$sex <- s
        } else {
            sa_df$sex <- NULL
        }
    }
    
    imputation_report_row <- data.frame(
        dataset = ds_cancer,
        total_samples = n_samples,
        purity_missing_rate = sprintf("%.1f%%", p_miss_cnt / n_samples * 100),
        purity_imputed_value = ifelse(is.na(p_imp_val), "Dropped (>40% NA)", round(p_imp_val, 4)),
        age_missing_rate = sprintf("%.1f%%", a_miss_cnt / n_samples * 100),
        age_imputed_value = ifelse(is.na(a_imp_val), "Dropped (>20% NA)", round(a_imp_val, 1)),
        sex_missing_rate = sprintf("%.1f%%", s_miss_cnt / n_samples * 100),
        sex_imputed_value = ifelse(is.na(s_imp_val), "Dropped (>20% NA)", s_imp_val),
        stringsAsFactors = FALSE
    )
    imputation_report_all[[ds_cancer]] <- imputation_report_row
    
    cat("\n  [Covariate Imputation Report]\n")
    cat(sprintf("    Purity: Missing %s | Imputed value: %s\n", imputation_report_row$purity_missing_rate, imputation_report_row$purity_imputed_value))
    cat(sprintf("    Age   : Missing %s | Imputed value: %s\n", imputation_report_row$age_missing_rate, imputation_report_row$age_imputed_value))
    cat(sprintf("    Sex   : Missing %s | Imputed value: %s\n", imputation_report_row$sex_missing_rate, imputation_report_row$sex_imputed_value))
    cat("  --------------------------------------\n")

    # --- Get TP53 classification for this dataset ---
    tp53_ds <- tp53_all %>% filter(cancer_type == ds_cancer)

    # Create output subdirectory
    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    # Track summary for this dataset
    ds_summary <- data.frame(
        dataset = ds_folder,
        cancer_type = ds_cancer,
        stringsAsFactors = FALSE
    )

    # ========================================================================
    # Comparison 1: TP53mt vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 1: TP53mt vs TP53wt ---\n")

    group1 <- setNames(
        ifelse(tp53_ds$mt == 1, "TP53mt", "TP53wt"),
        tp53_ds$sample_id
    )
    group1 <- factor(group1[!is.na(group1)], levels = c("TP53wt", "TP53mt"))

    cat("    Group sizes: TP53wt=", sum(group1 == "TP53wt", na.rm = TRUE),
        " TP53mt=", sum(group1 == "TP53mt", na.rm = TRUE), "\n")

    deg1 <- run_limma_group_comparison(mat, group1, purity_vec, sa_df)

    if (!is.null(deg1)) {
        deg1 <- add_gene_symbol(deg1, ensembl2symbol)
        fwrite(deg1, file.path(ds_out, "DEG_TP53mt_vs_TP53wt.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "TP53mt_vs_TP53wt")
        writeData(wb, 1, deg1)
        saveWorkbook(wb, file.path(ds_out, "DEG_TP53mt_vs_TP53wt.xlsx"), overwrite = TRUE)

        ds_summary$TP53mt_vs_TP53wt_total <- nrow(deg1)
        ds_summary$TP53mt_vs_TP53wt_sig <- sum(deg1$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$TP53mt_vs_TP53wt_up <- sum(deg1$adj.P.Val < 0.05 & deg1$logFC > 0, na.rm = TRUE)
        ds_summary$TP53mt_vs_TP53wt_down <- sum(deg1$adj.P.Val < 0.05 & deg1$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$TP53mt_vs_TP53wt_total <- NA
        ds_summary$TP53mt_vs_TP53wt_sig <- NA
        ds_summary$TP53mt_vs_TP53wt_up <- NA
        ds_summary$TP53mt_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 2: MUT_GOF vs MUT_LOF
    # ========================================================================
    cat("\n  --- Comparison 2: MUT_GOF vs MUT_LOF ---\n")

    tp53_mut <- tp53_ds %>% filter(mt == 1)
    gof_samples <- tp53_mut$sample_id[tp53_mut$GOF == 1]
    lof_samples <- tp53_mut$sample_id[tp53_mut$LOF == 1]

    cat("    GOF samples:", length(gof_samples), " LOF samples:", length(lof_samples), "\n")

    group2_ids <- c(gof_samples, lof_samples)
    group2 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("MUT_LOF", length(lof_samples))),
        group2_ids
    )
    group2 <- factor(group2, levels = c("MUT_LOF", "MUT_GOF"))

    deg2 <- run_limma_group_comparison(mat, group2, purity_vec, sa_df)

    if (!is.null(deg2)) {
        deg2 <- add_gene_symbol(deg2, ensembl2symbol)
        fwrite(deg2, file.path(ds_out, "DEG_MUT_GOF_vs_MUT_LOF.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "GOF_vs_LOF")
        writeData(wb, 1, deg2)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_GOF_vs_MUT_LOF.xlsx"), overwrite = TRUE)

        ds_summary$GOF_vs_LOF_total <- nrow(deg2)
        ds_summary$GOF_vs_LOF_sig <- sum(deg2$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$GOF_vs_LOF_up <- sum(deg2$adj.P.Val < 0.05 & deg2$logFC > 0, na.rm = TRUE)
        ds_summary$GOF_vs_LOF_down <- sum(deg2$adj.P.Val < 0.05 & deg2$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$GOF_vs_LOF_total <- NA
        ds_summary$GOF_vs_LOF_sig <- NA
        ds_summary$GOF_vs_LOF_up <- NA
        ds_summary$GOF_vs_LOF_down <- NA
    }

    # ========================================================================
    # Comparison 3: Hotspots vs MUT_LOF
    # ========================================================================
    cat("\n  --- Comparison 3: Hotspots vs MUT_LOF ---\n")

    hotspot_samples <- tp53_mut$sample_id[tp53_mut$hotspot == 1]

    cat("    Hotspot samples:", length(hotspot_samples), " LOF samples:", length(lof_samples), "\n")

    group3_ids <- c(hotspot_samples, lof_samples)
    # Remove duplicates (a hotspot might also be LOF classified)
    group3_labels <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("MUT_LOF", length(lof_samples))),
        group3_ids
    )
    # If a sample appears in both, keep Hotspot classification
    group3_labels <- group3_labels[!duplicated(names(group3_labels))]
    group3 <- factor(group3_labels, levels = c("MUT_LOF", "Hotspot"))

    deg3 <- run_limma_group_comparison(mat, group3, purity_vec, sa_df)

    if (!is.null(deg3)) {
        deg3 <- add_gene_symbol(deg3, ensembl2symbol)
        fwrite(deg3, file.path(ds_out, "DEG_Hotspot_vs_MUT_LOF.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "Hotspot_vs_LOF")
        writeData(wb, 1, deg3)
        saveWorkbook(wb, file.path(ds_out, "DEG_Hotspot_vs_MUT_LOF.xlsx"), overwrite = TRUE)

        ds_summary$Hotspot_vs_LOF_total <- nrow(deg3)
        ds_summary$Hotspot_vs_LOF_sig <- sum(deg3$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$Hotspot_vs_LOF_up <- sum(deg3$adj.P.Val < 0.05 & deg3$logFC > 0, na.rm = TRUE)
        ds_summary$Hotspot_vs_LOF_down <- sum(deg3$adj.P.Val < 0.05 & deg3$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$Hotspot_vs_LOF_total <- NA
        ds_summary$Hotspot_vs_LOF_sig <- NA
        ds_summary$Hotspot_vs_LOF_up <- NA
        ds_summary$Hotspot_vs_LOF_down <- NA
    }

    # ========================================================================
    # Comparison 4: MUT_GOF vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 4: MUT_GOF vs TP53wt ---\n")

    wt_samples <- tp53_ds$sample_id[tp53_ds$wt == 1]
    cat("    GOF samples:", length(gof_samples), " TP53wt samples:", length(wt_samples), "\n")

    group4_ids <- c(gof_samples, wt_samples)
    group4 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("TP53wt", length(wt_samples))),
        group4_ids
    )
    group4 <- factor(group4, levels = c("TP53wt", "MUT_GOF"))

    deg4 <- run_limma_group_comparison(mat, group4, purity_vec, sa_df)

    if (!is.null(deg4)) {
        deg4 <- add_gene_symbol(deg4, ensembl2symbol)
        fwrite(deg4, file.path(ds_out, "DEG_MUT_GOF_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "GOF_vs_TP53wt")
        writeData(wb, 1, deg4)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_GOF_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$GOF_vs_TP53wt_total <- nrow(deg4)
        ds_summary$GOF_vs_TP53wt_sig <- sum(deg4$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$GOF_vs_TP53wt_up <- sum(deg4$adj.P.Val < 0.05 & deg4$logFC > 0, na.rm = TRUE)
        ds_summary$GOF_vs_TP53wt_down <- sum(deg4$adj.P.Val < 0.05 & deg4$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$GOF_vs_TP53wt_total <- NA
        ds_summary$GOF_vs_TP53wt_sig <- NA
        ds_summary$GOF_vs_TP53wt_up <- NA
        ds_summary$GOF_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 5: MUT_LOF vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 5: MUT_LOF vs TP53wt ---\n")

    cat("    LOF samples:", length(lof_samples), " TP53wt samples:", length(wt_samples), "\n")

    group5_ids <- c(lof_samples, wt_samples)
    group5 <- setNames(
        c(rep("MUT_LOF", length(lof_samples)), rep("TP53wt", length(wt_samples))),
        group5_ids
    )
    group5 <- factor(group5, levels = c("TP53wt", "MUT_LOF"))

    deg5 <- run_limma_group_comparison(mat, group5, purity_vec, sa_df)

    if (!is.null(deg5)) {
        deg5 <- add_gene_symbol(deg5, ensembl2symbol)
        fwrite(deg5, file.path(ds_out, "DEG_MUT_LOF_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "LOF_vs_TP53wt")
        writeData(wb, 1, deg5)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_LOF_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$LOF_vs_TP53wt_total <- nrow(deg5)
        ds_summary$LOF_vs_TP53wt_sig <- sum(deg5$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$LOF_vs_TP53wt_up <- sum(deg5$adj.P.Val < 0.05 & deg5$logFC > 0, na.rm = TRUE)
        ds_summary$LOF_vs_TP53wt_down <- sum(deg5$adj.P.Val < 0.05 & deg5$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$LOF_vs_TP53wt_total <- NA
        ds_summary$LOF_vs_TP53wt_sig <- NA
        ds_summary$LOF_vs_TP53wt_up <- NA
        ds_summary$LOF_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 6: Hotspots vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 6: Hotspots vs TP53wt ---\n")

    cat("    Hotspot samples:", length(hotspot_samples), " TP53wt samples:", length(wt_samples), "\n")

    group6_ids <- c(hotspot_samples, wt_samples)
    group6 <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("TP53wt", length(wt_samples))),
        group6_ids
    )
    group6 <- factor(group6, levels = c("TP53wt", "Hotspot"))

    deg6 <- run_limma_group_comparison(mat, group6, purity_vec, sa_df)

    if (!is.null(deg6)) {
        deg6 <- add_gene_symbol(deg6, ensembl2symbol)
        fwrite(deg6, file.path(ds_out, "DEG_Hotspot_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "Hotspot_vs_TP53wt")
        writeData(wb, 1, deg6)
        saveWorkbook(wb, file.path(ds_out, "DEG_Hotspot_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$Hotspot_vs_TP53wt_total <- nrow(deg6)
        ds_summary$Hotspot_vs_TP53wt_sig <- sum(deg6$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$Hotspot_vs_TP53wt_up <- sum(deg6$adj.P.Val < 0.05 & deg6$logFC > 0, na.rm = TRUE)
        ds_summary$Hotspot_vs_TP53wt_down <- sum(deg6$adj.P.Val < 0.05 & deg6$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$Hotspot_vs_TP53wt_total <- NA
        ds_summary$Hotspot_vs_TP53wt_sig <- NA
        ds_summary$Hotspot_vs_TP53wt_up <- NA
        ds_summary$Hotspot_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 7: DN vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 7: DN vs TP53wt ---\n")

    dn_samples <- tp53_mut$sample_id[tp53_mut$DN == 1]
    cat("    DN samples:", length(dn_samples), " TP53wt samples:", length(wt_samples), "\n")

    group7_ids <- c(dn_samples, wt_samples)
    group7 <- setNames(
        c(rep("DN", length(dn_samples)), rep("TP53wt", length(wt_samples))),
        group7_ids
    )
    group7 <- factor(group7, levels = c("TP53wt", "DN"))

    deg7 <- run_limma_group_comparison(mat, group7, purity_vec, sa_df)

    if (!is.null(deg7)) {
        deg7 <- add_gene_symbol(deg7, ensembl2symbol)
        fwrite(deg7, file.path(ds_out, "DEG_DN_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "DN_vs_TP53wt")
        writeData(wb, 1, deg7)
        saveWorkbook(wb, file.path(ds_out, "DEG_DN_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$DN_vs_TP53wt_total <- nrow(deg7)
        ds_summary$DN_vs_TP53wt_sig <- sum(deg7$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$DN_vs_TP53wt_up <- sum(deg7$adj.P.Val < 0.05 & deg7$logFC > 0, na.rm = TRUE)
        ds_summary$DN_vs_TP53wt_down <- sum(deg7$adj.P.Val < 0.05 & deg7$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$DN_vs_TP53wt_total <- NA
        ds_summary$DN_vs_TP53wt_sig <- NA
        ds_summary$DN_vs_TP53wt_up <- NA
        ds_summary$DN_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 8: Non-DN vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 8: Non-DN vs TP53wt ---\n")

    nondn_samples <- tp53_mut$sample_id[tp53_mut$non_DN == 1]
    cat("    Non-DN samples:", length(nondn_samples), " TP53wt samples:", length(wt_samples), "\n")

    group8_ids <- c(nondn_samples, wt_samples)
    group8 <- setNames(
        c(rep("NonDN", length(nondn_samples)), rep("TP53wt", length(wt_samples))),
        group8_ids
    )
    group8 <- factor(group8, levels = c("TP53wt", "NonDN"))

    deg8 <- run_limma_group_comparison(mat, group8, purity_vec, sa_df)

    if (!is.null(deg8)) {
        deg8 <- add_gene_symbol(deg8, ensembl2symbol)
        fwrite(deg8, file.path(ds_out, "DEG_NonDN_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "NonDN_vs_TP53wt")
        writeData(wb, 1, deg8)
        saveWorkbook(wb, file.path(ds_out, "DEG_NonDN_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$NonDN_vs_TP53wt_total <- nrow(deg8)
        ds_summary$NonDN_vs_TP53wt_sig <- sum(deg8$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$NonDN_vs_TP53wt_up <- sum(deg8$adj.P.Val < 0.05 & deg8$logFC > 0, na.rm = TRUE)
        ds_summary$NonDN_vs_TP53wt_down <- sum(deg8$adj.P.Val < 0.05 & deg8$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$NonDN_vs_TP53wt_total <- NA
        ds_summary$NonDN_vs_TP53wt_sig <- NA
        ds_summary$NonDN_vs_TP53wt_up <- NA
        ds_summary$NonDN_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 9: DN vs non-DN
    # ========================================================================
    cat("\n  --- Comparison 9: DN vs non-DN ---\n")

    cat("    DN samples:", length(dn_samples), " Non-DN samples:", length(nondn_samples), "\n")

    group9_ids <- c(dn_samples, nondn_samples)
    group9_labels <- setNames(
        c(rep("DN", length(dn_samples)), rep("NonDN", length(nondn_samples))),
        group9_ids
    )
    # Handle possible overlap: a sample could be both DN and non-DN (shouldn't happen, but safe)
    group9_labels <- group9_labels[!duplicated(names(group9_labels))]
    group9 <- factor(group9_labels, levels = c("NonDN", "DN"))

    deg9 <- run_limma_group_comparison(mat, group9, purity_vec, sa_df)

    if (!is.null(deg9)) {
        deg9 <- add_gene_symbol(deg9, ensembl2symbol)
        fwrite(deg9, file.path(ds_out, "DEG_DN_vs_NonDN.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "DN_vs_NonDN")
        writeData(wb, 1, deg9)
        saveWorkbook(wb, file.path(ds_out, "DEG_DN_vs_NonDN.xlsx"), overwrite = TRUE)
        ds_summary$DN_vs_NonDN_total <- nrow(deg9)
        ds_summary$DN_vs_NonDN_sig <- sum(deg9$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$DN_vs_NonDN_up <- sum(deg9$adj.P.Val < 0.05 & deg9$logFC > 0, na.rm = TRUE)
        ds_summary$DN_vs_NonDN_down <- sum(deg9$adj.P.Val < 0.05 & deg9$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$DN_vs_NonDN_total <- NA
        ds_summary$DN_vs_NonDN_sig <- NA
        ds_summary$DN_vs_NonDN_up <- NA
        ds_summary$DN_vs_NonDN_down <- NA
    }

    all_deg_summary[[ds_folder]] <- ds_summary
    cat("\n")
}

# ==============================================================================
# Section 5: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Summary Statistics\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(all_deg_summary)

# Save summary CSV
fwrite(summary_df, file.path(output_dir, "DEG_summary_statistics.csv"))

# Save summary XLSX with formatting
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
saveWorkbook(wb_sum, file.path(output_dir, "DEG_summary_statistics.xlsx"), overwrite = TRUE)

# Print summary
cat("\n=== Differential Gene Summary (padj < 0.05) ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")

# Save covariate imputation report
if (length(imputation_report_all) > 0) {
    imp_df <- do.call(rbind, imputation_report_all)
    write.csv(imp_df, file.path(output_dir, "summary_covariate_imputation.csv"), row.names = FALSE)
    cat("Summary of covariate imputation saved to summary_covariate_imputation.csv\n")
}


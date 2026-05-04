################################################################################
# Survival Analysis for TP53 Mutation Categories
#
# Description:
#   Statistical analysis of Overall Survival (OS) and Progression-Free Survival 
#   (PFS) differences between various TP53 mutation groups across datasets.
#
# Methodology:
#   - Extracts OS and PFS data from [dataset]_survival.txt
#   - Performs Kaplan-Meier survival analysis and Log-rank test.
#   - Calculates Hazard Ratios (HR) using univariate Cox proportional hazards regression.
#   - Adjusts Log-rank P-values using Benjamini-Hochberg (FDR).
#   - Now includes Subgroup Covariates (Age, Sex, Purity, Cancer-specific clinical features).
#
# Date: 2026-05-04
################################################################################

library(dplyr)
library(survival)
library(openxlsx)
library(tidyr)
library(readr)

# Set base path
base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"
output_dir <- file.path(base_path, "survival_analysis_subgroup_adjusted")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("=== Survival Analysis ===\n")
cat("Output directory:", output_dir, "\n\n")

cancer_types <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

# Define comparisons
comparisons <- list(
  "TP53mt_vs_TP53wt"    = list(g1_name="TP53mt",   g2_name="TP53wt",  g1_cond="mt == 1",               g2_cond="wt == 1"),
  "MUT_GOF_vs_MUT_LOF"  = list(g1_name="MUT_GOF",  g2_name="MUT_LOF", g1_cond="GOF == 1",              g2_cond="LOF == 1"),
  "hotspots_vs_MUT_LOF" = list(g1_name="hotspots", g2_name="MUT_LOF", g1_cond="hotspot == 1",          g2_cond="LOF == 1"),
  "MUT_GOF_vs_TP53wt"   = list(g1_name="MUT_GOF",  g2_name="TP53wt",  g1_cond="GOF == 1",              g2_cond="wt == 1"),
  "MUT_LOF_vs_TP53wt"   = list(g1_name="MUT_LOF",  g2_name="TP53wt",  g1_cond="LOF == 1",              g2_cond="wt == 1"),
  "hotspots_vs_TP53wt"  = list(g1_name="hotspots", g2_name="TP53wt",  g1_cond="hotspot == 1",          g2_cond="wt == 1"),
  "DN_vs_TP53wt"        = list(g1_name="DN",       g2_name="TP53wt",  g1_cond="DN == 1",               g2_cond="wt == 1"),
  "Non_DN_vs_TP53wt"    = list(g1_name="Non_DN",   g2_name="TP53wt",  g1_cond="non_DN == 1 & mt == 1", g2_cond="wt == 1"),
  "DN_vs_non_DN"        = list(g1_name="DN",       g2_name="Non_DN",  g1_cond="DN == 1",               g2_cond="non_DN == 1 & mt == 1")
)

#' Get cancer-specific subgroup covariates
get_cancer_specific_covariates <- function(ds_dir, ds_cancer, sample_ids) {
    res <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
    
    match_id <- function(clin_ids, sid) {
        idx <- which(clin_ids == sid | sapply(clin_ids, function(pid) grepl(pid, sid, fixed = TRUE)))
        if (length(idx) > 0) return(idx[1])
        return(NA)
    }

    if (ds_cancer == "BRCA") {
        fp <- file.path(ds_dir, "HS_CPTAC_BRCA_2018_CLI.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin)) {
                clin$Sample.ID <- sub("^X", "", clin$Sample.ID)
                for (col_in in c("ER.Updated.Clinical.Status", "Her2.Updated.Clinical.Status")) {
                    col_out <- ifelse(col_in == "ER.Updated.Clinical.Status", "ER_status", "HER2_status")
                    res[[col_out]] <- NA_character_
                    if (col_in %in% colnames(clin)) {
                        for (i in seq_along(sample_ids)) {
                            idx <- match_id(clin$Sample.ID, sample_ids[i])
                            if (!is.na(idx)) res[i, col_out] <- as.character(clin[[col_in]][idx])
                        }
                    }
                }
            }
        }
    } else if (ds_cancer == "COAD") {
        fp <- file.path(ds_dir, "COAD_phenotype.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin) && "MSI_H" %in% colnames(clin)) {
                res$MSI <- NA_character_
                id_col <- colnames(clin)[1]
                for (i in seq_along(sample_ids)) {
                    idx <- match_id(clin[[id_col]], sample_ids[i])
                    if (!is.na(idx)) res[i, "MSI"] <- as.character(clin[["MSI_H"]][idx])
                }
            }
        }
    } else if (ds_cancer == "GBM") {
        fp <- file.path(ds_dir, "GBM_IDH_mutant_list.csv")
        res$IDH_status <- "WT"
        if (file.exists(fp)) {
            clin <- tryCatch(read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
            if (!is.null(clin) && "Tumor_Sample_Barcode" %in% colnames(clin)) {
                for (i in seq_along(sample_ids)) {
                    idx <- match_id(clin$Tumor_Sample_Barcode, sample_ids[i])
                    if (!is.na(idx)) res[i, "IDH_status"] <- "Mutant"
                }
            }
        }
    } else if (ds_cancer == "HNSCC") {
        fp <- file.path(ds_dir, "HS_CPTAC_HNSCC_CLI.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin) && "tumor_site_curated" %in% colnames(clin)) {
                res$tumor_site <- NA_character_
                for (i in seq_along(sample_ids)) {
                    idx <- match_id(clin$case_id, sample_ids[i])
                    if (!is.na(idx)) res[i, "tumor_site"] <- as.character(clin$tumor_site_curated[idx])
                }
            }
        }
    } else if (ds_cancer == "LSCC") {
        fp <- file.path(ds_dir, "LSCC_meta.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin) && "Tobacco_smoking_history" %in% colnames(clin)) {
                res$smoking_history <- NA_character_
                id_col <- colnames(clin)[1]
                for (i in seq_along(sample_ids)) {
                    idx <- match_id(clin[[id_col]], sample_ids[i])
                    if (!is.na(idx)) res[i, "smoking_history"] <- as.character(clin$Tobacco_smoking_history[idx])
                }
            }
        }
    } else if (ds_cancer == "LUAD") {
        fp <- file.path(ds_dir, "LUAD_meta.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin)) {
                cols_to_map <- c("Tobacco_smoking_history" = "smoking_history", 
                                 "EGFR_mutation" = "EGFR_mutation", 
                                 "KRAS_mutation" = "KRAS_mutation", 
                                 "STK11_mutation" = "STK11_mutation")
                id_col <- colnames(clin)[1]
                for (col_in in names(cols_to_map)) {
                    col_out <- cols_to_map[[col_in]]
                    if (col_in %in% colnames(clin)) {
                        res[[col_out]] <- NA_character_
                        for (i in seq_along(sample_ids)) {
                            idx <- match_id(clin[[id_col]], sample_ids[i])
                            if (!is.na(idx)) res[i, col_out] <- as.character(clin[[col_in]][idx])
                        }
                    }
                }
            }
        }
    } else if (ds_cancer == "PDAC") {
        fp <- file.path(ds_dir, "PDAC_meta.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin)) {
                cols_to_map <- c("KRAS_mutation" = "KRAS_mutation", 
                                 "SMAD4_mutation" = "SMAD4_mutation", 
                                 "CDKN2A_mutation" = "CDKN2A_mutation")
                id_col <- colnames(clin)[1]
                for (col_in in names(cols_to_map)) {
                    col_out <- cols_to_map[[col_in]]
                    if (col_in %in% colnames(clin)) {
                        res[[col_out]] <- NA_character_
                        for (i in seq_along(sample_ids)) {
                            idx <- match_id(clin[[id_col]], sample_ids[i])
                            if (!is.na(idx)) res[i, col_out] <- as.character(clin[[col_in]][idx])
                        }
                    }
                }
            }
        }
    } else if (ds_cancer == "UCEC") {
        fp <- file.path(ds_dir, "UCEC_phenotype.txt")
        if (file.exists(fp)) {
            clin <- tryCatch(readr::read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(clin)) {
                cols_to_map <- c("MSI_H" = "MSI", "POLE" = "POLE")
                id_col <- colnames(clin)[1]
                for (col_in in names(cols_to_map)) {
                    col_out <- cols_to_map[[col_in]]
                    if (col_in %in% colnames(clin)) {
                        res[[col_out]] <- NA_character_
                        for (i in seq_along(sample_ids)) {
                            idx <- match_id(clin[[id_col]], sample_ids[i])
                            if (!is.na(idx)) res[i, col_out] <- as.character(clin[[col_in]][idx])
                        }
                    }
                }
            }
        }
    }
    
    return(res)
}

results_list <- list()
imputation_report_all <- list()

for (cancer in cancer_types) {
  cat("\nProcessing dataset:", cancer, "\n")
  
  # 1. Read classification file
  class_file <- file.path(base_path, "TP53_mutation_classification", paste0(cancer, "_TP53_classification.csv"))
  if (!file.exists(class_file)) {
    cat("  Classification file not found. Skipping.\n")
    next
  }
  tp53_class <- read.csv(class_file, stringsAsFactors = FALSE)
  
  # 2. Read Meta file for Age and Sex
  meta_file <- file.path(base_path, cancer, paste0(cancer, "_meta.txt"))
  meta_df <- data.frame(sample_id = character())
  if (file.exists(meta_file)) {
    tmp_meta <- read.delim(meta_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(tmp_meta) > 0 && tmp_meta[1, 1] == "data_type") {
      tmp_meta <- tmp_meta[-1, ]
    }
    id_col_meta <- colnames(tmp_meta)[1]
    tmp_meta <- tmp_meta %>% rename(sample_id = all_of(id_col_meta))
    
    if ("Age" %in% colnames(tmp_meta)) tmp_meta$Age <- as.numeric(tmp_meta$Age)
    
    cols_to_keep <- c("sample_id")
    if ("Age" %in% colnames(tmp_meta)) cols_to_keep <- c(cols_to_keep, "Age")
    if ("Sex" %in% colnames(tmp_meta)) cols_to_keep <- c(cols_to_keep, "Sex")
    
    meta_df <- tmp_meta[, cols_to_keep, drop = FALSE]
  }
  
  # 4. Read survival file
  surv_file <- file.path(base_path, cancer, paste0(cancer, "_survival.txt"))
  if (!file.exists(surv_file)) {
    cat("  Survival file not found. Skipping.\n")
    next
  }
  surv_df <- read.delim(surv_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  required_cols <- c("case_id", "OS_days", "OS_event", "PFS_days", "PFS_event")
  if (!all(required_cols %in% colnames(surv_df))) {
    cat("  Missing standard survival columns. Skipping.\n")
    next
  }
  
  surv_df$OS_days <- as.numeric(surv_df$OS_days)
  surv_df$OS_event <- as.numeric(surv_df$OS_event)
  surv_df$PFS_days <- as.numeric(surv_df$PFS_days)
  surv_df$PFS_event <- as.numeric(surv_df$PFS_event)
  
  # 5. Merge Data
  merged_df <- inner_join(surv_df, tp53_class, by = c("case_id" = "sample_id"))
  if (nrow(meta_df) > 0) {
    merged_df <- left_join(merged_df, meta_df, by = c("case_id" = "sample_id"))
  }

  
  # Add extra cancer-specific covariates
  extra_cov <- get_cancer_specific_covariates(file.path(base_path, cancer), cancer, merged_df$case_id)
  extra_cov_names <- character(0)
  if (ncol(extra_cov) > 1) {
      colnames(extra_cov)[1] <- "case_id"
      merged_df <- left_join(merged_df, extra_cov, by = "case_id")
      extra_cov_names <- setdiff(colnames(extra_cov), "case_id")
  }
  
  # --- Covariate Imputation at Dataset Level ---
  n_samples <- nrow(merged_df)
  
  # Age
  a_miss_cnt <- 0; a_imp_val <- NA
  if ("Age" %in% colnames(merged_df)) {
      a <- merged_df$Age
      a_miss_cnt <- sum(is.na(a))
      if (sum(is.finite(a)) >= n_samples * 0.8) {
          a_imp_val <- median(a, na.rm = TRUE)
          a[is.na(a)] <- a_imp_val
          merged_df$Age <- a
      } else {
          merged_df$Age <- NULL
      }
  }
  
  # Categorical (Sex + extra_cov_names)
  cat_covs <- c("Sex", extra_cov_names)
  cat_covs <- intersect(cat_covs, colnames(merged_df))
  cat_miss_info <- list()
  
  for (cov in cat_covs) {
      val <- merged_df[[cov]]
      miss_idx <- is.na(val) | val == ""
      miss_cnt <- sum(miss_idx)
      
      if (sum(!miss_idx) >= n_samples * 0.5) {
          val[miss_idx] <- "Unknown"
          # Convert to factor with safe valid names to avoid modeling issues
          merged_df[[cov]] <- as.factor(make.names(val))
          cat_miss_info[[cov]] <- miss_cnt
      } else {
          merged_df[[cov]] <- NULL
          cat_miss_info[[cov]] <- "Dropped (>50% missing)"
      }
  }
  
  # Report imputation
  imp_report <- data.frame(
      Dataset = cancer,
      Total_Samples = n_samples,
      Age_Missing = a_miss_cnt,
      Age_Imputed_Value = ifelse(is.na(a_imp_val), "Dropped", round(a_imp_val, 1)),
      stringsAsFactors = FALSE
  )
  for (cov in cat_covs) {
      imp_report[[paste0(cov, "_Missing")]] <- if (!is.null(cat_miss_info[[cov]])) as.character(cat_miss_info[[cov]]) else "NA"
  }
  imputation_report_all[[cancer]] <- imp_report
  
  cat("  [Covariate Imputation Report]\n")
  cat("    Age missing:", a_miss_cnt, "-> Imputed with:", ifelse(is.na(a_imp_val), "Dropped", round(a_imp_val, 1)), "\n")
  for (cov in cat_covs) {
      cat("   ", cov, "missing:", cat_miss_info[[cov]], "-> Imputed with: Unknown\n")
  }
  
  # --- Iterate over comparisons ---
  for (comp_name in names(comparisons)) {
    comp <- comparisons[[comp_name]]
    
    idx_g1 <- eval(parse(text = comp$g1_cond), merged_df)
    idx_g2 <- eval(parse(text = comp$g2_cond), merged_df)
    
    comp_data <- merged_df[idx_g1 | idx_g2, ]
    
    if (nrow(comp_data) == 0) next
    
    comp_data$Group <- NA
    comp_data$Group[eval(parse(text = comp$g1_cond), comp_data)] <- "G1"
    comp_data$Group[eval(parse(text = comp$g2_cond), comp_data)] <- "G2"
    comp_data <- comp_data[!is.na(comp_data$Group), ]
    
    if (length(unique(comp_data$Group)) < 2) next
    
    # --- Function to determine valid covariates ---
    get_valid_covariates <- function(data) {
      valid_cov <- c()
      if ("Age" %in% colnames(data) && var(data$Age, na.rm = TRUE) > 0) valid_cov <- c(valid_cov, "Age")
      for (cov in cat_covs) {
          if (cov %in% colnames(data) && length(unique(data[[cov]])) > 1) {
              valid_cov <- c(valid_cov, cov)
          }
      }
      return(valid_cov)
    }

    # --- Overall Survival (OS) Analysis ---
    os_data <- comp_data[!is.na(comp_data$OS_days) & !is.na(comp_data$OS_event), ]
    os_data <- os_data[os_data$OS_days >= 0, ]
    
    if (nrow(os_data) >= 10 && length(unique(os_data$Group)) == 2) {
      os_data$Group <- factor(os_data$Group, levels = c("G2", "G1")) 
      surv_obj <- Surv(time = os_data$OS_days, event = os_data$OS_event)
      
      p_logrank <- tryCatch({
        sdiff <- survdiff(surv_obj ~ Group, data = os_data)
        1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
      }, error = function(e) NA)
      
      valid_cov_os <- get_valid_covariates(os_data)
      formula_str <- "surv_obj ~ Group"
      if (length(valid_cov_os) > 0) formula_str <- paste(formula_str, "+", paste(valid_cov_os, collapse = " + "))
      form <- as.formula(formula_str)
      
      cox_fit <- tryCatch(coxph(form, data = os_data), error=function(e) NULL)
      hr <- NA; hr_lower <- NA; hr_upper <- NA; p_cox <- NA
      if (!is.null(cox_fit)) {
        s <- summary(cox_fit)
        if ("GroupG1" %in% rownames(s$conf.int)) {
          hr <- s$conf.int["GroupG1", "exp(coef)"]
          hr_lower <- s$conf.int["GroupG1", "lower .95"]
          hr_upper <- s$conf.int["GroupG1", "upper .95"]
          p_cox <- s$coefficients["GroupG1", "Pr(>|z|)"]
        }
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Dataset = cancer,
        Comparison = comp_name,
        Group1 = comp$g1_name,
        Group2 = comp$g2_name,
        Survival_Type = "OS",
        Covariates_Adjusted = paste(valid_cov_os, collapse = ", "),
        N_Total = nrow(os_data),
        N_G1 = sum(os_data$Group == "G1"),
        N_G2 = sum(os_data$Group == "G2"),
        Events_G1 = sum(os_data$OS_event[os_data$Group == "G1"]),
        Events_G2 = sum(os_data$OS_event[os_data$Group == "G2"]),
        HR_Adjusted = hr,
        HR_CI95_Lower = hr_lower,
        HR_CI95_Upper = hr_upper,
        P_Value_Cox = p_cox,
        P_Value_LogRank_Unadjusted = p_logrank,
        stringsAsFactors = FALSE
      )
    }
    
    # --- Progression-Free Survival (PFS) Analysis ---
    pfs_data <- comp_data[!is.na(comp_data$PFS_days) & !is.na(comp_data$PFS_event), ]
    pfs_data <- pfs_data[pfs_data$PFS_days >= 0, ]
    
    if (nrow(pfs_data) >= 10 && length(unique(pfs_data$Group)) == 2) {
      pfs_data$Group <- factor(pfs_data$Group, levels = c("G2", "G1")) 
      surv_obj <- Surv(time = pfs_data$PFS_days, event = pfs_data$PFS_event)
      
      p_logrank <- tryCatch({
        sdiff <- survdiff(surv_obj ~ Group, data = pfs_data)
        1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
      }, error = function(e) NA)
      
      valid_cov_pfs <- get_valid_covariates(pfs_data)
      formula_str <- "surv_obj ~ Group"
      if (length(valid_cov_pfs) > 0) formula_str <- paste(formula_str, "+", paste(valid_cov_pfs, collapse = " + "))
      form <- as.formula(formula_str)
      
      cox_fit <- tryCatch(coxph(form, data = pfs_data), error=function(e) NULL)
      hr <- NA; hr_lower <- NA; hr_upper <- NA; p_cox <- NA
      if (!is.null(cox_fit)) {
        s <- summary(cox_fit)
        if ("GroupG1" %in% rownames(s$conf.int)) {
          hr <- s$conf.int["GroupG1", "exp(coef)"]
          hr_lower <- s$conf.int["GroupG1", "lower .95"]
          hr_upper <- s$conf.int["GroupG1", "upper .95"]
          p_cox <- s$coefficients["GroupG1", "Pr(>|z|)"]
        }
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Dataset = cancer,
        Comparison = comp_name,
        Group1 = comp$g1_name,
        Group2 = comp$g2_name,
        Survival_Type = "PFS",
        Covariates_Adjusted = paste(valid_cov_pfs, collapse = ", "),
        N_Total = nrow(pfs_data),
        N_G1 = sum(pfs_data$Group == "G1"),
        N_G2 = sum(pfs_data$Group == "G2"),
        Events_G1 = sum(pfs_data$PFS_event[pfs_data$Group == "G1"]),
        Events_G2 = sum(pfs_data$PFS_event[pfs_data$Group == "G2"]),
        HR_Adjusted = hr,
        HR_CI95_Lower = hr_lower,
        HR_CI95_Upper = hr_upper,
        P_Value_Cox = p_cox,
        P_Value_LogRank_Unadjusted = p_logrank,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(imputation_report_all) > 0) {
  imp_df <- bind_rows(imputation_report_all)
  write.csv(imp_df, file.path(output_dir, "summary_covariate_imputation.csv"), row.names = FALSE)
}

if (length(results_list) > 0) {
  cat("\nAggregating results and calculating FDR...\n")
  final_results <- bind_rows(results_list)
  
  final_results <- final_results %>%
    group_by(Dataset, Survival_Type) %>%
    mutate(FDR_Cox = p.adjust(P_Value_Cox, method = "BH")) %>%
    ungroup() %>%
    arrange(Dataset, Survival_Type, P_Value_Cox)
    
  csv_file <- file.path(output_dir, "survival_analysis_summary.csv")
  write.csv(final_results, csv_file, row.names = FALSE)
  
  xlsx_file <- file.path(output_dir, "survival_analysis_summary.xlsx")
  wb <- createWorkbook()
  addWorksheet(wb, "Survival")
  writeData(wb, "Survival", final_results)
  
  header_style <- createStyle(textDecoration = "bold", fgFill = "#4F81BD", fontColour = "white")
  addStyle(wb, "Survival", header_style, rows = 1, cols = 1:ncol(final_results), gridExpand = TRUE)
  
  sig_style <- createStyle(fontColour = "#C00000", textDecoration = "bold")
  
  if (any(!is.na(final_results$FDR_Cox) & final_results$FDR_Cox < 0.05)) {
    sig_idx <- which(!is.na(final_results$FDR_Cox) & final_results$FDR_Cox < 0.05) + 1
    addStyle(wb, "Survival", sig_style, rows = sig_idx, cols = which(colnames(final_results) == "FDR_Cox"), gridExpand = FALSE)
  }
  
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  cat("Completed successfully! Files saved to:\n  -", csv_file, "\n  -", xlsx_file, "\n")
} else {
  cat("\nNo valid statistical comparisons could be performed.\n")
}

################################################################################
# Phenotype Association Analysis for TP53 Mutation Categories (Covariate Adjusted)
#
# Description:
#   Statistical analysis of phenotypic differences between various TP53 
#   mutation groups across datasets, adjusting for Age, Sex, and WES_purity.
#
# Methodology:
#   - Continuous phenotypes: Multiple Linear Regression (lm)
#     Phenotype ~ TP53_Group + Age + Sex + WES_purity
#   - Binary phenotypes: Logistic Regression (glm)
#     Phenotype ~ TP53_Group + Age + Sex + WES_purity
#   - Multi-class phenotypes: Fisher's exact test (fallback, no covariates)
#   P-values are adjusted using Benjamini-Hochberg (FDR) within each dataset, 
#   comparison pair, and phenotype category.
#
# Date: 2026-04-24
################################################################################

library(dplyr)
library(readr)
library(openxlsx)
library(tidyr)

# Set base path
base_path <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF"
output_dir <- file.path(base_path, "phenotype_association_subgroup_adjusted")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("=== Phenotype Association Analysis (Covariate Adjusted) ===\n")
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

all_results_list <- list()
imputation_report_all <- list()

for (cancer in cancer_types) {
  cat("Processing dataset:", cancer, "\n")
  
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
  
  # 3. Read phenotype file
  pheno_file <- file.path(base_path, cancer, paste0(cancer, "_phenotype.txt"))
  if (!file.exists(pheno_file)) {
    cat("  Phenotype file not found. Skipping.\n")
    next
  }
  pheno_df <- read.delim(pheno_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- colnames(pheno_df)[1]
  
  # Ensure WES_purity is numeric if it exists
  if ("WES_purity" %in% colnames(pheno_df)) {
    pheno_df$WES_purity <- as.numeric(pheno_df$WES_purity)
  }
  
  # 4. Merge Data
  merged_df <- inner_join(pheno_df, tp53_class, by = setNames("sample_id", id_col))
  if (nrow(meta_df) > 0) {
    merged_df <- left_join(merged_df, meta_df, by = setNames("sample_id", id_col))
  }
  
  # Add extra covariates
  extra_cov <- get_cancer_specific_covariates(file.path(base_path, cancer), cancer, merged_df[[id_col]])
  extra_cov_names <- character(0)
  if (ncol(extra_cov) > 1) {
      colnames(extra_cov)[1] <- id_col
      merged_df <- left_join(merged_df, extra_cov, by = id_col)
      extra_cov_names <- setdiff(colnames(extra_cov), id_col)
  }
  
  # Identify phenotype columns (excluding clinical covariates and id)
  pheno_cols <- setdiff(colnames(pheno_df), c(id_col, "Age", "Sex", "WES_purity", extra_cov_names))
  
  # --- Covariate Imputation at Dataset Level ---
  n_samples <- nrow(merged_df)
  
  # Purity
  p_miss_cnt <- 0; p_imp_val <- NA
  if ("WES_purity" %in% colnames(merged_df)) {
      p <- merged_df$WES_purity
      p_miss_cnt <- sum(is.na(p))
      if (sum(is.finite(p)) >= n_samples * 0.6) {
          p_imp_val <- median(p, na.rm = TRUE)
          p[is.na(p)] <- p_imp_val
          merged_df$WES_purity <- p
      } else {
          merged_df$WES_purity <- NULL
      }
  }
  
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
  
  # Sex
  s_miss_cnt <- 0; s_imp_val <- NA
  if ("Sex" %in% colnames(merged_df)) {
      s <- merged_df$Sex
      s_miss_cnt <- sum(is.na(s) | s == "Unknown")
      if (sum(!(is.na(s) | s == "Unknown")) >= n_samples * 0.8) {
          s_imp_val <- "Unknown"
          if (s_miss_cnt > 0) {
              s[is.na(s)] <- s_imp_val
              merged_df$Sex <- s
          }
      } else {
          merged_df$Sex <- NULL
      }
  }
  
  imputation_report_row <- data.frame(
      dataset = cancer,
      total_samples = n_samples,
      purity_missing_rate = sprintf("%.1f%%", p_miss_cnt / n_samples * 100),
      purity_imputed_value = as.character(ifelse(is.na(p_imp_val), "Dropped (>40% NA)", ifelse(p_miss_cnt == 0, "None", round(p_imp_val, 4)))),
      age_missing_rate = sprintf("%.1f%%", a_miss_cnt / n_samples * 100),
      age_imputed_value = as.character(ifelse(is.na(a_imp_val), "Dropped (>20% NA)", ifelse(a_miss_cnt == 0, "None", round(a_imp_val, 1)))),
      sex_missing_rate = sprintf("%.1f%%", s_miss_cnt / n_samples * 100),
      sex_imputed_value = as.character(ifelse(is.na(s_imp_val), "Dropped (>20% NA)", ifelse(s_miss_cnt == 0, "None", s_imp_val))),
      stringsAsFactors = FALSE
  )
  
  cat("\n  [Covariate Imputation Report]\n")
  cat(sprintf("    Purity: Missing %s | Imputed value: %s\n", imputation_report_row$purity_missing_rate, imputation_report_row$purity_imputed_value))
  cat(sprintf("    Age   : Missing %s | Imputed value: %s\n", imputation_report_row$age_missing_rate, imputation_report_row$age_imputed_value))
  cat(sprintf("    Sex   : Missing %s | Imputed value: %s\n", imputation_report_row$sex_missing_rate, imputation_report_row$sex_imputed_value))
  
  # Extra covariates imputation
  for (ec in extra_cov_names) {
      v <- merged_df[[ec]]
      miss_cnt <- sum(is.na(v) | v == "Unknown")
      imp_val <- NA
      if (sum(!(is.na(v) | v == "Unknown")) >= n_samples * 0.8) {
          imp_val <- "Unknown"
          if (miss_cnt > 0) {
              v[is.na(v)] <- imp_val
              merged_df[[ec]] <- v
          }
      } else {
          merged_df[[ec]] <- NULL
          extra_cov_names <- setdiff(extra_cov_names, ec)
      }
      
      miss_rate_str <- sprintf("%.1f%%", miss_cnt / n_samples * 100)
      imp_val_str <- as.character(ifelse(is.na(imp_val), "Dropped (>20% NA)", ifelse(miss_cnt == 0, "None", imp_val)))
      
      imputation_report_row[[paste0(ec, "_missing_rate")]] <- miss_rate_str
      imputation_report_row[[paste0(ec, "_imputed_value")]] <- imp_val_str
      
      cat(sprintf("    %s: Missing %s | Imputed value: %s\n", ec, miss_rate_str, imp_val_str))
  }
  
  imputation_report_all[[cancer]] <- imputation_report_row
  cat("  --------------------------------------\n")
  
  # Perform tests for each comparison
  for (comp_name in names(comparisons)) {
    comp <- comparisons[[comp_name]]
    
    g1_df <- subset(merged_df, eval(parse(text = comp$g1_cond)))
    g2_df <- subset(merged_df, eval(parse(text = comp$g2_cond)))
    
    if (nrow(g1_df) < 3 || nrow(g2_df) < 3) {
      cat("  - Skipping", comp_name, "due to low sample size\n")
      next
    }
    
    for (feature in pheno_cols) {
      # Prepare dataset for regression
      cov_cols <- intersect(colnames(merged_df), c("Age", "Sex", "WES_purity", extra_cov_names))
      
      ana_df <- bind_rows(
        g1_df %>% select(all_of(feature), all_of(cov_cols)) %>% mutate(Group = "G1"),
        g2_df %>% select(all_of(feature), all_of(cov_cols)) %>% mutate(Group = "G2")
      )
      
      # Remove NA for the target feature
      ana_df <- ana_df[!is.na(ana_df[[feature]]) & ana_df[[feature]] != "", ]
      
      n_g1 <- sum(ana_df$Group == "G1")
      n_g2 <- sum(ana_df$Group == "G2")
      if (n_g1 < 3 || n_g2 < 3) next
      
      ana_df$Group <- factor(ana_df$Group, levels = c("G2", "G1")) # G2 is reference
      
      # Determine valid covariates for this specific subset
      valid_covariates <- c()
      if ("Age" %in% colnames(ana_df) && var(ana_df$Age, na.rm = TRUE) > 0) {
        valid_covariates <- c(valid_covariates, "Age")
      }
      if ("Sex" %in% colnames(ana_df) && length(unique(na.omit(ana_df$Sex))) > 1) {
        valid_covariates <- c(valid_covariates, "Sex")
      }
      if ("WES_purity" %in% colnames(ana_df) && var(ana_df$WES_purity, na.rm = TRUE) > 0) {
        valid_covariates <- c(valid_covariates, "WES_purity")
      }
      for (ec in extra_cov_names) {
          if (ec %in% colnames(ana_df) && length(unique(na.omit(ana_df[[ec]]))) > 1) {
              valid_covariates <- c(valid_covariates, ec)
          }
      }
      
      # Construct formula
      formula_str <- paste0("`", feature, "` ~ Group")
      if (length(valid_covariates) > 0) {
        formula_str <- paste0(formula_str, " + ", paste(valid_covariates, collapse = " + "))
      }
      form <- as.formula(formula_str)
      
      v1 <- ana_df[[feature]][ana_df$Group == "G1"]
      v2 <- ana_df[[feature]][ana_df$Group == "G2"]
      is_num <- is.numeric(v1) && is.numeric(v2)
      
      mean1 <- if(is_num) mean(v1, na.rm=TRUE) else NA
      mean2 <- if(is_num) mean(v2, na.rm=TRUE) else NA
      diff_val <- if(is_num) mean1 - mean2 else NA
      
      p_val <- NA
      
      if (is_num) {
        if (var(ana_df[[feature]], na.rm = TRUE) == 0) next
        
        # Multiple Linear Regression
        tryCatch({
          fit <- lm(form, data = ana_df)
          coef_table <- coef(summary(fit))
          if ("GroupG1" %in% rownames(coef_table)) {
            p_val <- coef_table["GroupG1", "Pr(>|t|)"]
          }
        }, error = function(e) NULL)
        
        data_type <- "Numeric (lm)"
        
      } else {
        val_levels <- unique(as.character(ana_df[[feature]]))
        if (length(val_levels) == 2) {
          # Binary Logistic Regression
          ana_df$feature_bin <- factor(ana_df[[feature]], levels = val_levels)
          form_bin <- as.formula(paste("feature_bin ~ Group", 
                                       if(length(valid_covariates)>0) paste("+", paste(valid_covariates, collapse=" + ")) else ""))
          tryCatch({
            fit <- glm(form_bin, data = ana_df, family = binomial)
            coef_table <- coef(summary(fit))
            if ("GroupG1" %in% rownames(coef_table)) {
              p_val <- coef_table["GroupG1", "Pr(>|z|)"]
            }
          }, error = function(e) NULL)
          data_type <- "Binary (glm)"
          
        } else {
          # Multi-class Categorical -> Fallback to Fisher/Chi-sq without covariates
          tab <- table(ana_df$Group, ana_df[[feature]])
          tryCatch({
            p_val <- fisher.test(tab, simulate.p.value = TRUE)$p.value
          }, error = function(e) {
            tryCatch({ p_val <- chisq.test(tab)$p.value }, error = function(e2) NULL)
          })
          data_type <- "Multi-class (Fisher/Chi-sq)"
        }
      }
      
      if (!is.na(p_val)) {
        all_results_list[[length(all_results_list) + 1]] <- data.frame(
          Dataset = cancer,
          Comparison = comp_name,
          Group1 = comp$g1_name,
          Group2 = comp$g2_name,
          Phenotype = feature,
          Data_Type = data_type,
          Covariates_Adjusted = paste(valid_covariates, collapse = ", "),
          N_Group1 = n_g1,
          N_Group2 = n_g2,
          Mean_Group1 = mean1,
          Mean_Group2 = mean2,
          Difference = diff_val,
          P_Value = p_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Aggregate results and perform Multiple Testing Correction
if (length(all_results_list) > 0) {
  cat("\nAggregating results and calculating FDR...\n")
  final_results <- bind_rows(all_results_list)
  
  # Assign category based on phenotype prefix
  final_results <- final_results %>%
    mutate(
      Phenotype_Category = case_when(
        grepl("^CIBERSORT_", Phenotype) ~ "CIBERSORT",
        grepl("^ESTIMATE_", Phenotype) ~ "ESTIMATE",
        grepl("^PROGENy_", Phenotype) ~ "PROGENy",
        grepl("^HALLMARK_", Phenotype) ~ "HALLMARK",
        grepl("^xCell_", Phenotype) ~ "xCell",
        grepl("^PTM_SEA_", Phenotype) ~ "PTM_SEA",
        grepl("^CNV_index_", Phenotype) ~ "CNV_index",
        grepl("^Mutation_signature_", Phenotype) ~ "Mutation_Signature",
        TRUE ~ "Other_Clinical"
      )
    ) %>%
    relocate(Phenotype_Category, .after = Phenotype)

  # Calculate False Discovery Rate (FDR) using Benjamini-Hochberg procedure within each Category
  final_results <- final_results %>%
    group_by(Dataset, Comparison, Phenotype_Category) %>%
    mutate(FDR = p.adjust(P_Value, method = "BH")) %>%
    ungroup() %>%
    arrange(Dataset, Comparison, Phenotype_Category, P_Value)
    
  # Write output files
  csv_file <- file.path(output_dir, "phenotype_association_summary.csv")
  write.csv(final_results, csv_file, row.names = FALSE)
  
  xlsx_file <- file.path(output_dir, "phenotype_association_summary.xlsx")
  wb <- createWorkbook()
  addWorksheet(wb, "Associations")
  writeData(wb, "Associations", final_results)
  
  # Formatting styles for the Excel workbook
  header_style <- createStyle(textDecoration = "bold", fgFill = "#4F81BD", fontColour = "white")
  addStyle(wb, "Associations", header_style, rows = 1, cols = 1:ncol(final_results), gridExpand = TRUE)
  
  sig_style <- createStyle(fontColour = "#C00000", textDecoration = "bold")
  
  # Highlight significant rows (FDR < 0.05)
  if (any(!is.na(final_results$FDR) & final_results$FDR < 0.05)) {
    sig_idx <- which(!is.na(final_results$FDR) & final_results$FDR < 0.05) + 1
    addStyle(wb, "Associations", sig_style, rows = sig_idx, cols = which(colnames(final_results) == "FDR"), gridExpand = FALSE)
  }
  
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  cat("Completed successfully! Files saved to:\n  -", csv_file, "\n  -", xlsx_file, "\n")
  
} else {
  cat("\nNo valid statistical comparisons could be performed.\n")
}

# Save covariate imputation report
if (length(imputation_report_all) > 0) {
    imp_df <- bind_rows(imputation_report_all)
    write.csv(imp_df, file.path(output_dir, "summary_covariate_imputation.csv"), row.names = FALSE)
    cat("Summary of covariate imputation saved to summary_covariate_imputation.csv\n")
}

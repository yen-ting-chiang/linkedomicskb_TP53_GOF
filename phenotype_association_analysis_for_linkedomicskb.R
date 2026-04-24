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
output_dir <- file.path(base_path, "phenotype_association")

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

all_results_list <- list()

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
  
  # Identify phenotype columns (excluding clinical covariates if they exist in pheno)
  pheno_cols <- setdiff(colnames(pheno_df), c(id_col, "Age", "Sex"))
  
  # Ensure WES_purity is numeric if it exists
  if ("WES_purity" %in% colnames(pheno_df)) {
    pheno_df$WES_purity <- as.numeric(pheno_df$WES_purity)
  }
  
  # 4. Merge Data
  merged_df <- inner_join(pheno_df, tp53_class, by = setNames("sample_id", id_col))
  if (nrow(meta_df) > 0) {
    merged_df <- left_join(merged_df, meta_df, by = setNames("sample_id", id_col))
  }
  
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
      cov_cols <- intersect(colnames(merged_df), c("Age", "Sex", "WES_purity"))
      
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

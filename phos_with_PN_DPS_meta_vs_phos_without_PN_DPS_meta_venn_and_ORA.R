# phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta_venn_and_ORA.R
#
# Purpose:
#   Generate Venn diagrams and Over-Representation Analysis (ORA) for the
#   integrated phosphosite-level meta-analysis results comparing:
#     - Phosphoprotein DPS with protein normalization (phos_with_PN_DPS_meta)
#     - Phosphoprotein DPS without protein normalization (phos_without_PN_DPS_meta)
#
#   Role mapping (adapted from mRNA_vs_protein_venn_diagram_and_ORA.R):
#     mRNA role   -> phos_without_PN_DPS_meta (without_PN_Z_meta, without_PN_padj)
#     protein role -> phos_with_PN_DPS_meta   (with_PN_Z_meta, with_PN_padj)
#
#   Venn diagrams are constructed at the phosphosite level using the
#   "phosphosite" column. For ORA, the parent gene is extracted from the
#   "gene_symbol_phosphosite" column (e.g., "MCM6_S762" -> "MCM6"), duplicates
#   are removed, and the resulting unique gene list is used for enrichment
#   analysis with MSigDB Hallmark gene sets.
#
# Input:
#   - phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta/
#       META_integrated_TP53mt_vs_TP53wt.csv
#     (output from phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta.R)
#
# Output:
#   - phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta/venn_diagram/
#       Venn_UP.pdf, Venn_UP.tiff
#       Venn_DOWN.pdf, Venn_DOWN.tiff
#       Venn_UP_regions.csv, Venn_UP_regions.xlsx
#       Venn_DOWN_regions.csv, Venn_DOWN_regions.xlsx
#       ORA_Hallmark_*.csv, ORA_Hallmark_*.xlsx
#       ORA_Hallmark_*_bubble.pdf, ORA_Hallmark_*_bubble.tiff
#
# References:
#   - Liberzon A, et al. The Molecular Signatures Database (MSigDb) hallmark
#     gene set collection. Cell Syst. 2015;1(6):417-425. PMID: 26771021
#   - Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for
#     interpreting omics data. Innovation (Camb). 2021;2(3):100141.
#     PMID: 34557778

# Load necessary libraries
library(dplyr)
library(openxlsx)
library(VennDiagram)
library(grid)
library(futile.logger)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)

# Suppress VennDiagram logger output to prevent cluttering the console
flog.threshold(ERROR, name = "VennDiagramLogger")

# Define paths
input_dir <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF/phos_with_PN_DPS_meta_vs_phos_without_PN_DPS_meta"
output_dir <- file.path(input_dir, "venn_diagram")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read data
data_file <- file.path(input_dir, "META_integrated_TP53mt_vs_TP53wt.csv")
if (!file.exists(data_file)) {
  stop("Input file META_integrated_TP53mt_vs_TP53wt.csv not found.")
}
meta_data <- read.csv(data_file, stringsAsFactors = FALSE)

# Filter out rows without gene_symbol_phosphosite
meta_data <- meta_data %>% filter(!is.na(gene_symbol_phosphosite) & gene_symbol_phosphosite != "")

# Load MSigDB Hallmark gene sets for ORA
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Helper function: extract parent gene from gene_symbol_phosphosite
# e.g., "MCM6_S762" -> "MCM6", "TP53BP1_S1758" -> "TP53BP1"
extract_parent_gene <- function(gene_symbol_phosphosite) {
  sub("_[^_]*$", "", gene_symbol_phosphosite)
}

# Function to perform ORA and plot results
# Input: phosphosite_list = vector of gene_symbol_phosphosite values
#        The parent genes are extracted and deduplicated before ORA
perform_ora_and_plot <- function(phosphosite_list, region_name, out_dir) {
  phosphosite_list <- phosphosite_list[!is.na(phosphosite_list) & phosphosite_list != ""]
  if (length(phosphosite_list) == 0) return(NULL)
  
  # Extract parent genes and remove duplicates
  gene_list <- unique(extract_parent_gene(phosphosite_list))
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list) == 0) return(NULL)
  
  cat("  ORA for", region_name, ": ", length(phosphosite_list), "phosphosites ->",
      length(gene_list), "unique parent genes\n")
  
  # Perform ORA using clusterProfiler
  # Set cutoffs to 1 to ensure we can pad up to 5 even if non-significant
  ora_res <- enricher(gene = gene_list, TERM2GENE = msig_h, pvalueCutoff = 1, qvalueCutoff = 1)
  
  if (!is.null(ora_res) && nrow(as.data.frame(ora_res)) > 0) {
    res_df <- as.data.frame(ora_res)
    
    # Calculate numeric GeneRatio for plotting
    res_df$GeneRatioNum <- sapply(res_df$GeneRatio, function(x) {
      parts <- as.numeric(strsplit(x, "/")[[1]])
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
    
    if (nrow(plot_df) > 0) {
      # Save tables
      write_df <- res_df[res_df$p.adjust < 0.05, ]
      if (nrow(write_df) < 5) write_df <- head(res_df, 5)
      
      write.csv(write_df, file.path(out_dir, paste0("ORA_Hallmark_", region_name, ".csv")), row.names = FALSE)
      write.xlsx(write_df, file.path(out_dir, paste0("ORA_Hallmark_", region_name, ".xlsx")), row.names = FALSE)
      
      # Clean Pathway names: remove "HALLMARK_" or "Hallmark "
      plot_df$Description <- gsub("(?i)^hallmark[_ ]", "", plot_df$Description)
      
      # Order Description factor by p.adjust (descending so smallest padj is at the top)
      plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
      
      # Control plot title visibility (TRUE to show, FALSE to hide)
      show_plot_title <- TRUE
      plot_title_text <- if(show_plot_title) paste("Hallmark ORA:", region_name) else NULL
      
      # Bubble plot (custom ggplot)
      p <- ggplot(plot_df, aes(x = GeneRatioNum, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        # scale_color_gradient controls the p.adjust color legend, use breaks to specify tick marks
        #scale_color_viridis_c(option = "viridis", direction = -1, name = "p.adjust")+
        scale_color_gradient(low = "#E64B35", high = "#4DBBD5", name = "p.adjust")+
        labs(title = plot_title_text, x = "GeneRatio", y = "") +
        # Expand axis limits to prevent bubbles at edges from being cut off or overlapping the legend
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
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5, color = "black"),
          # legend size and position
          legend.position = "right",                                             # Ensure legend is placed outside on the right
          legend.box.margin = margin(0, 0, 0, 10),                               # Add space between main plot and legend to prevent overlapping
          legend.title = element_text(size = 7, face = "bold", color = "black"), # Legend title size
          legend.text = element_text(size = 7, color = "black"),                 # Legend text size
          legend.key.size = unit(0.3, "cm")                                      # Legend key icon size
        )
      
      if(!show_plot_title) {
        p <- p + theme(plot.title = element_blank())
      }
      
      pdf(file.path(out_dir, paste0("ORA_Hallmark_", region_name, "_bubble.pdf")), width = 3.8, height = 2.5)
      print(p)
      dev.off()
      
      tiff(file.path(out_dir, paste0("ORA_Hallmark_", region_name, "_bubble.tiff")), width = 3.8, height = 2.5, units = "in", res = 300, compression = "lzw")
      print(p)
      dev.off()
      
      cat("ORA results and bubble plots generated for", region_name, "\n")
    }
  } else {
    cat("No significant ORA results for", region_name, "\n")
  }
}

# ==============================================================================
# 1. Process UP-regulated phosphosites
# ==============================================================================

# Extract lists for Venn diagram 1 (UP)
# mRNA role -> phos_without_PN_DPS_meta (without_PN_Z_meta > 0 & without_PN_padj < 0.05)
phos_without_PN_up_list <- meta_data %>%
  filter(without_PN_Z_meta > 0 & without_PN_padj < 0.05) %>%
  pull(gene_symbol_phosphosite) %>%
  unique()

# protein role -> phos_with_PN_DPS_meta (with_PN_Z_meta > 0 & with_PN_padj < 0.05)
phos_with_PN_up_list <- meta_data %>%
  filter(with_PN_Z_meta > 0 & with_PN_padj < 0.05) %>%
  pull(gene_symbol_phosphosite) %>%
  unique()

# Function to get Venn intersections and save as table
get_venn_regions <- function(listA, listB, nameA, nameB) {
  onlyA <- setdiff(listA, listB)
  onlyB <- setdiff(listB, listA)
  intersection <- intersect(listA, listB)
  
  max_len <- max(length(onlyA), length(onlyB), length(intersection))
  
  length(onlyA) <- max_len
  length(onlyB) <- max_len
  length(intersection) <- max_len
  
  res <- data.frame(
    Only_A = onlyA,
    Intersection = intersection,
    Only_B = onlyB,
    stringsAsFactors = FALSE
  )
  colnames(res) <- c(paste0("Only_", nameA), "Intersection", paste0("Only_", nameB))
  return(res)
}

up_regions <- get_venn_regions(phos_without_PN_up_list, phos_with_PN_up_list,
                               "phos_without_PN", "phos_with_PN")
write.csv(up_regions, file.path(output_dir, "Venn_UP_regions.csv"), row.names = FALSE, na = "")
write.xlsx(up_regions, file.path(output_dir, "Venn_UP_regions.xlsx"), row.names = FALSE)

# Plot UP Venn Diagram
# Using NPG (Nature Publishing Group) colors for high-impact journal style
color_up <- c("#E64B35FF", "#4DBBD5FF") # Nature Red and Nature Blue

venn.plot.up <- venn.diagram(
  x = list(phos_without_PN_UP = phos_without_PN_up_list,
           phos_with_PN_UP = phos_with_PN_up_list),
  filename = NULL,
  lwd = 2,
  lty = 'solid',
  col = color_up,
  fill = color_up,
  alpha = 0.6,
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-45, 45),
  cat.dist = c(0.15, 0.15),
  cex = 0.4,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  scaled = FALSE, # Ensure both circles are the same size
  main = "Up-regulated Phosphosites\n(phos_without_PN vs phos_with_PN)",
  main.cex = 0.4,
  main.fontface = "bold",
  margin = 0.1
)

# Save UP Venn to PDF
pdf(file.path(output_dir, "Venn_UP.pdf"), width = 1.5, height = 1.5)
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# Save UP Venn to TIFF
tiff(file.path(output_dir, "Venn_UP.tiff"), width = 1.5, height = 1.5, units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# ==============================================================================
# 1.5. ORA for UP-regulated regions (using parent genes, deduplicated)
# ==============================================================================
perform_ora_and_plot(up_regions$Only_phos_without_PN, "UP_Only_phos_without_PN", output_dir)
perform_ora_and_plot(up_regions$Intersection, "UP_Intersection", output_dir)
perform_ora_and_plot(up_regions$Only_phos_with_PN, "UP_Only_phos_with_PN", output_dir)

# ==============================================================================
# 2. Process DOWN-regulated phosphosites
# ==============================================================================

# Extract lists for Venn diagram 2 (DOWN)
# mRNA role -> phos_without_PN_DPS_meta (without_PN_Z_meta < 0 & without_PN_padj < 0.05)
phos_without_PN_down_list <- meta_data %>%
  filter(without_PN_Z_meta < 0 & without_PN_padj < 0.05) %>%
  pull(gene_symbol_phosphosite) %>%
  unique()

# protein role -> phos_with_PN_DPS_meta (with_PN_Z_meta < 0 & with_PN_padj < 0.05)
phos_with_PN_down_list <- meta_data %>%
  filter(with_PN_Z_meta < 0 & with_PN_padj < 0.05) %>%
  pull(gene_symbol_phosphosite) %>%
  unique()

down_regions <- get_venn_regions(phos_without_PN_down_list, phos_with_PN_down_list,
                                 "phos_without_PN", "phos_with_PN")
write.csv(down_regions, file.path(output_dir, "Venn_DOWN_regions.csv"), row.names = FALSE, na = "")
write.xlsx(down_regions, file.path(output_dir, "Venn_DOWN_regions.xlsx"), row.names = FALSE)

# Plot DOWN Venn Diagram
# Using another set of NPG colors
color_down <- c("#00A087FF", "#3C5488FF") # Nature Green and Nature Dark Blue

venn.plot.down <- venn.diagram(
  x = list(phos_without_PN_DOWN = phos_without_PN_down_list,
           phos_with_PN_DOWN = phos_with_PN_down_list),
  filename = NULL,
  lwd = 2,
  lty = 'solid',
  col = color_down,
  fill = color_down,
  alpha = 0.6,
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-45, 45),
  cat.dist = c(0.15, 0.15),
  cex = 0.4,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  scaled = FALSE, # Ensure both circles are the same size
  main = "Down-regulated Phosphosites\n(phos_without_PN vs phos_with_PN)",
  main.cex = 0.4,
  main.fontface = "bold",
  margin = 0.1
)

# Save DOWN Venn to PDF
pdf(file.path(output_dir, "Venn_DOWN.pdf"), width = 1.5, height = 1.5)
grid.newpage()
grid.draw(venn.plot.down)
dev.off()

# Save DOWN Venn to TIFF
tiff(file.path(output_dir, "Venn_DOWN.tiff"), width = 1.5, height = 1.5, units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.down)
dev.off()

# ==============================================================================
# 2.5. ORA for DOWN-regulated regions (using parent genes, deduplicated)
# ==============================================================================
perform_ora_and_plot(down_regions$Only_phos_without_PN, "DOWN_Only_phos_without_PN", output_dir)
perform_ora_and_plot(down_regions$Intersection, "DOWN_Intersection", output_dir)
perform_ora_and_plot(down_regions$Only_phos_with_PN, "DOWN_Only_phos_with_PN", output_dir)

print("Venn diagrams, region lists, and ORA analysis have been successfully generated.")

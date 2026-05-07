# mRNA_vs_protein_venn_diagram_and_ORA.R

# Load necessary libraries
library(dplyr)
library(openxlsx)
library(VennDiagram)
library(grid)
library(futile.logger)

# Suppress VennDiagram logger output to prevent cluttering the console
flog.threshold(ERROR, name = "VennDiagramLogger")

# Define paths
input_dir <- "C:/Users/danny/Documents/R_project/linkedomicskb_TP53_GOF/mRNA_DEG_vs_protein_DEG_meta"
output_dir <- file.path(input_dir, "venn_diagram")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read data
data_file <- file.path(input_dir, "META_DEG_TP53mt_vs_TP53wt.csv")
if (!file.exists(data_file)) {
  stop("Input file META_DEG_TP53mt_vs_TP53wt.csv not found.")
}
meta_data <- read.csv(data_file, stringsAsFactors = FALSE)

# Filter out rows without gene_symbol
meta_data <- meta_data %>% filter(!is.na(gene_symbol) & gene_symbol != "")

# ==============================================================================
# 1. Process UP-regulated genes
# ==============================================================================

# Extract lists for Venn diagram 1 (UP)
mRNA_up_list <- meta_data %>%
  filter(mRNA_Z_meta > 0 & mRNA_padj < 0.05) %>%
  pull(gene_symbol) %>%
  unique()

protein_up_list <- meta_data %>%
  filter(protein_Z_meta > 0 & protein_padj < 0.05) %>%
  pull(gene_symbol) %>%
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

up_regions <- get_venn_regions(mRNA_up_list, protein_up_list, "mRNA", "protein")
write.csv(up_regions, file.path(output_dir, "Venn_UP_regions.csv"), row.names = FALSE, na = "")
write.xlsx(up_regions, file.path(output_dir, "Venn_UP_regions.xlsx"), row.names = FALSE)

# Plot UP Venn Diagram
# Using NPG (Nature Publishing Group) colors for high-impact journal style
color_up <- c("#E64B35FF", "#4DBBD5FF") # Nature Red and Nature Blue

venn.plot.up <- venn.diagram(
  x = list(mRNA_UP = mRNA_up_list, Protein_UP = protein_up_list),
  filename = NULL,
  lwd = 2,
  lty = 'solid',
  col = color_up,
  fill = color_up,
  alpha = 0.6,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  cex = 1.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  scaled = FALSE, # Ensure both circles are the same size
  main = "Up-regulated Genes\n(mRNA vs Protein)",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.1
)

# Save UP Venn to PDF
pdf(file.path(output_dir, "Venn_UP.pdf"), width = 6, height = 6)
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# Save UP Venn to TIFF
tiff(file.path(output_dir, "Venn_UP.tiff"), width = 6, height = 6, units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# ==============================================================================
# 2. Process DOWN-regulated genes
# ==============================================================================

# Extract lists for Venn diagram 2 (DOWN)
mRNA_down_list <- meta_data %>%
  filter(mRNA_Z_meta < 0 & mRNA_padj < 0.05) %>%
  pull(gene_symbol) %>%
  unique()

protein_down_list <- meta_data %>%
  filter(protein_Z_meta < 0 & protein_padj < 0.05) %>%
  pull(gene_symbol) %>%
  unique()

down_regions <- get_venn_regions(mRNA_down_list, protein_down_list, "mRNA", "protein")
write.csv(down_regions, file.path(output_dir, "Venn_DOWN_regions.csv"), row.names = FALSE, na = "")
write.xlsx(down_regions, file.path(output_dir, "Venn_DOWN_regions.xlsx"), row.names = FALSE)

# Plot DOWN Venn Diagram
# Using another set of NPG colors
color_down <- c("#00A087FF", "#3C5488FF") # Nature Green and Nature Dark Blue

venn.plot.down <- venn.diagram(
  x = list(mRNA_DOWN = mRNA_down_list, Protein_DOWN = protein_down_list),
  filename = NULL,
  lwd = 2,
  lty = 'solid',
  col = color_down,
  fill = color_down,
  alpha = 0.6,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  cex = 1.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  scaled = FALSE, # Ensure both circles are the same size
  main = "Down-regulated Genes\n(mRNA vs Protein)",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.1
)

# Save DOWN Venn to PDF
pdf(file.path(output_dir, "Venn_DOWN.pdf"), width = 6, height = 6)
grid.newpage()
grid.draw(venn.plot.down)
dev.off()

# Save DOWN Venn to TIFF
tiff(file.path(output_dir, "Venn_DOWN.tiff"), width = 6, height = 6, units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.down)
dev.off()

print("Venn diagrams and region lists have been successfully generated.")

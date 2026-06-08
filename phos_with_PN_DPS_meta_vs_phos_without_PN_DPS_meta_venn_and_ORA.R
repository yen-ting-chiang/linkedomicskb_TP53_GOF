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
#   "gene_symbol_phosphosite" column. For ORA, the parent gene is extracted
#   (e.g., "MCM6_S762" -> "MCM6"), duplicates are removed, and the resulting
#   unique gene list is used for enrichment analysis.
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
#       ORA_{collection}_{region}.csv, .xlsx, _bubble.pdf, _bubble.tiff
#
# References:
#   - Liberzon A, et al. The Molecular Signatures Database (MSigDb) hallmark
#     gene set collection. Cell Syst. 2015;1(6):417-425. PMID: 26771021
#   - Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes.
#     Nucleic Acids Res. 2000;28(1):27-30. PMID: 10592173
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

# ==============================================================================
# USER CONFIGURATION: Venn Diagram Color Palettes
# ==============================================================================
# Select a color palette for UP and DOWN Venn diagrams by changing the index.
# Each palette contains two colors: c(left_circle, right_circle).
#
# Available palettes for UP-regulated Venn diagram:
#   1 = NPG (Nature Publishing Group) Red & Blue    : "#E64B35FF", "#4DBBD5FF"
#   2 = Lancet Red & Blue                           : "#ED0000FF", "#0099B4FF"
#   3 = JAMA Blue & Orange                          : "#374E55FF", "#DF8F44FF"
#   4 = NEJM Warm Red & Deep Cerulean               : "#BC3C29FF", "#0072B5FF"
#   5 = JCO Cerulean & Crimson                      : "#0073C2FF", "#EFC000FF"
#   6 = Science Orange & Teal                       : "#FF6F00FF", "#009688FF"
#   7 = Cell Coral & Steel Blue                     : "#FF6B6BFF", "#4682B4FF"
#   8 = PNAS Blue & Gold                            : "#3B4992FF", "#E6A024FF"
#
# Available palettes for DOWN-regulated Venn diagram:
#   1 = NPG Green & Dark Blue                       : "#00A087FF", "#3C5488FF"
#   2 = Lancet Green & Purple                       : "#00468BFF", "#925E9FFF"
#   3 = JAMA Teal & Mauve                           : "#45818EFF", "#96567DFF"
#   4 = NEJM Dark Cyan & Muted Purple               : "#20854EFF", "#7876B1FF"
#   5 = JCO Forest Green & Slate                    : "#20854EFF", "#6F99ADFF"
#   6 = Science Emerald & Indigo                    : "#2E7D32FF", "#4A148CFF"
#   7 = Cell Sage & Navy                            : "#78A878FF", "#2C3E50FF"
#   8 = PNAS Forest & Plum                          : "#008B45FF", "#7B3A96FF"

palette_up_choice   <- 1  # Change this number to select UP palette (1-8)
palette_down_choice <- 1  # Change this number to select DOWN palette (1-8)

up_palettes <- list(
  c("#E64B35FF", "#4DBBD5FF"),  # 1: NPG Red & Blue
  c("#ED0000FF", "#0099B4FF"),  # 2: Lancet Red & Blue
  c("#374E55FF", "#DF8F44FF"),  # 3: JAMA Blue & Orange
  c("#BC3C29FF", "#0072B5FF"),  # 4: NEJM Red & Cerulean
  c("#0073C2FF", "#EFC000FF"),  # 5: JCO Cerulean & Gold
  c("#FF6F00FF", "#009688FF"),  # 6: Science Orange & Teal
  c("#FF6B6BFF", "#4682B4FF"),  # 7: Cell Coral & Steel Blue
  c("#3B4992FF", "#E6A024FF")   # 8: PNAS Blue & Gold
)

down_palettes <- list(
  c("#00A087FF", "#3C5488FF"),  # 1: NPG Green & Dark Blue
  c("#00468BFF", "#925E9FFF"),  # 2: Lancet Blue & Purple
  c("#45818EFF", "#96567DFF"),  # 3: JAMA Teal & Mauve
  c("#20854EFF", "#7876B1FF"),  # 4: NEJM Cyan & Purple
  c("#20854EFF", "#6F99ADFF"),  # 5: JCO Forest & Slate
  c("#2E7D32FF", "#4A148CFF"),  # 6: Science Emerald & Indigo
  c("#78A878FF", "#2C3E50FF"),  # 7: Cell Sage & Navy
  c("#008B45FF", "#7B3A96FF")   # 8: PNAS Forest & Plum
)

color_up   <- up_palettes[[palette_up_choice]]
color_down <- down_palettes[[palette_down_choice]]

cat("Venn UP palette:", palette_up_choice, "->", color_up, "\n")
cat("Venn DOWN palette:", palette_down_choice, "->", color_down, "\n")

# ==============================================================================
# USER CONFIGURATION: ORA Gene Set Collections
# ==============================================================================
# Set each collection to TRUE or FALSE to include/exclude from ORA analysis.
# Each enabled collection will produce its own set of output files:
#   ORA_{collection_label}_{region}.csv, .xlsx, _bubble.pdf, _bubble.tiff
#
# Reference for MSigDB collections:
#   Subramanian A, et al. Gene set enrichment analysis: a knowledge-based
#   approach for interpreting genome-wide expression profiles. Proc Natl Acad
#   Sci U S A. 2005;102(43):15545-50. PMID: 16199517

run_Hallmark   <- TRUE   # MSigDB H: Hallmark gene sets (50 sets)
run_KEGG       <- TRUE   # MSigDB C2:CP:KEGG_LEGACY (186 KEGG pathway sets)

# Build the list of gene set collections to run
# Each entry: list(label, category, subcategory, prefix_to_strip)
#   - label: short name used in output filenames
#   - category: MSigDB category (e.g., "H", "C2")
#   - subcategory: MSigDB subcategory (NULL for none, e.g., "CP:KEGG_LEGACY")
#   - prefix_to_strip: regex to remove from pathway names for cleaner plotting
ora_collections <- list()

if (run_Hallmark) {
  ora_collections[["Hallmark"]] <- list(
    label           = "Hallmark",
    category        = "H",
    subcategory     = NULL,
    prefix_to_strip = "(?i)^hallmark[_ ]"
  )
}

if (run_KEGG) {
  ora_collections[["KEGG"]] <- list(
    label           = "KEGG",
    category        = "C2",
    subcategory     = "CP:KEGG_LEGACY",
    prefix_to_strip = "(?i)^kegg[_ ]"
  )
}

cat("ORA collections enabled:",
    paste(names(ora_collections), collapse = ", "), "\n\n")

# ==============================================================================
# USER CONFIGURATION: Bubble Plot Color Palettes
# ==============================================================================
# Select a color gradient for bubble plot p.adjust mapping.
# The gradient runs from low p.adjust (most significant) to high p.adjust.
# Change palette_bubble_choice to switch palette.
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

palette_bubble_choice <- 6  # Change this number to select bubble palette (1-8)

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

cat("Bubble plot palette:", palette_bubble_choice, "->",
    bubble_color_low, "(low) to", bubble_color_high, "(high)\n")

# ==============================================================================
# Setup: Paths and Data
# ==============================================================================

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

# Helper function: extract parent gene from gene_symbol_phosphosite
# e.g., "MCM6_S762" -> "MCM6", "TP53BP1_S1758" -> "TP53BP1"
extract_parent_gene <- function(gene_symbol_phosphosite) {
  sub("_[^_]*$", "", gene_symbol_phosphosite)
}

# Define ORA background universe: all parent genes detected in the integrated
# dataset (deduplicated). This ensures the Fisher's exact test denominator
# reflects the experimentally observable gene space, rather than defaulting to
# all genes present in MSigDB.
# Reference: Timmons JA, et al. Multiple sources of bias confound functional
#   enrichment analysis of global -omics data. Genome Biol. 2015;16:186.
#   PMID: 26346307
# Reference: Wijesooriya K, et al. Urgent need for consistent standards in
#   functional enrichment analysis. PLoS Comput Biol. 2022;18(3):e1009935.
#   PMID: 35263338
ora_universe <- unique(extract_parent_gene(meta_data$gene_symbol_phosphosite))
ora_universe <- ora_universe[!is.na(ora_universe) & ora_universe != ""]
cat("ORA background universe:", length(ora_universe), "unique parent genes\n")

# ==============================================================================
# Load MSigDB Gene Set Collections
# ==============================================================================

msig_db_list <- list()

for (coll_name in names(ora_collections)) {
  coll <- ora_collections[[coll_name]]
  cat("Loading MSigDB collection:", coll_name, "\n")
  
  if (is.null(coll$subcategory)) {
    msig_db_list[[coll_name]] <- msigdbr(species = "Homo sapiens",
                                          category = coll$category) %>%
      dplyr::select(gs_name, gene_symbol)
  } else {
    msig_db_list[[coll_name]] <- msigdbr(species = "Homo sapiens",
                                          category = coll$category,
                                          subcategory = coll$subcategory) %>%
      dplyr::select(gs_name, gene_symbol)
  }
  
  cat("  Loaded", length(unique(msig_db_list[[coll_name]]$gs_name)),
      "gene sets\n")
}

# ==============================================================================
# ORA Function (supports multiple collections)
# ==============================================================================

# Function to perform ORA and plot results for a single collection
# Input: phosphosite_list = vector of gene_symbol_phosphosite values
#        The parent genes are extracted and deduplicated before ORA
perform_ora_and_plot <- function(phosphosite_list, region_name, out_dir,
                                 universe, term2gene, coll_label,
                                 prefix_to_strip) {
  phosphosite_list <- phosphosite_list[!is.na(phosphosite_list) & phosphosite_list != ""]
  if (length(phosphosite_list) == 0) return(NULL)
  
  # Extract parent genes and remove duplicates
  gene_list <- unique(extract_parent_gene(phosphosite_list))
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list) == 0) return(NULL)
  
  cat("  ORA [", coll_label, "] for", region_name, ":",
      length(phosphosite_list), "phosphosites ->",
      length(gene_list), "unique parent genes\n")
  
  # Perform ORA using clusterProfiler
  # Set cutoffs to 1 to ensure we can pad up to 5 even if non-significant
  ora_res <- enricher(gene = gene_list, TERM2GENE = term2gene,
                      universe = universe,
                      pvalueCutoff = 1, qvalueCutoff = 1)
  
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
      
      file_prefix <- paste0("ORA_", coll_label, "_", region_name)
      
      write.csv(write_df, file.path(out_dir, paste0(file_prefix, ".csv")),
                row.names = FALSE)
      write.xlsx(write_df, file.path(out_dir, paste0(file_prefix, ".xlsx")),
                 row.names = FALSE)
      
      # Clean pathway names: remove collection prefix for cleaner plotting
      plot_df$Description <- gsub(prefix_to_strip, "", plot_df$Description)
      
      # Order Description factor by p.adjust (descending so smallest padj is at the top)
      plot_df$Description <- factor(plot_df$Description,
                                     levels = rev(plot_df$Description))
      
      # Control plot title visibility (TRUE to show, FALSE to hide)
      show_plot_title <- TRUE
      plot_title_text <- if(show_plot_title) {
        paste(coll_label, "ORA:", region_name)
      } else {
        NULL
      }
      
      # Bubble plot (custom ggplot)
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
      
      if(!show_plot_title) {
        p <- p + theme(plot.title = element_blank())
      }
      
      pdf(file.path(out_dir, paste0(file_prefix, "_bubble.pdf")),
          width = 3.8, height = 2.5)
      print(p)
      dev.off()
      
      tiff(file.path(out_dir, paste0(file_prefix, "_bubble.tiff")),
           width = 3.8, height = 2.5, units = "in", res = 300,
           compression = "lzw")
      print(p)
      dev.off()
      
      cat("  ORA results and bubble plots generated for [",
          coll_label, "]", region_name, "\n")
    }
  } else {
    cat("  No significant ORA results for [", coll_label, "]",
        region_name, "\n")
  }
}

# Wrapper: run ORA across all enabled collections for a given region
run_all_ora <- function(phosphosite_list, region_name, out_dir, universe) {
  for (coll_name in names(ora_collections)) {
    coll <- ora_collections[[coll_name]]
    perform_ora_and_plot(
      phosphosite_list = phosphosite_list,
      region_name      = region_name,
      out_dir          = out_dir,
      universe         = universe,
      term2gene        = msig_db_list[[coll_name]],
      coll_label       = coll$label,
      prefix_to_strip  = coll$prefix_to_strip
    )
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
  colnames(res) <- c(paste0("Only_", nameA), "Intersection",
                      paste0("Only_", nameB))
  return(res)
}

# Helper function: convert phosphosite Venn regions to parent gene regions
# Extracts the parent gene from each gene_symbol_phosphosite in each column,
# removes duplicates, and pads columns to equal length.
get_parent_gene_regions <- function(venn_regions_df) {
  gene_cols <- lapply(venn_regions_df, function(col) {
    col <- col[!is.na(col) & col != ""]
    unique(extract_parent_gene(col))
  })
  max_len <- max(sapply(gene_cols, length))
  gene_cols <- lapply(gene_cols, function(x) {
    length(x) <- max_len
    x
  })
  res <- as.data.frame(gene_cols, stringsAsFactors = FALSE)
  colnames(res) <- colnames(venn_regions_df)
  return(res)
}

up_regions <- get_venn_regions(phos_without_PN_up_list, phos_with_PN_up_list,
                               "phos_without_PN", "phos_with_PN")
write.csv(up_regions, file.path(output_dir, "Venn_UP_regions.csv"),
          row.names = FALSE, na = "")
write.xlsx(up_regions, file.path(output_dir, "Venn_UP_regions.xlsx"),
           row.names = FALSE)

# Save parent gene (deduplicated) version of UP regions
up_regions_genes <- get_parent_gene_regions(up_regions)
write.csv(up_regions_genes, file.path(output_dir, "Venn_UP_regions_parent_gene.csv"),
          row.names = FALSE, na = "")
write.xlsx(up_regions_genes, file.path(output_dir, "Venn_UP_regions_parent_gene.xlsx"),
           row.names = FALSE)

# Plot UP Venn Diagram
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
  scaled = FALSE,
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
tiff(file.path(output_dir, "Venn_UP.tiff"), width = 1.5, height = 1.5,
     units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# ==============================================================================
# 1.5. ORA for UP-regulated regions (using parent genes, deduplicated)
# ==============================================================================
run_all_ora(up_regions$Only_phos_without_PN, "UP_Only_phos_without_PN",
            output_dir, ora_universe)
run_all_ora(up_regions$Intersection, "UP_Intersection",
            output_dir, ora_universe)
run_all_ora(up_regions$Only_phos_with_PN, "UP_Only_phos_with_PN",
            output_dir, ora_universe)

# ORA for union of two adjacent Venn regions (full circles)
# Left circle:  Only_phos_without_PN + Intersection = all phos_without_PN significant
# Right circle: Only_phos_with_PN   + Intersection = all phos_with_PN significant
up_union_without_PN <- unique(c(
  up_regions$Only_phos_without_PN[!is.na(up_regions$Only_phos_without_PN)],
  up_regions$Intersection[!is.na(up_regions$Intersection)]
))
up_union_with_PN <- unique(c(
  up_regions$Only_phos_with_PN[!is.na(up_regions$Only_phos_with_PN)],
  up_regions$Intersection[!is.na(up_regions$Intersection)]
))

run_all_ora(up_union_without_PN, "UP_Union_phos_without_PN",
            output_dir, ora_universe)
run_all_ora(up_union_with_PN, "UP_Union_phos_with_PN",
            output_dir, ora_universe)

# ORA for total union of all three Venn regions
# (all phosphosites significant in at least one method)
up_union_all <- unique(c(
  up_regions$Only_phos_without_PN[!is.na(up_regions$Only_phos_without_PN)],
  up_regions$Intersection[!is.na(up_regions$Intersection)],
  up_regions$Only_phos_with_PN[!is.na(up_regions$Only_phos_with_PN)]
))

run_all_ora(up_union_all, "UP_Union_all",
            output_dir, ora_universe)

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
write.csv(down_regions, file.path(output_dir, "Venn_DOWN_regions.csv"),
          row.names = FALSE, na = "")
write.xlsx(down_regions, file.path(output_dir, "Venn_DOWN_regions.xlsx"),
           row.names = FALSE)

# Save parent gene (deduplicated) version of DOWN regions
down_regions_genes <- get_parent_gene_regions(down_regions)
write.csv(down_regions_genes, file.path(output_dir, "Venn_DOWN_regions_parent_gene.csv"),
          row.names = FALSE, na = "")
write.xlsx(down_regions_genes, file.path(output_dir, "Venn_DOWN_regions_parent_gene.xlsx"),
           row.names = FALSE)

# Plot DOWN Venn Diagram
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
  scaled = FALSE,
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
tiff(file.path(output_dir, "Venn_DOWN.tiff"), width = 1.5, height = 1.5,
     units = "in", res = 300, compression = "lzw")
grid.newpage()
grid.draw(venn.plot.down)
dev.off()

# ==============================================================================
# 2.5. ORA for DOWN-regulated regions (using parent genes, deduplicated)
# ==============================================================================
run_all_ora(down_regions$Only_phos_without_PN, "DOWN_Only_phos_without_PN",
            output_dir, ora_universe)
run_all_ora(down_regions$Intersection, "DOWN_Intersection",
            output_dir, ora_universe)
run_all_ora(down_regions$Only_phos_with_PN, "DOWN_Only_phos_with_PN",
            output_dir, ora_universe)

# ORA for union of two adjacent Venn regions (full circles)
down_union_without_PN <- unique(c(
  down_regions$Only_phos_without_PN[!is.na(down_regions$Only_phos_without_PN)],
  down_regions$Intersection[!is.na(down_regions$Intersection)]
))
down_union_with_PN <- unique(c(
  down_regions$Only_phos_with_PN[!is.na(down_regions$Only_phos_with_PN)],
  down_regions$Intersection[!is.na(down_regions$Intersection)]
))

run_all_ora(down_union_without_PN, "DOWN_Union_phos_without_PN",
            output_dir, ora_universe)
run_all_ora(down_union_with_PN, "DOWN_Union_phos_with_PN",
            output_dir, ora_universe)

# ORA for total union of all three Venn regions
down_union_all <- unique(c(
  down_regions$Only_phos_without_PN[!is.na(down_regions$Only_phos_without_PN)],
  down_regions$Intersection[!is.na(down_regions$Intersection)],
  down_regions$Only_phos_with_PN[!is.na(down_regions$Only_phos_with_PN)]
))

run_all_ora(down_union_all, "DOWN_Union_all",
            output_dir, ora_universe)

print("Venn diagrams, region lists, and ORA analysis have been successfully generated.")

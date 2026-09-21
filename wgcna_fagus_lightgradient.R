# =============================================================================
# WGCNA — Weighted Gene Co-expression Network Analysis
# Fagus sylvatica transcriptome — leaf-position-responsive genes
#
# Input:  15,126 leaf-position-responsive genes (from GAMM)
#         VST-normalised expression matrix
#         Sample metadata (physiological traits + drought group)
#
# Output: Gene modules, module-trait correlations, hub genes
#
# Pipeline:
#   1. Prepare expression matrix
#   2. Check sample quality (outlier detection)
#   3. Choose soft-thresholding power (scale-free topology)
#   4. Build weighted co-expression network
#   5. Detect modules
#   6. Module-trait correlation
#   7. Hub gene identification
#   8. Export for visualisation
#
# Data: ENA accession PRJEB64934
#       Place all required input files in the working directory before running.
#       See README.md for the full list of input files.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------

# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("WGCNA")

library(WGCNA)
library(dplyr)
library(ggplot2)
library(pheatmap)

# WGCNA requires multi-threading to be enabled
allowWGCNAThreads()

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
message("--- Step 1: Load data ---")

# Load VST expression matrix
expr_full <- as.matrix(
  read.table("expr_full_vst.txt", header = TRUE, row.names = 1)
)

# Load GAMM results
gamm_results <- read.csv("gamm_leafpos_results.csv")

# Load sample metadata
sample_info_clean <- read.table(
  "sample_info_clean.txt",
  header    = TRUE,
  row.names = 1
)

# Load response shape gene lists
# leafpos_responsive_genes.txt:
#   - Predictor: leaf position ONLY (ordered factor 1-7)
#   - Controls: kinship PC1 + PC2 (fixed), genotype (random)
#   - Identifies: genes with developmental positional gradient
#   - Does NOT identify: genes responding specifically to Lux,
#     PIabs or Fv/Fm independently of position
#   - N = 15,126 genes
leafpos_genes <- scan("leafpos_responsive_genes.txt",
                       what = character(), quiet = TRUE)

message("Total VST genes:            ", nrow(expr_full))
message("Leaf-position responsive:   ", length(leafpos_genes))

# Total VST genes:            43,649
# Leaf-position responsive:   15,126

# -----------------------------------------------------------------------------
# 2. Prepare expression matrix for WGCNA
#    WGCNA expects: rows = samples, columns = genes
#    Use all 15,126 responsive genes
# -----------------------------------------------------------------------------
message("\n--- Step 2: Prepare WGCNA expression matrix ---")

# Subset to responsive genes and align samples
common_samples <- intersect(colnames(expr_full),
                             rownames(sample_info_clean))
expr_wgcna <- t(expr_full[leafpos_genes, common_samples])
# Now: rows = samples (185), columns = genes (15,126)

message("WGCNA matrix: ", nrow(expr_wgcna), " samples × ",
        ncol(expr_wgcna), " genes")

# Align metadata
meta_wgcna <- sample_info_clean[common_samples, ]
stopifnot(all(rownames(meta_wgcna) == rownames(expr_wgcna)))
message("Sample alignment: OK")

# -----------------------------------------------------------------------------
# 3. Sample quality check — outlier detection
#    Samples with very different expression profiles may be technical outliers
# -----------------------------------------------------------------------------
message("\n--- Step 3: Sample quality check ---")

# Hierarchical clustering of samples
sample_tree <- hclust(dist(expr_wgcna), method = "average")

# Plot sample dendrogram
pdf("wgcna_sample_dendrogram.pdf", width = 14, height = 6)
plot(sample_tree,
     main  = "Sample clustering — check for outliers",
     sub   = "",
     xlab  = "",
     cex   = 0.7,
     cex.main = 1)
abline(h = 80, col = "red", lty = 2)   # adjust height if needed
dev.off()
message("Saved: wgcna_sample_dendrogram.pdf")
message("Check the dendrogram for outlier samples before proceeding")

# Check for genes or samples with too many missing values
gsg <- goodSamplesGenes(expr_wgcna, verbose = 3)
if (!gsg$allOK) {
  message("Removing flagged genes/samples...")
  expr_wgcna <- expr_wgcna[gsg$goodSamples, gsg$goodGenes]
  message("Retained: ", nrow(expr_wgcna), " samples × ",
          ncol(expr_wgcna), " genes")
} else {
  message("All samples and genes passed quality check.")
}

# Check wgcna_sample_dendrogram.pdf before proceeding:
# Any sample that branches off much earlier than the others is a potential outlier.

# B9_7 and B7_7 branch off at height ~120, clearly separated from the main cluster.
# B1_6 branches at a lower height and is closer to the main group —> it is a borderline case and should be retained.
# check what makes the other two unusual:

# Check their expression profiles
# Are they outliers in terms of total expression?
sample_totals <- rowSums(expr_wgcna)
sort(sample_totals)[1:10]   # lowest total expression samples

# Check their metadata — are they from unusual positions or trees?
sample_info_clean[c("B9_7", "B7_7", "B1_6"), ]

# Check their distance to the nearest cluster member
dist_mat <- as.matrix(dist(expr_wgcna))
# Mean distance to all other samples
mean_dist <- rowMeans(dist_mat)
sort(mean_dist, decreasing = TRUE)[1:5]

#     R6_4     R8_6     B8_7     B7_7     B7_6     B9_3     R7_6     B5_2
# 153007.9 153100.7 153119.7 153133.0 153148.2 153176.0 153184.2 153194.8
#     R6_3     B5_3
# 153196.3 153246.6
#      genotype condition drought_reaction   lux  rel_expo mean_pi_abs mean_fvfm
# B9_7       B9         7               nd 20000 1.0000000      0.7725    0.6625
# B7_7       B7         7               nd 20000 1.0000000      0.3435    0.5860
# B1_6       B1         6              res  1205 0.3362165      4.4025    0.8060
#     B9_7     B7_7     C4_7     B4_4     C4_6
# 138.4234 121.6127 113.6291 111.2721 110.9658

# Lux = 20,000 for both B9_7 and B7_7 -> sensor was saturated, actual light intensity might have been even higher
# Fv/Fm + PIabs is very low -> severe photosynthetic stress
# Mean distance is the lowest -> they have genuinely unusual transcriptomes -> remove
# B1_6 has Lux = 1,205 (normal), Fv/Fm = 0.806 (healthy), and a mean distance of 110.9 — only slightly above average.
# It is not a true outlier, just slightly more variable -> Keep

# Remove B9_7 and B7_7
samples_to_remove <- c("B9_7", "B7_7")
expr_wgcna <- expr_wgcna[!rownames(expr_wgcna) %in% samples_to_remove, ]
meta_wgcna <- meta_wgcna[!rownames(meta_wgcna) %in% samples_to_remove, ]

message("Samples after outlier removal: ", nrow(expr_wgcna))
# Should be 183

stopifnot(all(rownames(expr_wgcna) == rownames(meta_wgcna)))
message("Alignment OK — proceed to Step 4 (soft-thresholding power)")

# -----------------------------------------------------------------------------
# 4. Soft-thresholding power selection
#    WGCNA transforms the correlation matrix to a weighted adjacency matrix
#    using a power function. The power is chosen so the network approximates
#    scale-free topology (characteristic of biological networks).
#    Choose the lowest power where R² > 0.85
# -----------------------------------------------------------------------------
message("\n--- Step 4: Soft-thresholding power selection ---")

powers     <- c(1:10, seq(12, 30, by = 2))
sft        <- pickSoftThreshold(
  expr_wgcna,
  powerVector  = powers,
  verbose      = 5,
  networkType  = "signed"    # signed: preserves direction of correlation
)

# Plot scale-free topology fit
pdf("wgcna_softthreshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# R² plot
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)",
     ylab = "Scale-free topology R²",
     main = "Scale-free topology fit",
     type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red", cex = 0.8)
abline(h = 0.85, col = "blue", lty = 2)

# Mean connectivity plot
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity",
     main = "Mean connectivity",
     type = "n")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, col = "red", cex = 0.8)
dev.off()
message("Saved: wgcna_softthreshold.pdf")

# Automatically select power where R² first exceeds 0.85
r2_values  <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
soft_power <- powers[which(r2_values > 0.85)[1]]

if (is.na(soft_power)) {
  message("WARNING: R² never exceeded 0.85 — using power = 12 as fallback")
  soft_power <- 12
}
message("Selected soft-thresholding power: ", soft_power)
# Selected soft-thresholding power: 9

# Check wgcna_softthreshold.pdf -> You want the lowest power where the R² line (left plot) crosses 0.85 (blue dashed line).
# For plant transcriptomes with signed networks this is typically between 12 and 20.
# The script selects it automatically but check the plot matches the selected value -> decided on power 9

# -----------------------------------------------------------------------------
# 5. Build network and detect modules
#    blockwiseModules handles large gene sets by splitting into blocks
#    minModuleSize = 30 is standard for transcriptome data
#    mergeCutHeight = 0.25 merges modules with >75% similarity
# -----------------------------------------------------------------------------

### This is the slowest step — with 15,126 genes and 185 samples expect 20–40 minutes.
### The maxBlockSize = 20000 setting means it will run in one block
### If it runs out of memory reduce to maxBlockSize = 10,000.

message("\n--- Step 5: Build network and detect modules ---")
message("This is the most computationally intensive step...")

bwnet <- blockwiseModules(
  expr_wgcna,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,    # use colour labels
  pamRespectsDendro = FALSE,
  saveTOMs          = TRUE,
  saveTOMFileBase   = "wgcna_TOM",
  verbose           = 3,
  maxBlockSize      = 20000,    # adjust if memory is limited
  nThreads          = 4         # parallel threads
)

# Save network object
saveRDS(bwnet, "wgcna_network.rds")
message("Network saved: wgcna_network.rds")

# Module summary
module_colours <- bwnet$colors
module_table   <- table(module_colours)
n_modules      <- sum(names(module_table) != "grey")

message("\nModule summary:")
message("  Total modules detected: ", n_modules,
        " (grey = unassigned genes)")
message("  Unassigned genes (grey): ",
        sum(module_colours == "grey"))
print(sort(module_table, decreasing = TRUE))

# Plot module dendrogram
pdf("wgcna_module_dendrogram.pdf", width = 14, height = 6)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  module_colours[bwnet$blockGenes[[1]]],
  "Module colours",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene dendrogram and module colours"
)
dev.off()
message("Saved: wgcna_module_dendrogram.pdf")

# Module summary:
#   Total modules detected: 9 (grey = unassigned genes)
#   Unassigned genes (grey): 3390
# module_colours
# turquoise      grey      blue     brown    yellow     green       red     black
#      4903      3390      2006      1392       978       918       851       246
#      pink   magenta
#       238       204

# The network is dominated by two large modules. The low module count (9) suggests mergeCutHeight = 0.25 may be merging modules that should stay separate.
# Check the module eigengene correlation: if turquoise and blue are highly correlated with each other, they should probably be kept separate (lower mergeCutHeight).
# If they are anticorrelated or uncorrelated, the current structure is fine:

# -----------------------------------------------------------------------------
# 6. Module eigengenes
#    The module eigengene (ME) is the first PC of all genes in a module
#    It summarises the overall expression profile of the module
# -----------------------------------------------------------------------------
# Correlation between module eigengenes
MEs <- moduleEigengenes(expr_wgcna, bwnet$colors)$eigengenes
MEs <- orderMEs(MEs)

# Correlation between the two largest modules
cor(MEs$MEturquoise, MEs$MEblue)
# If > 0.7: they are similar → lower mergeCutHeight
# If < 0.5: they are distinct → current structure is fine

# Plot all ME correlations
pdf("wgcna_ME_correlation.pdf", width=7, height=6)
plotEigengeneNetworks(MEs, "",
                      marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2))
dev.off()

# -0.7075225

# the two large modules represent opposing transcriptional programmes along the light gradient,
# not two halves of the same response.
# One module likely increases expression from shade to sun while the other decreases.
# They should absolutely remain as separate modules — do not lower mergeCutHeight to split them further,
# the current structure is biologically meaningful.

# ME correlation dendrogram:
# From the plot the modules cluster into two clear groups:
# Positively correlated cluster on the left -> likely represent genes that respond in the same direction to the light gradient
# Second cluster on the right and are anticorrelated with the first cluster -> opposing transcriptional programme.
# => some genes are induced by high light (shade-to-sun upregulation) while others are repressed (sun-to-shade upregulation)

# -----------------------------------------------------------------------------
# 7. Module-trait correlation
#    Correlate each module eigengene with physiological and
#    environmental traits
#    Identifies which modules respond to which aspect of the light gradient
# -----------------------------------------------------------------------------
message("\n--- Step 7: Module-trait correlation ---")

# Prepare trait matrix
# Use numeric traits only — convert drought_reaction to binary
traits <- data.frame(
  condition       = as.numeric(meta_wgcna$condition),
  lux             = as.numeric(meta_wgcna$lux),
  rel_expo        = as.numeric(meta_wgcna$rel_expo),
  mean_pi_abs     = as.numeric(meta_wgcna$mean_pi_abs),
  mean_fvfm       = as.numeric(meta_wgcna$mean_fvfm),
  drought_res     = as.integer(meta_wgcna$drought_reaction == "res"),
  drought_sus     = as.integer(meta_wgcna$drought_reaction == "sus"),
  row.names       = rownames(meta_wgcna)
)

# Remove drought columns if all NA (only 8 trees have drought info)
# Replace NA with 0 for binary drought traits
traits$drought_res[is.na(traits$drought_res)] <- 0
traits$drought_sus[is.na(traits$drought_sus)] <- 0

message("Trait matrix: ", nrow(traits), " samples × ",
        ncol(traits), " traits")

# Correlation between module eigengenes and traits
n_samples_mt   <- nrow(MEs)
module_trait_r <- cor(MEs, traits, use = "pairwise.complete.obs")
module_trait_p <- corPvalueStudent(module_trait_r, n_samples_mt)

# Adjust p-values for multiple testing
module_trait_p_adj <- matrix(
  p.adjust(as.vector(module_trait_p), method = "BH"),
  nrow = nrow(module_trait_p),
  ncol = ncol(module_trait_p),
  dimnames = dimnames(module_trait_p)
)

# Create text matrix for heatmap (r and p combined)
text_matrix <- paste0(
  round(module_trait_r, 2), "\n(",
  round(module_trait_p_adj, 3), ")"
)
dim(text_matrix) <- dim(module_trait_r)

# Module-trait heatmap
pdf("wgcna_module_trait_heatmap.pdf",
    width  = 10,
    height = max(6, nrow(module_trait_r) * 0.35))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix       = module_trait_r,
  xLabels      = colnames(traits),
  yLabels      = rownames(module_trait_r),
  ySymbols     = rownames(module_trait_r),
  colorMatrix  = module_trait_p_adj < 0.05,   # highlight significant
  colors       = blueWhiteRed(50),
  textMatrix   = text_matrix,
  setStdMargins = FALSE,
  cex.text     = 0.5,
  zlim         = c(-1, 1),
  main         = "Module-trait correlations\n(r value, FDR-adjusted p)"
)
dev.off()
message("Saved: wgcna_module_trait_heatmap.pdf")

# Save correlation tables
write.csv(as.data.frame(module_trait_r),
          "wgcna_module_trait_r.csv", row.names = TRUE)
write.csv(as.data.frame(module_trait_p_adj),
          "wgcna_module_trait_p_adj.csv", row.names = TRUE)

# Print top module-trait associations
message("\nTop module-trait associations (FDR < 0.05):")
mt_df <- as.data.frame(as.table(module_trait_r)) %>%
  rename(module = Var1, trait = Var2, r = Freq) %>%
  mutate(p_adj = as.vector(module_trait_p_adj)) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)
print(head(mt_df, 20))

# produced an error, try more robust conversion:
mt_long <- data.frame(
  module = rownames(module_trait_r)[row(module_trait_r)],
  trait  = colnames(module_trait_r)[col(module_trait_r)],
  r      = as.vector(module_trait_r),
  p_adj  = as.vector(module_trait_p_adj)
) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

message("Strongest module-trait associations (FDR < 0.05):")
write.csv(as.data.frame(mt_long),
          "mt_long.csv", row.names = TRUE)


# -----------------------------------------------------------------------------
# 8. Add response shape annotation to modules
#    Check whether modules are enriched for linear vs non-linear genes
# -----------------------------------------------------------------------------
message("\n--- Step 8: Response shape per module ---")

shape_df <- gamm_results %>%
  dplyr::select(gene_id, response_shape, edf, delta_aic) %>%
  filter(gene_id %in% names(module_colours))

module_shape_df <- data.frame(
  gene_id = names(module_colours),
  module  = module_colours
) %>%
  left_join(shape_df, by = "gene_id")

# Shape distribution per module
shape_summary <- module_shape_df %>%
  filter(module != "grey") %>%
  group_by(module, response_shape) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(module, desc(n))

message("Response shape distribution per module (top modules):")
shape_summary %>%
  filter(module %in% names(sort(module_table,
                                 decreasing = TRUE)[1:5])) %>%
  print()

write.csv(shape_summary, "wgcna_module_shape_summary.csv",
          row.names = FALSE)

#    module    response_shape       n   pct
#    <chr>     <chr>            <int> <dbl>
#  1 blue      linear            1096  54.6
#  2 blue      weakly_nonlinear   800  39.9
#  3 blue      nonlinear          110   5.5
#  4 brown     linear            1003  72.1
#  5 brown     weakly_nonlinear   357  25.6
#  6 brown     nonlinear           32   2.3
#  7 turquoise linear            3542  72.2
#  8 turquoise weakly_nonlinear  1251  25.5
#  9 turquoise nonlinear          110   2.2
# 10 yellow    weakly_nonlinear   852  87.1
# 11 yellow    linear              98  10
# 12 yellow    nonlinear           28   2.9

# -----------------------------------------------------------------------------
# 9. Hub gene identification
#    Hub genes have the highest module membership (kME)
#    kME = correlation between gene expression and module eigengene
#    High kME = gene is central to the module's co-expression pattern
# -----------------------------------------------------------------------------
library(topGO)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load saved outputs — no need to rerun WGCNA
# -----------------------------------------------------------------------------
gene_module_df <- read.csv("wgcna_gene_module_assignments.csv")

# Load topGO annotation
# File: Bhaga_vs_Interpro5.69-101_GO_clean_TopGOinput
# Source: Fagus sylvatica Bhaga reference genome annotation (InterPro 5.69–101)
# Place this file in the working directory before running this section.
# See README.md for details on obtaining this file.
topgo <- read.delim(
  "Bhaga_vs_Interpro5.69-101_GO_clean_TopGOinput",
  header           = FALSE,
  stringsAsFactors = FALSE
)
colnames(topgo) <- c("gene_id", "GO_terms")
topgo$GO_terms  <- gsub("\\s+", "", topgo$GO_terms)

# Resolve dplyr/topGO namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# -----------------------------------------------------------------------------
# 2. Define gene universe
#    Use all genes tested in WGCNA (including grey)
# -----------------------------------------------------------------------------
gene_universe <- as.character(gene_module_df$gene_id)
message("Gene universe: ", length(gene_universe))

# Build gene-to-GO mapping
geneID2GO <- with(topgo,
                  setNames(strsplit(GO_terms, ","), gene_id))
geneID2GO <- geneID2GO[names(geneID2GO) %in% gene_universe]

pct_annotated <- round(length(geneID2GO) / length(gene_universe) * 100, 1)
message("Genes with GO annotation: ", length(geneID2GO),
        " (", pct_annotated, "%)")

# -----------------------------------------------------------------------------
# 3. TopGO function — runs for one module, one ontology
# -----------------------------------------------------------------------------
run_topgo_module <- function(module_name, gene_module_df,
                              gene_universe, geneID2GO,
                              ontology  = "BP",
                              node_size = 5,
                              top_nodes = 50) {

  message("  Running: ", module_name, " (", ontology, ")")

  # Genes in this module
  sig_genes <- as.character(
    gene_module_df$gene_id[gene_module_df$module == module_name]
  )

  sig_annotated <- sum(sig_genes %in% names(geneID2GO))
  if (sig_annotated < 5) {
    message("    Too few annotated genes (", sig_annotated, ") — skipping")
    return(NULL)
  }

  # Named factor: 1 = in module, 0 = not in module
  geneList        <- factor(as.integer(gene_universe %in% sig_genes))
  names(geneList) <- gene_universe

  # Build topGOdata object
  GOdata <- new("topGOdata",
                ontology  = ontology,
                allGenes  = geneList,
                annot     = annFUN.gene2GO,
                gene2GO   = geneID2GO,
                nodeSize  = node_size)

  # Run classic and elim algorithms
  result_classic <- runTest(GOdata, algorithm = "classic",
                            statistic = "fisher")
  result_elim    <- runTest(GOdata, algorithm = "elim",
                            statistic = "fisher")

  n_terms  <- length(score(result_classic))
  go_table <- GenTable(
    GOdata,
    classic_p = result_classic,
    elim_p    = result_elim,
    orderBy   = "elim_p",
    topNodes  = min(top_nodes, n_terms)
  )

  go_table$classic_p <- as.numeric(go_table$classic_p)
  go_table$elim_p    <- as.numeric(go_table$elim_p)
  go_table$elim_fdr  <- p.adjust(go_table$elim_p, method = "BH")
  go_table$module    <- module_name
  go_table$ontology  <- ontology
  go_table$n_module_genes <- length(sig_genes)

  return(go_table)
}

# -----------------------------------------------------------------------------
# 4. Run for all modules × three ontologies
#    Skip grey (unassigned genes)
# -----------------------------------------------------------------------------
modules    <- unique(gene_module_df$module)
modules    <- modules[modules != "grey"]
ontologies <- c("BP", "MF", "CC")

message("\nRunning GO enrichment for ", length(modules),
        " modules × 3 ontologies...")

all_results <- list()

for (mod in modules) {
  for (ont in ontologies) {
    result <- run_topgo_module(
      module_name    = mod,
      gene_module_df = gene_module_df,
      gene_universe  = gene_universe,
      geneID2GO      = geneID2GO,
      ontology       = ont
    )
    if (!is.null(result)) {
      all_results[[paste0(mod, "_", ont)]] <- result
    }
  }
}

# -----------------------------------------------------------------------------
# 5. Combine and save
# -----------------------------------------------------------------------------
results_df <- bind_rows(all_results)

# Significant results
sig_results <- results_df %>%
  filter(elim_p < 0.05) %>%
  arrange(module, ontology, elim_p)

message("\nTotal GO terms tested: ", nrow(results_df))
message("Significant (elim p < 0.05): ", nrow(sig_results))

write.csv(results_df,  "wgcna_topGO_all_results.csv",  row.names = FALSE)
write.csv(sig_results, "wgcna_topGO_significant.csv",  row.names = FALSE)
message("Saved: wgcna_topGO_all_results.csv")
message("Saved: wgcna_topGO_significant.csv")

# -----------------------------------------------------------------------------
# 6. Summary — top 5 BP terms per module
# -----------------------------------------------------------------------------
message("\n--- Top 5 BP terms per module (elim p < 0.05) ---")
sig_results %>%
  filter(ontology == "BP") %>%
  group_by(module) %>%
  slice_head(n = 5) %>%
  dplyr::select(module, GO.ID, Term, Significant,
                Expected, elim_p, elim_fdr) %>%
  print(n = 100)
# -----------------------------------------------------------------------------
# 10. Output summary
# -----------------------------------------------------------------------------
message("\n============================================================")
message("WGCNA complete. Output files:")
message("  wgcna_sample_dendrogram.pdf      — sample quality check")
message("  wgcna_softthreshold.pdf          — power selection plot")
message("  wgcna_module_dendrogram.pdf      — gene dendrogram + modules")
message("  wgcna_module_trait_heatmap.pdf   — module-trait correlations")
message("  wgcna_network.rds                — full network object")
message("  wgcna_module_trait_r.csv         — correlation matrix")
message("  wgcna_module_trait_p_adj.csv     — FDR-adjusted p-values")
message("  wgcna_module_shape_summary.csv   — EDF shape per module")
message("  wgcna_hub_genes.csv              — top hub genes per module")
message("  wgcna_gene_module_assignments.csv — gene → module table")
message("============================================================")

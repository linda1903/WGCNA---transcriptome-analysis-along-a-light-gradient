# WGCNA — Weighted Gene Co-expression Network Analysis  
## *Fagus sylvatica* transcriptome along a natural light gradient

This repository contains the R pipeline used to construct and interpret a
weighted gene co-expression network from RNA-seq data of European beech
(*Fagus sylvatica*) leaves sampled across a natural light gradient.
The analysis is part of a doctoral project at the Senckenberg Biodiversity
and Climate Research Centre (SBiK-F), Frankfurt.

**Manuscript in preparation.** Code and data are shared to support
transparency and reproducibility.

---

## Data availability

Raw sequencing reads and sample metadata are publicly available on the
European Nucleotide Archive (ENA) under accession
**[PRJEB64934](https://www.ebi.ac.uk/ena/browser/view/PRJEB64934)**.

The dataset comprises 185 mRNA-seq samples from 28 *F. sylvatica* trees
sampled at seven canopy positions (leaf positions 1–7, shade to full sun)
along a natural light gradient, together with paired photosynthetic
performance measurements (PI~abs~, Fv/Fm) and light intensity (lux,
relative exposure).

---

## Analysis overview

Starting from 15,126 leaf-position-responsive genes (identified by
Generalised Additive Mixed Models, GAMM), the pipeline:

1. Prepares the VST-normalised expression matrix (185 samples × 15,126 genes)
2. Detects and removes two technical outliers (B9_7, B7_7; saturated light
   sensors, severely impaired photosynthesis)
3. Selects soft-thresholding power for scale-free topology (power = 9)
4. Builds a signed weighted co-expression network using `blockwiseModules`
   (minModuleSize = 30, mergeCutHeight = 0.25)
5. Detects **9 co-expression modules** (11,736 assigned genes)
6. Correlates module eigengenes with physiological and environmental traits
7. Annotates modules by response shape (linear / weakly non-linear / non-linear)
8. Runs GO enrichment per module × ontology (BP, MF, CC) using topGO
   (classic + elim algorithms)

**Key result:** The two largest modules — turquoise (4,903 genes, enriched
for photosynthesis and light harvesting) and blue (2,006 genes, enriched for
stress response and protein quality control) — are strongly anticorrelated
(ME~turquoise~ vs ME~blue~: r = −0.71), representing opposing transcriptional
programmes along the light gradient.

---

## Repository contents

```
wgcna_fagus_lightgradient.R          Main analysis pipeline
README.md                            This file
```

---

## Input files

The following files must be present in the working directory before running
the script. All expression data and metadata are available from ENA
(PRJEB64934).

| File | Description |
|------|-------------|
| `expr_full_vst.txt` | VST-normalised expression matrix (genes × samples; 43,649 genes, 185 samples) |
| `sample_info_clean.txt` | Sample metadata (genotype, leaf position, lux, rel_expo, PI~abs~, Fv/Fm, drought group) |
| `gamm_leafpos_results.csv` | GAMM output: gene_id, response_shape, edf, delta_aic for all 15,126 responsive genes |
| `leafpos_responsive_genes.txt` | List of 15,126 leaf-position-responsive gene IDs (one per line) |
| `wgcna_gene_module_assignments.csv` | Gene → module table (output of Step 5; required as input for GO enrichment in Step 9) |
| `Bhaga_vs_Interpro5.69-101_GO_clean_TopGOinput` | Gene-to-GO term mapping for the *F. sylvatica* Bhaga reference genome (InterPro 5.69–101); tab-delimited, two columns: gene_id and comma-separated GO terms |

> **Note on the GO annotation file:** This file was generated from the
> *F. sylvatica* Bhaga reference genome annotation using InterPro release
> 5.69–101. It is not included in this repository due to size. Please contact
> the corresponding author or obtain it from the reference genome project if
> you wish to reproduce the GO enrichment step.

---

## Output files

| File | Description |
|------|-------------|
| `wgcna_sample_dendrogram.pdf` | Sample clustering dendrogram for outlier detection |
| `wgcna_softthreshold.pdf` | Scale-free topology fit and mean connectivity vs. soft-threshold power |
| `wgcna_module_dendrogram.pdf` | Gene dendrogram with module colour assignments |
| `wgcna_ME_correlation.pdf` | Module eigengene correlation dendrogram and heatmap |
| `wgcna_module_trait_heatmap.pdf` | Module–trait correlation heatmap (r, FDR-adjusted p) |
| `wgcna_network.rds` | Full `blockwiseModules` network object (R binary) |
| `wgcna_gene_module_assignments.csv` | Gene ID → module colour table |
| `wgcna_module_trait_r.csv` | Module–trait Pearson correlation matrix |
| `wgcna_module_trait_p_adj.csv` | BH-adjusted p-values for module–trait correlations |
| `wgcna_module_shape_summary.csv` | Response shape (linear/non-linear) distribution per module |
| `wgcna_topGO_all_results.csv` | All GO terms tested across modules and ontologies |
| `wgcna_topGO_significant.csv` | Significant GO terms (elim p < 0.05) |
| `wgcna_hub_genes.csv` | Top hub genes per module (highest kME) |
| `mt_long.csv` | Significant module–trait associations in long format |

Large intermediate files generated during the run (`wgcna_TOM-block.1.RData`)
are not committed to version control (see `.gitignore`).

---

## Dependencies

All analyses were run in R. Install the required packages before running:

```r
# CRAN
install.packages(c("dplyr", "ggplot2", "pheatmap"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("WGCNA", "topGO"))
```

| Package | Version used | Purpose |
|---------|-------------|---------|
| WGCNA | 1.72-5 | Network construction and module detection |
| topGO | 2.54.0 | GO enrichment analysis |
| dplyr | 1.1.4 | Data manipulation |
| ggplot2 | 3.5.1 | Visualisation |
| pheatmap | 1.0.12 | Heatmaps |

> The full session info (R version, all package versions) will be added
> on manuscript submission.

---

## Usage

1. Clone this repository
2. Download raw reads from ENA (PRJEB64934) and generate the VST-normalised
   expression matrix with your preferred RNA-seq pipeline (e.g. STAR + DESeq2
   `vst()`)
3. Place all input files listed above in the working directory
4. Run the script end-to-end or section by section:

```r
source("wgcna_fagus_lightgradient.R")
```

Steps 1–8 are sequential. Step 9 (GO enrichment) requires
`wgcna_gene_module_assignments.csv` from Step 5, so run Step 5 first or
provide a pre-computed assignment table.

**Runtime:** Steps 1–8 take approximately 20–40 minutes on a standard
workstation (4 cores, 16 GB RAM) for a 185-sample × 15,126-gene matrix.
Memory usage peaks at approximately 8 GB during TOM computation. If memory
is limited, reduce `maxBlockSize` to 10,000 in Step 5.

---

## Contact

Linda Eberhardt  
Senckenberg Biodiversity and Climate Research Centre (SBiK-F), Frankfurt  
linda.eberhardt@senckenberg.de 
GitHub: [github.com/leberhardt91](https://github.com/linda1903)

---

## Licence

Code is released under the MIT Licence. See `LICENSE` for details.


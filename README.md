# Single-Cell PBMC Analysis Pipeline

Welcome to the single-cell RNA-seq analysis repository! This project contains the R-based data processing, differential expression, and visualization pipelines used to compare patient samples (e.g., CD59-deficient and other neuroimmunological conditions) against healthy controls.

## 📁 Repository Structure

Here is a quick guide to what each file does. If you are a new developer or analyst joining the project, this is how the architecture is organized:

*   **`data_preparation.R`**
    The starting point. This script handles the initial ingestion, filtering, and cleaning of the raw single-cell datasets before they are fed into the main analytical pipelines.
*   **`functions.R`**
    The workhorse of the repository. To keep the analysis scripts readable, all custom logic—helper functions, complex plot generators (like custom forest plots or heatmaps), and data wrangling wrappers—are stored here. You will see `source("functions.R")` at the top of the other scripts.
*   **`multi_cell_analysis.R`**
    The master pipeline. This large script runs integration, meta-analyses, and pathway scoring across *all* immune cell types simultaneously. It generates the big-picture outputs, such as multi-cell UpSet plots, Over-Representation Analysis (ORA) comparisons, and broad pathway AUC heatmaps.
*   **`same_celltype-comparisons/` (Subfolder)**
    Contains cell-type-specific R Markdown (`.Rmd`) notebooks (e.g., `B_cells_analysis.Rmd`, `NK_anlysis.Rmd`, `CD4+_T_cells_analysis.Rmd`). While `multi_cell_analysis.R` looks at the whole system, these scripts zoom in. They generate granular, publication-ready figures (like UMAPs, specific GSEA Lollipops, and targeted forest plots) isolated to a single immune cell population.

## 🧠 Jargon Buster

If you are coming from a traditional data science background rather than a strict biology background, here are a few niche bioinformatics terms used heavily in this codebase:

*   **Seurat Object (`seu`)**: The core R data structure used heavily in this code. It stores single-cell data—including raw counts, metadata, and clustering data—all in one variable.
*   **SCT (SCTransform)**: A normalization method that uses a regularized negative binomial Generalized Linear Model (GLM) to remove technical noise (like sequencing depth variations) from the data without squashing actual biological signals.
*   **DE / DEGs (Differentially Expressed Genes)**: Finding which features (genes) are statistically turned "up" or "down" when comparing the patient group to the healthy group.
*   **GSEA & ORA**: Gene Set Enrichment Analysis & Over-Representation Analysis. Two statistical methods used to translate a massive list of differentially expressed genes (DEGs) into understandable high-level groups (e.g., instead of listing 50 individual genes, these tools tell us the "Metabolism" or "Inflammatory" pathways are active).
*   **AUC Score**: Area Under the Curve. In this specific context, it is used as a numerical score to quantify the overall activity level of a specific pathway in a given sample.
*   **`.qs` Files**: You'll see `qread()` and `qsave()` everywhere. The `qs` (Quick Serialization) package is just a highly optimized, multi-threaded alternative to R's native `saveRDS()`. It is used to quickly read/write massive data objects to the disk.

## 🚀 Getting Started

1.  **Hardware Requirements**: Single-cell data is massive. Many scripts here explicitly set `future.globals.maxSize` to handle objects upwards of 100-320 GB. Ensure you are running this on an appropriately provisioned cluster or server.
2.  **Environment Setup**: Check the top of `multi_cell_analysis.R` or the `.Rmd` files for the full list of required R libraries (e.g., `Seurat`, `ggplot2`, `ComplexHeatmap`, `qs`). 
3.  **Paths**: The scripts currently rely on absolute paths pointing to a `/bigdata/` directory. You may need to update these paths or set up environment variables depending on your local compute cluster.

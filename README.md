# Noncanonical Peptide Analysis Pipeline

This repository contains a comprehensive analysis pipeline for identifying and characterizing noncanonical peptides using proteomics and RNA-seq data.

## Project Structure

```
├── 0. Data processing.r            # Initial data processing and filtering
├── MHC affinity analysis.R        # MHC binding affinity prediction
├── NA, Survival Analysis.R        # Survival analysis for peptides
├── NMF clustering analysis.R      # NMF-based clustering analysis
├── Peptide_analysis.R            # Core peptide analysis functions
├── Translation Region check.R     # Translation region classification
├── Figure 2. RNA seq.R           # RNA-seq analysis pipeline
└── Figure 2. RNA seq_visualization.R  # RNA-seq visualization
```

## Dependencies

### R Packages
- Data Processing: tidyverse, dplyr, readxl, writexl
- Statistical Analysis: survival, survminer, ComplexHeatmap
- Visualization: ggplot2, pheatmap, viridis, gridExtra
- Bioinformatics: Biostrings, GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38
- Machine Learning: NMF, impute
- Others: Cairo, circlize, grid

## Installation

```R
# Install required packages
required_packages <- c(
  "tidyverse", "dplyr", "readxl", "writexl",
  "survival", "survminer", "ComplexHeatmap",
  "ggplot2", "pheatmap", "viridis", "gridExtra"
)

for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# BiocManager installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Biostrings",
  "GenomicFeatures",
  "BSgenome.Hsapiens.UCSC.hg38"
))
```



## Key Features

- Comprehensive peptide identification and characterization
- MHC binding affinity prediction
- Translation region classification
- Survival analysis
- RNA-seq integration
- Advanced visualization capabilities


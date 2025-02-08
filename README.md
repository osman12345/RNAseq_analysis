# RNAseq_analysis
# edgeR & limma Differential Expression Analysis Shiny App

This Shiny app performs RNA-seq differential expression analysis using the edgeR and limma packages. The app allows users to upload count data and sample information, run the analysis (including data filtering, normalization, voom transformation, and model fitting), and interactively view diagnostic plots and differential expression results for multiple contrasts.

## Features

- **Data Upload:** Upload CSV files for count data and sample metadata.
- **Data Processing:** Uses edgeR for data normalization and filtering and limma's voom transformation for modeling.
- **Contrast Selection:** Includes pairwise comparisons (Treated vs Control, Placebo vs Control, Treated vs Placebo) and an overall F-test (omnibus test) for the three groups.
- **Interactive Visualizations:** Displays library sizes, MD plots, boxplots, MDS plots, and heatmaps.
- **Results Download:** Download differential expression results for the selected contrast.

## Requirements

- **R (version 3.6 or later recommended)**
- **R Packages:**
  - [shiny](https://cran.r-project.org/package=shiny)
  - [edgeR](https://bioconductor.org/packages/edgeR/)
  - [limma](https://bioconductor.org/packages/limma/)
  - [org.Ss.eg.db](https://bioconductor.org/packages/org.Ss.eg.db/) *(for annotation)*
  - [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)
  - [gplots](https://cran.r-project.org/package=gplots)
  - [tidyverse](https://www.tidyverse.org/)
  - [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/)

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/your-repo-name.git
   cd your-repo-name


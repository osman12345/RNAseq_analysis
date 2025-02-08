# RNAseq Differential Expression Analysis Shiny App

This Shiny app performs RNAseq differential expression analysis using **edgeR** and **limma**. It allows users to upload count data and sample metadata, choose an annotation package to map gene symbols to ENTREZ IDs, and then run the analysis. The app generates diagnostic plots and tables for visualization and lets users download the differential expression results.

## Features

- **Data Upload:**  
  Upload CSV files for:
  - **Count Data:** Genes as rows and samples as columns (the first column contains gene symbols).
  - **Sample Metadata:** Must include a `Condition` column (e.g., "Control", "Treated").

- **Annotation Package Selection:**  
  Choose from multiple annotation packages:
  - Human (`org.Hs.eg.db`)
  - Mouse (`org.Mm.eg.db`)
  - Rat (`org.Rn.eg.db`)
  - Pig (`org.Ss.eg.db`)

- **Data Processing:**  
  - Filtering of low-expressed genes using **edgeR**.
  - Normalization of count data.
  - Voom transformation with quality weights.
  - Linear modeling and contrast definition using **limma**.

- **Differential Expression Analysis:**  
  The app defines contrasts (e.g., Treated vs Control) and outputs the differentially expressed genes (DEG) for further investigation.

- **Diagnostic Plots:**  
  - Library sizes
  - MD plot
  - Boxplots (unnormalized and normalized logCPM)
  - MDS plot
  - Heatmap of top variable genes

- **Results Download:**  
  Download the DEG table as a CSV file.

## Requirements

- **R** (version 3.6 or later recommended)

- **R Packages:**
  - [shiny](https://cran.r-project.org/package=shiny)
  - [edgeR](https://bioconductor.org/packages/edgeR/)
  - [limma](https://bioconductor.org/packages/limma/)
  - [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/)
  - [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/)
  - [org.Mm.eg.db](https://bioconductor.org/packages/org.Mm.eg.db/)
  - [org.Rn.eg.db](https://bioconductor.org/packages/org.Rn.eg.db/)
  - [org.Ss.eg.db](https://bioconductor.org/packages/org.Ss.eg.db/)
  - [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)
  - [gplots](https://cran.r-project.org/package=gplots)
  - [tidyverse](https://www.tidyverse.org/)

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/osman12345/RNAseq_analysis.git
   cd RNAseq_analysis

## Install Required Packages:
Open R or RStudio and run:
install.packages(c("shiny", "RColorBrewer", "gplots", "tidyverse"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "AnnotationDbi", 
                         "org.Hs.eg.db", "org.Mm.eg.db", 
                         "org.Rn.eg.db", "org.Ss.eg.db"))
## Data Requirements
** Count Data File:
A CSV file containing the gene expression counts.
Format:

- Rows: Genes (the first column contains gene symbols).
- Columns: Samples (each column header represents a sample).
** Sample Information File:
A CSV file containing sample metadata.
Format:
- Must include a column named Condition.
- Example values: "Control", "Treated".
## License
This project is licensed under the [MIT License](https://opensource.org/license/mit). 



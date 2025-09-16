# Mouse Heart snRNA-seq Analysis

This repository contains R scripts used to analyze mouse heart single-nucleus RNA-seq (snRNA-seq) data, focusing on cardiomyoblast and cardiomyocyte populations.  
Scripts are organized to allow reproduction of the analysis presented in the manuscript.

## Data
We analyzed four publicly available GEO datasets:
- GSM5355657_nonCM_Sham_1 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5355657))
- GSM5355658_nonCM_Sham_2 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5355658))
- GSM5943175_CM_Sham_1 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5943175))
- GSM5943176_CM_Sham_2 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5943176))

## Scripts
- `01_QC.R` – Quality control of raw data  
- `02_Merge.R` – Merging datasets and doublet removal  
- `03_Integration.R` – Normalization and integration (SCTransform)  
- `04_Clustering.R` – Dimensionality reduction, clustering, UMAP  
- `05_annotation.R` – Cell type annotation  
- `06_CM_subset.R` – Subsetting cardiac populations and re-clustering  
- `07_StemID2.R` – StemID2 stemness analysis (RaceID)  
- `08_GSEA_FigS2.R` – UpSet plot and GSEA/ssGSEA (Fig. S2)  
- `Fig1.R` – Integrated pipeline to generate Figure 1  

## Environment
- R version 4.4.3  
- Key packages: Seurat 5.3.0, sctransform 0.4.2, tidyverse 2.0.0, RaceID 0.3.9, GSVA 2.0.7, AUCell 1.28.0, escape 1.0.2  
- Full session information is provided in [sessionInfo_2025-09-16.txt](sessionInfo_2025-09-16.txt)

## Usage
1. Obtain raw data from GEO and preprocess as described.  
2. Run scripts sequentially (`01_QC.R` → `02_Merge.R` → ... → `Fig1.R`).  
3. Adjust the variable `base_dir` in each script to match your local directory structure (default: `./Heart`).  
4. Figures and intermediate RDS files will be saved under each script’s output folder.  

## License
MIT License (see `LICENSE`).

## Citation
This repository will be archived on Zenodo with a DOI upon manuscript submission.  
Please cite the corresponding DOI when using these scripts.

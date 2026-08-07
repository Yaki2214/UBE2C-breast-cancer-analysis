# Example Data and Input Data Sources

The full datasets used in the associated study are not redistributed
in this GitHub repository because of their file sizes and because
the original data are available from public repositories.

This directory provides information about the input files recognized
by `analysis.py` and, where applicable, small example files illustrating
the required data format.

## 1. Gene_expression.csv

### Source
cBioPortal for Cancer Genomics

Study:
Breast Invasive Carcinoma (TCGA, PanCancer Atlas)

Expression dataset:
RNA_Seq_v2_mRNA_median_all_sample_Zscores

### Description
This file contains gene-expression Z-scores used for gene-level
correlation, stratification, and screening analyses.

The original expression dataset should be downloaded from cBioPortal
and formatted as `Gene_expression.csv` before running the program.

The analysis program expects genes and samples to be arranged in the
format demonstrated in `Gene_expression_example.csv`.


## 2. Clinical_data.csv

### Source
cBioPortal for Cancer Genomics

Study:
Breast Invasive Carcinoma (TCGA, PanCancer Atlas)

Source files:
- clinical_patient
- clinical_sample

### Description
Patient-level and sample-level clinical annotations downloaded from
cBioPortal were combined into a single table for analysis.

The resulting table was saved as:

`Clinical_data.csv`

Sample identifiers must correspond to those used in the expression
matrix.


## 3. Signaling_pathway.csv

### Source
Derived dataset.

### Description
Pathway activity scores were calculated from the gene-expression
matrix using single-sample Gene Set Enrichment Analysis (ssGSEA).

Gene sets were defined using the pathway-gene mappings contained in:

`Full_Gene_List.csv`

The resulting sample-by-pathway score matrix was saved as:

`Signaling_pathway.csv`

This derived dataset was used for pathway-level correlation,
stratification, and screening analyses.


## 4. Full_Gene_List.csv

This file defines the genes belonging to each signaling pathway used
for ssGSEA-based pathway scoring.

The pathway definitions used for the manuscript-associated analyses
should be retained with the corresponding code release whenever
redistribution is permitted.


## 5. interactions.csv

### Source
Drug-Gene Interaction Database (DGIdb)

### Description
Drug-gene interaction annotations were downloaded from DGIdb and
formatted as `interactions.csv`.

These annotations are used by `analysis.py` to annotate genes with
known drug interactions.

The annotation file is not required for statistical calculations when
drug-target annotation is not used.


## 6. Drug annotation data

Drug and drug-gene interaction information used by the analysis tool
was obtained from DGIdb.

Because DGIdb is periodically updated, the database version and/or
download date used for the manuscript-associated analysis should be
reported whenever available.


## 7. Knockdown_Dependency_Score.csv

### Source
DepMap Portal

### Description
CRISPR gene-dependency data were obtained from the DepMap Portal and
formatted for use by the dependency-analysis module of `analysis.py`.

The DepMap release used for the manuscript-associated analysis should
be specified because DepMap datasets are periodically updated.

The complete DepMap dataset is not redistributed in this repository.
Users should obtain the appropriate release directly from the
DepMap Portal.


## Reproducing the input datasets

The recommended workflow is:

1. Download TCGA-BRCA expression and clinical data from cBioPortal.
2. Prepare `Gene_expression.csv`.
3. Merge the cBioPortal patient- and sample-level clinical annotations
   to generate `Clinical_data.csv`.
4. Calculate ssGSEA pathway scores using `Full_Gene_List.csv` to
   generate `Signaling_pathway.csv`.
5. Obtain drug-gene interaction information from DGIdb when
   drug-target annotation is required.
6. Obtain the appropriate CRISPR dependency dataset from DepMap
   when dependency analysis is required.
7. Place the required files in the same working directory as
   `analysis.py`.

Not all input files are required for every analysis module.

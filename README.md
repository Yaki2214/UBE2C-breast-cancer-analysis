# Integrated Gene Analysis Tool

This repository contains the Python analysis code used for selected
bioinformatic and statistical analyses in the study of
UBE2C-associated genomic instability and therapeutic vulnerability
in breast cancer.

The analysis workflow is implemented as a single integrated,
menu-driven Python program (`analysis.py`).

## Main functions

The program provides the following analysis modules:

### A. Correlation and grouping analysis
1. Two-target correlation scatter plots
2. Batch correlation matrices
3. Patient co-expression heatmaps
4. Patient waterfall plots
5. High/low group comparisons

### B. Clinical and survival analysis
6. Clinical feature distribution analysis
7. Kaplan-Meier survival analysis
8. Dual-target Kaplan-Meier analysis
9. Cox proportional hazards analysis

### C. Large-scale screening
10. Correlation screening
11. Differential expression screening

### D. Feature engineering
12. Virtual gene-signature generation

### E. External-data visualization
13. Bar/lollipop plots
14. Dual comparison plots
15. Volcano plots
16. Waterfall plots

### F. CRISPR dependency analysis
17. Dependency-score waterfall plots


## Installation

Python 3.14.0 is required.

Install the required Python packages using:

    pip install -r requirements.txt


## Running the program

Run:

    python analysis.py

The program will automatically detect available input files and
display the corresponding analysis options.


## Input files

Depending on the analysis being performed, the program recognizes
the following input files:

- `Gene_expression.csv`
- `Signaling_pathway.csv`
- `Clinical_data.csv`
- `protein_zscores_TCGA.csv`
- `Signaling_pathway_Protein.csv`
- `Knockdown_Dependency_Score.csv`
- `interactions.csv`
- `Full_Gene_List.csv`

Not all files are required for every analysis module.

Input files should be placed in the same working directory in which
the program is executed unless otherwise specified.


## Output

Results are automatically written to a timestamped directory:

    Desktop/GeneAnalysis_YYYYMMDD_HHMMSS/

Depending on the selected analysis, output may include:

- PNG/PDF/TIFF/SVG figures
- Excel workbooks
- CSV files
- correlation statistics
- survival-analysis results
- Cox regression results
- differential-screening results


## Statistical methods

Depending on the selected module, statistical methods implemented
in the program include:

- Pearson correlation
- Welch's t-test
- one-way ANOVA
- Benjamini-Hochberg false-discovery-rate correction
- log-rank tests
- Cox proportional hazards regression
- median-based high/low stratification

Please refer to the associated manuscript for the specific
statistical methods used for each analysis.


## Data availability

The program is designed to analyze user-supplied or publicly
available datasets. Public datasets used in the associated
manuscript should be obtained from their original repositories.

No patient-identifiable information is included in this repository.


## Reproducibility

The manuscript-associated version of this software is archived as
release [VERSION].

The source code contained in that release corresponds to the
version used for the analyses reported in the manuscript.

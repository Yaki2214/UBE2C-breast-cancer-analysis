# Data Processing

This directory contains scripts used to generate derived input data
for the integrated analysis tool.

### Software environment

The manuscript-associated ssGSEA analysis was performed using:

- R version 4.5.2
- GSVA version 2.4.0

## ssGSEA pathway scoring

Pathway activity scores were generated from the gene-expression
matrix using single-sample Gene Set Enrichment Analysis (ssGSEA)
implemented in the R package GSVA.

The workflow is implemented in:

`calculate_ssgsea.R`

### Input files

The script requires:

- `Gene_expression.csv`
- `Full_Gene_List.csv`

`Gene_expression.csv` contains gene-expression values, with genes
represented as rows and samples represented as columns.

`Full_Gene_List.csv` contains the pathway-specific gene sets used
for ssGSEA analysis.

### Analysis workflow

1. The expression matrix is imported into R.
2. Pathway-specific gene sets are converted into an R list.
3. ssGSEA scores are calculated using `ssgseaParam()` and `gsva()`
   from the GSVA package.
4. For each pathway, ssGSEA scores are min-max normalized across
   samples to a range from 0 to 1.
5. The normalized pathway score matrix is exported as
   `Signaling_pathway.csv`.

The resulting `Signaling_pathway.csv` can then be used directly by
`analysis.py`.

### Running the script

From R, set the working directory to the directory containing the
required input files and run:

    source("calculate_ssgsea.R")

Alternatively, from a command line with Rscript available:

    Rscript calculate_ssgsea.R

### Output

The script generates:

`Signaling_pathway.csv`

This file contains the normalized ssGSEA pathway activity scores
used for subsequent pathway-level analyses.

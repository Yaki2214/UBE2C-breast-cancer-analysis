# ============================================================
# ssGSEA pathway score calculation
# ============================================================

library(GSVA)

# ------------------------------------------------------------
# 1. Load gene expression matrix
# ------------------------------------------------------------
# Rows: genes
# Columns: samples

expression_matrix <- as.matrix(
    read.csv(
        "Gene_expression.csv",
        row.names = 1,
        check.names = FALSE
    )
)

storage.mode(expression_matrix) <- "numeric"


# ------------------------------------------------------------
# 2. Load pathway gene sets
# ------------------------------------------------------------
# Each column represents one pathway.
# Non-missing entries in each column represent genes belonging
# to that pathway.

gene_sets <- read.csv(
    "Full_Gene_List.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

gene_set_list <- lapply(
    gene_sets,
    function(x) as.character(x[!is.na(x) & x != ""])
)


# ------------------------------------------------------------
# 3. Create ssGSEA parameter object
# ------------------------------------------------------------

ssgsea_params <- ssgseaParam(
    expression_matrix,
    gene_set_list
)


# ------------------------------------------------------------
# 4. Calculate ssGSEA scores
# ------------------------------------------------------------

ssgsea_scores <- gsva(ssgsea_params)


# ------------------------------------------------------------
# 5. Min-max normalization
# ------------------------------------------------------------
# Scores for each pathway are independently scaled across
# samples to a range from 0 to 1.

ssgsea_scores_normalized <- t(
    apply(
        ssgsea_scores,
        1,
        function(x) {
            score_range <- max(x) - min(x)

            if (score_range == 0) {
                return(rep(0, length(x)))
            }

            (x - min(x)) / score_range
        }
    )
)


# ------------------------------------------------------------
# 6. Export results
# ------------------------------------------------------------

write.csv(
    ssgsea_scores_normalized,
    "Signaling_pathway.csv"
)

cat("ssGSEA pathway scoring completed.\n")
cat("Output: Signaling_pathway.csv\n")

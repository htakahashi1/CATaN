#' Prepare TF and transcriptome matrices for CCA
#'
#' Takes a filtered count matrix and a TF-GRN matrix, finds the
#' intersection of genes, removes duplicates, and performs TMM normalisation
#' followed by log2(CPM + 1) transformation on the transcriptome data.
#'
#' @param counts A numeric matrix of filtered raw counts (genes x samples),
#'   typically the output of \code{\link{filter_low_expression}}.
#' @param tf_matrix A numeric matrix of TF-gene connectivity scores
#'   (genes x TFs). Row names must be gene identifiers matching those in
#'   \code{counts}.
#' @param norm_method Character. Normalisation method passed to
#'   \code{edgeR::calcNormFactors()} (default \code{"TMM"}).
#' @param verbose Logical. Print summary messages (default TRUE).
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{tf}{Numeric matrix of TF binding data for intersecting genes
#'     (gene x TF).}
#'   \item{expr}{Numeric matrix of TMM-normalised log2(CPM + 1) values for
#'     intersecting genes (gene x sample).}
#' }
#'
#' @details
#' Gene identifiers are matched by row names. Genes duplicated in either
#' matrix are removed prior to intersection. The transcriptome matrix is
#' normalised using edgeR's TMM method and transformed to log2(CPM + 1)
#' scale.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rnbinom(500, mu = 50, size = 2), nrow = 50, ncol = 10)
#' rownames(counts) <- paste0("gene", seq_len(50))
#' tf_mat <- matrix(sample(0:1, 50 * 5, replace = TRUE), nrow = 50, ncol = 5)
#' rownames(tf_mat) <- paste0("gene", seq_len(50))
#' colnames(tf_mat) <- paste0("TF", seq_len(5))
#' result <- prepare_matrices(counts, tf_mat)
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @export
prepare_matrices <- function(counts,
                             tf_matrix,
                             norm_method = "TMM",
                             verbose     = TRUE) {

    ## Input validation
    if (!is.matrix(counts) || !is.numeric(counts)) {
        stop("counts must be a numeric matrix")
    }
    if (!is.matrix(tf_matrix) || !is.numeric(tf_matrix)) {
        stop("tf_matrix must be a numeric matrix")
    }
    if (is.null(rownames(counts)) || is.null(rownames(tf_matrix))) {
        stop("Both counts and tf_matrix must have row names (gene identifiers)")
    }

    ## Remove duplicated genes from TF matrix
    dup_tf <- duplicated(rownames(tf_matrix)) |
              duplicated(rownames(tf_matrix), fromLast = TRUE)
    if (any(dup_tf)) {
        if (verbose) {
            message(sprintf("Removing %d duplicated genes from TF matrix",
                            sum(dup_tf)))
        }
        tf_matrix <- tf_matrix[!dup_tf, , drop = FALSE]
    }

    ## Remove duplicated genes from count matrix
    dup_count <- duplicated(rownames(counts)) |
                 duplicated(rownames(counts), fromLast = TRUE)
    if (any(dup_count)) {
        if (verbose) {
            message(sprintf("Removing %d duplicated genes from count matrix",
                            sum(dup_count)))
        }
        counts <- counts[!dup_count, , drop = FALSE]
    }

    ## Intersect genes
    common_genes <- intersect(rownames(counts), rownames(tf_matrix))
    if (length(common_genes) == 0) {
        stop("No common genes found between counts and tf_matrix. ",
             "Check that row names use the same gene identifier format.")
    }

    tf_out    <- tf_matrix[common_genes, , drop = FALSE]
    counts_out <- counts[common_genes, , drop = FALSE]

    ## TMM normalisation -> log2(CPM + 1)
    dge    <- edgeR::DGEList(counts = counts_out)
    dge    <- edgeR::calcNormFactors(dge, method = norm_method)
    logcpm <- log2(edgeR::cpm(dge) + 1)

    if (verbose) {
        message(sprintf(
            "Prepared matrices: %d common genes, %d TFs, %d samples",
            length(common_genes), ncol(tf_out), ncol(logcpm)))
    }

    list(tf = tf_out, expr = logcpm)
}

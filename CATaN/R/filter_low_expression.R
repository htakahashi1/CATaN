#' Filter low-expression genes from a count matrix
#'
#' Removes genes that are not expressed above specified thresholds in a
#' sufficient proportion of samples. Applies two filters: one on raw counts
#' and one on CPM values. Only genes passing both filters are retained.
#'
#' @param counts A numeric matrix of raw counts (genes x samples) or a
#'   \code{SummarizedExperiment} object.
#' @param min_count Numeric. Minimum raw count threshold (default 10).
#' @param min_cpm Numeric. Minimum CPM threshold (default 1).
#' @param min_prop Numeric. Minimum proportion of samples that must exceed both
#'   thresholds (default 0.15).
#' @param verbose Logical. Print filtering summary (default TRUE).
#'
#' @return A numeric matrix of filtered counts (genes x samples).
#'
#' @details
#' This function replicates the filtering strategy used in the CATaN pipeline:
#' \enumerate{
#'   \item Compute CPM using \code{edgeR::cpm()}.
#'   \item Identify genes where at least \code{min_prop} of samples have
#'         raw count >= \code{min_count}.
#'   \item Identify genes where at least \code{min_prop} of samples have
#'         CPM >= \code{min_cpm}.
#'   \item Retain genes passing both criteria.
#' }
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rnbinom(1000, mu = 10, size = 1), nrow = 100, ncol = 10)
#' rownames(counts) <- paste0("gene", seq_len(100))
#' filtered <- filter_low_expression(counts)
#'
#' @importFrom edgeR cpm
#' @export
filter_low_expression <- function(counts,
                                  min_count = 10,
                                  min_cpm   = 1,
                                  min_prop  = 0.15,
                                  verbose   = TRUE) {

    ## Handle SummarizedExperiment input
    if (is(counts, "SummarizedExperiment")) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("SummarizedExperiment package required for SE input")
        }
        counts <- SummarizedExperiment::assay(counts, "counts")
    }

    ## Input validation
    if (!is.matrix(counts) || !is.numeric(counts)) {
        stop("counts must be a numeric matrix")
    }
    if (is.null(rownames(counts))) {
        stop("counts must have row names (gene identifiers)")
    }

    n_genes_orig <- nrow(counts)
    n_samples    <- ncol(counts)
    min_samples  <- ceiling(n_samples * min_prop)

    ## Filter 1: raw counts
    pass_count <- rowSums(counts >= min_count) >= min_samples

    ## Filter 2: CPM
    cpm_mat    <- edgeR::cpm(counts)
    pass_cpm   <- rowSums(cpm_mat >= min_cpm) >= min_samples

    ## Intersection of both filters
    keep <- pass_count & pass_cpm
    filtered <- counts[keep, , drop = FALSE]

    if (verbose) {
        message(sprintf("Low expression filtering: %d -> %d genes (removed %d)",
                        n_genes_orig, nrow(filtered),
                        n_genes_orig - nrow(filtered)))
        message(sprintf(
            "  Criteria: count >= %g and CPM >= %g in >= %.0f%% samples",
            min_count, min_cpm, min_prop * 100))
    }

    filtered
}

#' Extract top and bottom percentile SNPs by CC score
#'
#' For each CC in a CATaNAnnotation object, the top and bottom 10% of
#' SNPs are extracted.
#'
#' @param annotation A \code{\linkS4class{CATaNAnnotation}} object with
#'   the \code{snp_scores} slot populated (output of
#'   \code{\link{align_cc_to_snps}}).
#' @param percentile Numeric. Fraction of SNPs to extract from each tail
#'   (default 0.1 = 10%%).
#'
#' @return A \code{\linkS4class{CATaNAnnotation}} object with the
#'   \code{top_snps} and \code{bottom_snps} slots populated.
#'
#' @examples
#' # Create a mock CATaNAnnotation with scores
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges::IRanges(seq(1, 100), width = 1))
#' mcols(gr)$score <- rnorm(100)
#' annot <- new("CATaNAnnotation",
#'     snp_scores = list(CC1_rotation = gr),
#'     percentile = 0.1, population = "test")
#' result <- extract_top_bottom_snps(annot)
#' length(topSnps(result)$CC1_rotation)
#'
#' @details
#' For each CC, SNPs are sorted by their aggregated CC score. The top
#' \code{percentile} fraction (by count) with the highest scores and the
#' bottom \code{percentile} fraction with the lowest scores are extracted.
#' This matches the original shell script logic:
#' \code{sort -g -k 4,4 | tail -n nrow/10} (top) and
#' \code{sort -g -k 4,4 | head -n nrow/10} (bottom).
#'
#' @importFrom GenomicRanges mcols
#' @importFrom methods is
#' @export
extract_top_bottom_snps <- function(annotation, percentile = 0.1) {

    ## ---- Input validation ----
    if (!is(annotation, "CATaNAnnotation")) {
        stop("annotation must be a CATaNAnnotation object")
    }
    if (length(annotation@snp_scores) == 0) {
        stop("snp_scores is empty. Run align_cc_to_snps() first.")
    }
    if (percentile <= 0 || percentile >= 1) {
        stop("percentile must be between 0 and 1 (exclusive)")
    }

    top_list    <- list()
    bottom_list <- list()

    for (cc_name in names(annotation@snp_scores)) {
        gr <- annotation@snp_scores[[cc_name]]
        n  <- length(gr)

        if (n == 0) {
            top_list[[cc_name]]    <- GenomicRanges::GRanges()
            bottom_list[[cc_name]] <- GenomicRanges::GRanges()
            next
        }

        ## Number of SNPs in each tail (integer division, matching original)
        n_tail <- n %/% as.integer(1 / percentile)

        ## Sort by score
        scores   <- GenomicRanges::mcols(gr)$score
        sort_idx <- order(scores)

        ## Bottom = lowest scores (head of sorted)
        bottom_idx <- sort_idx[seq_len(n_tail)]
        bottom_list[[cc_name]] <- sort(gr[bottom_idx])

        ## Top = highest scores (tail of sorted)
        top_idx <- sort_idx[seq(n - n_tail + 1, n)]
        top_list[[cc_name]] <- sort(gr[top_idx])
    }

    ## Update annotation object
    annotation@top_snps    <- top_list
    annotation@bottom_snps <- bottom_list
    annotation@percentile  <- percentile

    annotation
}

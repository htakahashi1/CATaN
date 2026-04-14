#' @title CATaNResult class
#'
#' @description S4 class to store results from CCA of TF-binding and
#'   transcriptome matrices. This is the output of \code{\link{run_cca}} and
#'   serves as input to \code{\link{align_cc_to_snps}}.
#'
#' @slot tf_matrix Filtered TF binding matrix (gene x TF) used in CCA.
#' @slot transcriptome Filtered and normalised logCPM matrix (gene x sample)
#'   used in CCA.
#' @slot hvg Character vector of highly variable gene names selected for CCA.
#' @slot svd_result List containing the raw SVD output (\code{u}, \code{d},
#'   \code{v}).
#' @slot n_cc Integer, the number of canonical components retained.
#' @slot tf_sample_loading DataFrame of TF weights (U matrix), TF x CC.
#' @slot transcriptome_sample_loading DataFrame of sample weights (V matrix), sample x CC.
#' @slot parameters DataFrame of CCA summary parameters including variance
#'   explained, SCF, and canonical correlations.
#' @slot sample_metadata DataFrame of sample annotations (optional).
#' @slot call The matched call used to produce this object.
#'
#' @importFrom methods setClass new validObject
#' @importFrom S4Vectors DataFrame
#' @export
setClass("CATaNResult",
    slots = list(
        tf_matrix       = "matrix",
        transcriptome   = "matrix",
        hvg             = "character",
        svd_result      = "list",
        n_cc            = "integer",
        tf_sample_loading  = "DataFrame",
        transcriptome_sample_loading  = "DataFrame",
        parameters      = "DataFrame",
        sample_metadata = "DataFrame",
        call            = "ANY"
    ),
    prototype = list(
        tf_matrix       = matrix(0, 0, 0),
        transcriptome   = matrix(0, 0, 0),
        hvg             = character(0),
        svd_result      = list(),
        n_cc            = 10L,
        tf_sample_loading      = S4Vectors::DataFrame(),
        transcriptome_sample_loading  = S4Vectors::DataFrame(),
        parameters      = S4Vectors::DataFrame(),
        sample_metadata = S4Vectors::DataFrame(),
        call            = NULL
    ),
    validity = function(object) {
        errors <- character()

        ## n_cc must be positive
        if (object@n_cc < 1L) {
            errors <- c(errors, "n_cc must be a positive integer")
        }

        ## SVD components must be present if svd_result is non-empty
        if (length(object@svd_result) > 0) {
            required <- c("u", "d", "v")
            missing <- setdiff(required, names(object@svd_result))
            if (length(missing) > 0) {
                errors <- c(errors,
                    paste("svd_result must contain:",
                          paste(required, collapse = ", ")))
            }
        }

        ## tf_sample_loading columns must match n_cc
        if (ncol(object@tf_sample_loading) > 0 &&
            ncol(object@tf_sample_loading) != object@n_cc) {
            errors <- c(errors,
                sprintf("tf_sample_loading has %d columns but n_cc is %d",
                        ncol(object@tf_sample_loading), object@n_cc))
        }

        ## transcriptome_sample_loading columns must match n_cc
        if (ncol(object@transcriptome_sample_loading) > 0 &&
            ncol(object@transcriptome_sample_loading) != object@n_cc) {
            errors <- c(errors,
                sprintf("transcriptome_sample_loading has %d columns but n_cc is %d",
                        ncol(object@transcriptome_sample_loading), object@n_cc))
        }

        if (length(errors) == 0) TRUE else errors
    }
)


#' @title CATaNAnnotation class
#'
#' @description S4 class to store SNP-level CC score annotations for sLDSC.
#'   This is the output of \code{\link{align_cc_to_snps}} and
#'   \code{\link{extract_top_bottom_snps}}.
#'
#' @slot snp_scores Named list of GRanges, one per CC. Each GRanges contains
#'   all SNPs with their aggregated CC score in \code{mcols(.)$score}.
#' @slot top_snps Named list of GRanges, one per CC. SNPs in the top
#'   percentile of CC scores.
#' @slot bottom_snps Named list of GRanges, one per CC. SNPs in the bottom
#'   percentile of CC scores.
#' @slot percentile Numeric, the fraction used for top/bottom extraction
#'   (default 0.1).
#' @slot population Character describing the reference population
#'   (e.g. "EUR", "EAS", "custom").
#' @slot catan_result The CATaNResult object used to generate this annotation.
#'
#' @importFrom methods setClass new
#' @importFrom GenomicRanges GRanges
#' @export
setClass("CATaNAnnotation",
    slots = list(
        snp_scores    = "list",
        top_snps      = "list",
        bottom_snps   = "list",
        percentile    = "numeric",
        population    = "character",
        catan_result  = "CATaNResult"
    ),
    prototype = list(
        snp_scores    = list(),
        top_snps      = list(),
        bottom_snps   = list(),
        percentile    = 0.1,
        population    = "custom"
    ),
    validity = function(object) {
        errors <- character()

        if (object@percentile <= 0 || object@percentile >= 1) {
            errors <- c(errors,
                "percentile must be between 0 and 1 (exclusive)")
        }

        ## All elements of snp_scores must be GRanges
        if (length(object@snp_scores) > 0) {
            is_gr <- vapply(object@snp_scores, is, logical(1), "GRanges")
            if (!all(is_gr)) {
                errors <- c(errors,
                    "All elements of snp_scores must be GRanges")
            }
        }

        if (length(errors) == 0) TRUE else errors
    }
)


## ---- Show methods ----

#' @rdname CATaNResult-class
#' @importFrom methods show
#' @return Prints a summary of the object to the console
#'   (invisibly returns NULL).
#' @examples
#' res <- new("CATaNResult")
#' show(res)
#' @export
setMethod("show", "CATaNResult", function(object) {
    cat("CATaNResult object\n")
    cat(sprintf("  Genes (HVG): %d\n", length(object@hvg)))
    cat(sprintf("  TFs: %d\n", ncol(object@tf_matrix)))
    cat(sprintf("  Samples: %d\n", ncol(object@transcriptome)))
    cat(sprintf("  Canonical components: %d\n", object@n_cc))
    if (nrow(object@parameters) > 0) {
        cat("  Parameters available: ",
            paste(rownames(object@parameters), collapse = ", "), "\n")
    }
})

#' @rdname CATaNAnnotation-class
#' @return Prints a summary of the object to the console
#'   (invisibly returns NULL).
#' @examples
#' annot <- new("CATaNAnnotation")
#' show(annot)
#' @export
setMethod("show", "CATaNAnnotation", function(object) {
    cat("CATaNAnnotation object\n")
    cat(sprintf("  CCs annotated: %d\n", length(object@snp_scores)))
    if (length(object@snp_scores) > 0) {
        n_snps <- length(object@snp_scores[[1]])
        cat(sprintf("  SNPs per CC: %d\n", n_snps))
    }
    cat(sprintf("  Percentile: %.0f%%\n", object@percentile * 100))
    cat(sprintf("  Population: %s\n", object@population))
    has_top <- length(object@top_snps) > 0
    cat(sprintf("  Top/bottom SNPs: %s\n",
                ifelse(has_top, "extracted", "not yet")))
})

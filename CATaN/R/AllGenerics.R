#' @title Accessor generics for CATaN classes
#' @name CATaN-accessors
#' @description Generic accessor functions for CATaNResult and
#'   CATaNAnnotation objects.
#'
#' @param object A CATaNResult or CATaNAnnotation object.
#' @return A DataFrame or list depending on the accessor.
#'
#' @examples
#' set.seed(42)
#' tf <- matrix(sample(0:1, 500, replace = TRUE), 50, 10)
#' rownames(tf) <- paste0("g", seq_len(50))
#' colnames(tf) <- paste0("TF", seq_len(10))
#' ex <- matrix(rnorm(250), 50, 5)
#' rownames(ex) <- rownames(tf)
#' colnames(ex) <- paste0("s", seq_len(5))
#' res <- run_cca(tf, ex, n_cc = 3L, n_hvg = 30L,
#'                verbose = FALSE)
#' tfsampleLoading(res)
#' trasampleLoading(res)
#' ccParameters(res)
NULL

#' @rdname CATaN-accessors
#' @export
setGeneric("tfsampleLoading", function(object) standardGeneric("tfsampleLoading"))

#' @rdname CATaN-accessors
#' @export
setGeneric("trasampleLoading", function(object) standardGeneric("trasampleLoading"))

#' @rdname CATaN-accessors
#' @export
setGeneric("ccParameters", function(object) standardGeneric("ccParameters"))

#' @rdname CATaN-accessors
#' @export
setGeneric("snpScores", function(object) standardGeneric("snpScores"))

#' @rdname CATaN-accessors
#' @export
setGeneric("topSnps", function(object) standardGeneric("topSnps"))

#' @rdname CATaN-accessors
#' @export
setGeneric("bottomSnps", function(object) standardGeneric("bottomSnps"))


## ---- Accessor methods for CATaNResult ----

#' @rdname CATaN-accessors
#' @export
setMethod("tfsampleLoading", "CATaNResult", function(object) object@tf_sample_loading)

#' @rdname CATaN-accessors
#' @export
setMethod("trasampleLoading", "CATaNResult",
    function(object) object@transcriptome_sample_loading)

#' @rdname CATaN-accessors
#' @export
setMethod("ccParameters", "CATaNResult", function(object) object@parameters)


## ---- Accessor methods for CATaNAnnotation ----

#' @rdname CATaN-accessors
#' @export
setMethod("snpScores", "CATaNAnnotation", function(object) object@snp_scores)

#' @rdname CATaN-accessors
#' @export
setMethod("topSnps", "CATaNAnnotation", function(object) object@top_snps)

#' @rdname CATaN-accessors
#' @export
setMethod("bottomSnps", "CATaNAnnotation", function(object) object@bottom_snps)

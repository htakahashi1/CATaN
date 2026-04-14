#' Run Canonical Correlation Analysis of TF-binding and transcriptome data
#'
#' Performs CCA between a TF-gene binding matrix (X) and a transcriptome
#' matrix (Y) via SVD of the cross-product matrix X'Y. Returns a
#' \code{\linkS4class{CATaNResult}} object containing canonical component
#' weights and summary parameters.
#'
#' @param tf_matrix Numeric matrix of TF-gene binding associations
#'   (gene x TF). Typically from \code{\link{prepare_matrices}}.
#' @param transcriptome Numeric matrix of normalised expression values
#'   (gene x sample). Typically log2(CPM + 1) from
#'   \code{\link{prepare_matrices}}.
#' @param n_cc Integer. Number of canonical components to retain
#'   (default 10L).
#' @param n_hvg Integer. Number of top highly variable genes to select from
#'   each matrix before taking their intersection (default 10000L).
#' @param sample_metadata A \code{data.frame} or \code{DataFrame} of sample
#'   annotations (optional). Row names or a column \code{"id"} should match
#'   column names of \code{transcriptome}.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A \code{\linkS4class{CATaNResult}} object.
#'
#' @details
#' The analysis proceeds as follows:
#' \enumerate{
#'   \item \strong{HVG selection}: For each matrix, gene-wise variance is
#'     computed and the top \code{n_hvg} genes are selected. The intersection
#'     of these two sets defines the gene set used for CCA.
#'   \item \strong{Zero-TF removal}: TF columns where all values are zero
#'     across the selected genes are removed.
#'   \item \strong{Bidirectional scaling}: Both matrices are z-score
#'     normalised first across genes (row-wise) then across features
#'     (column-wise). Genes producing NA (zero variance) are removed.
#'   \item \strong{SVD}: The cross-product matrix \eqn{X^T Y} is computed
#'     and decomposed via SVD. The left singular vectors (U) give TF weights
#'     and the right singular vectors (V) give sample weights for each CC.
#'   \item \strong{Parameter computation}: Variance explained, Squared
#'     Canonical Fraction (SCF), and canonical correlations are computed.
#' }
#'
#' @examples
#' set.seed(42)
#' ngenes <- 200; ntf <- 20; nsamp <- 30
#' tf_mat <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
#'                  nrow = ngenes, ncol = ntf)
#' rownames(tf_mat) <- paste0("gene", seq_len(ngenes))
#' colnames(tf_mat) <- paste0("TF", seq_len(ntf))
#' expr_mat <- matrix(rnorm(ngenes * nsamp), nrow = ngenes, ncol = nsamp)
#' rownames(expr_mat) <- paste0("gene", seq_len(ngenes))
#' colnames(expr_mat) <- paste0("sample", seq_len(nsamp))
#' res <- run_cca(tf_mat, expr_mat, n_cc = 5L, n_hvg = 100L)
#'
#' @importFrom stats var cor
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame
#' @export
run_cca <- function(tf_matrix,
                    transcriptome,
                    n_cc            = 10L,
                    n_hvg           = 10000L,
                    sample_metadata = NULL,
                    verbose         = TRUE) {

    mc <- match.call()
    n_cc  <- as.integer(n_cc)
    n_hvg <- as.integer(n_hvg)

    ## ---- Input validation ----
    if (!is.matrix(tf_matrix) || !is.numeric(tf_matrix)) {
        stop("tf_matrix must be a numeric matrix")
    }
    if (!is.matrix(transcriptome) || !is.numeric(transcriptome)) {
        stop("transcriptome must be a numeric matrix")
    }
    if (is.null(rownames(tf_matrix)) || is.null(rownames(transcriptome))) {
        stop("Both matrices must have row names (gene identifiers)")
    }

    common_genes <- intersect(rownames(tf_matrix), rownames(transcriptome))
    if (length(common_genes) == 0) {
        stop("No common genes between tf_matrix and transcriptome")
    }

    ## ---- Step 1: HVG selection ----
    if (verbose) message("Selecting highly variable genes...")

    tf_var   <- apply(tf_matrix[common_genes, , drop = FALSE], 1, var)
    expr_var <- apply(transcriptome[common_genes, , drop = FALSE], 1, var)

    n_hvg_actual <- min(n_hvg, length(common_genes))
    top_tf   <- names(sort(tf_var, decreasing = TRUE))[seq_len(n_hvg_actual)]
    top_expr <- names(sort(expr_var, decreasing = TRUE))[seq_len(n_hvg_actual)]
    hvg      <- intersect(top_tf, top_expr)

    if (length(hvg) < n_cc) {
        stop(sprintf("Only %d HVGs found, fewer than n_cc = %d. ",
                     length(hvg), n_cc),
             "Decrease n_cc or increase n_hvg.")
    }

    if (verbose) message(sprintf("  %d HVGs selected", length(hvg)))

    X <- tf_matrix[hvg, , drop = FALSE]
    Y <- transcriptome[hvg, , drop = FALSE]

    ## ---- Step 2: Remove all-zero TF columns ----
    tf_max <- apply(X, 2, max)
    allzero_tfs <- names(tf_max[tf_max == 0])
    if (length(allzero_tfs) > 0) {
        if (verbose) {
            message(sprintf("  Removing %d all-zero TF columns",
                            length(allzero_tfs)))
        }
        X <- X[, !colnames(X) %in% allzero_tfs, drop = FALSE]
    }

    if (ncol(X) < n_cc) {
        stop(sprintf("Only %d non-zero TFs, fewer than n_cc = %d",
                     ncol(X), n_cc))
    }

    ## ---- Step 3: Bidirectional scaling ----
    if (verbose) message("Performing bidirectional scaling...")

    ## Row-wise (gene-wise) z-score -> transpose -> column-wise z-score
    scaled_X <- .bidirectional_scale(X)
    scaled_Y <- .bidirectional_scale(Y)

    ## Remove NA columns (zero-variance features) from X
    na_cols_x <- apply(scaled_X, 2, function(col) any(is.na(col)))
    if (any(na_cols_x)) {
        if (verbose) {
            message(sprintf("  Removing %d zero-variance TF columns",
                            sum(na_cols_x)))
        }
        scaled_X <- scaled_X[, !na_cols_x, drop = FALSE]
    }

    ## Remove NA columns from Y
    na_cols_y <- apply(scaled_Y, 2, function(col) any(is.na(col)))
    if (any(na_cols_y)) {
        if (verbose) {
            message(sprintf("  Removing %d zero-variance gene columns from Y",
                            sum(na_cols_y)))
        }
        scaled_Y <- scaled_Y[, !na_cols_y, drop = FALSE]
    }

    ## Align genes between X and Y after column removal
    ## After bidirectional scaling, rows = features (TFs or samples),
    ## cols = genes. We need to keep common genes (columns).
    common_scaled <- intersect(colnames(scaled_X), colnames(scaled_Y))
    if (length(common_scaled) < n_cc) {
        stop("Insufficient genes remaining after scaling and NA removal")
    }
    scaled_X <- scaled_X[, common_scaled, drop = FALSE]
    scaled_Y <- scaled_Y[, common_scaled, drop = FALSE]

    if (verbose) {
        message(sprintf(
            "  After scaling: %d TFs x %d genes, %d samples x %d genes",
            nrow(scaled_X), ncol(scaled_X),
            nrow(scaled_Y), ncol(scaled_Y)))
    }

    ## ---- Step 4: SVD ----
    if (verbose) message("Computing SVD...")

    ## X is TF x gene, Y is sample x gene
    ## Cross-product: X %*% t(Y) -> but original code does X_t %*% Y
    ## In original code: XtY = scaled_X3_t %*% scaled_Y3
    ##   where scaled_X3_t = TF x gene, scaled_Y3 = gene x sample
    ## So XtY = TF x sample
    XtY <- scaled_X %*% t(scaled_Y)

    svd_res <- svd(XtY)

    ## Extract top n_cc components
    n_cc <- min(n_cc, length(svd_res$d))
    u <- svd_res$u[, seq_len(n_cc), drop = FALSE]
    v <- svd_res$v[, seq_len(n_cc), drop = FALSE]
    d <- svd_res$d

    ## Build TF weights DataFrame
    cc_names <- paste0("CC", seq_len(n_cc), "_rotation")
    tf_sample_loading <- S4Vectors::DataFrame(u)
    colnames(tf_sample_loading) <- cc_names
    rownames(tf_sample_loading) <- rownames(scaled_X)

    ## Build sample weights DataFrame
    transcriptome_sample_loading <- S4Vectors::DataFrame(v)
    colnames(transcriptome_sample_loading) <- cc_names
    rownames(transcriptome_sample_loading) <- rownames(scaled_Y)

    ## ---- Step 5: Compute parameters ----
    if (verbose) message("Computing CCA parameters...")
    params <- .compute_cca_parameters(
        scaled_X = scaled_X,
        scaled_Y = scaled_Y,
        u = u, v = v, d = d,
        n_cc = n_cc
    )

    ## ---- Build sample metadata DataFrame ----
    if (is.null(sample_metadata)) {
        smeta <- S4Vectors::DataFrame(row.names = colnames(transcriptome))
    } else {
        smeta <- S4Vectors::DataFrame(sample_metadata)
    }

    ## ---- Construct result ----
    if (verbose) message("Done.")

    new("CATaNResult",
        tf_matrix       = X,
        transcriptome   = Y,
        hvg             = hvg,
        svd_result      = list(u = svd_res$u, d = svd_res$d, v = svd_res$v),
        n_cc            = n_cc,
        tf_sample_loading = tf_sample_loading,
        transcriptome_sample_loading  = transcriptome_sample_loading,
        parameters      = params,
        sample_metadata = smeta,
        call            = mc)
}


## ---- Internal helpers ----

#' Bidirectional z-score scaling
#'
#' @param mat Numeric matrix (gene x feature).
#' @return Numeric matrix (feature x gene) after scaling.
#' @noRd
.bidirectional_scale <- function(mat) {
    ## Step 1: row-wise z-score (across features for each gene)
    ## apply(mat, 1, scale_fn) returns feature x gene
    scaled_1 <- apply(mat, 1, function(x) {
        s <- sd(x)
        if (s == 0) return(rep(NA_real_, length(x)))
        (x - mean(x)) / s
    })
    ## scaled_1 is now feature x gene

    ## Step 2: column-wise z-score (across genes for each feature)
    ## = row-wise on scaled_1
    scaled_2 <- apply(scaled_1, 1, function(x) {
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
        (x - mean(x, na.rm = TRUE)) / s
    })
    ## scaled_2 is now gene x feature -> transpose back to feature x gene
    t(scaled_2)
}


#' Compute CCA summary parameters
#' @noRd
.compute_cca_parameters <- function(scaled_X, scaled_Y, u, v, d, n_cc) {

    cc_names <- paste0("CC", seq_len(n_cc))

    ## --- Transcriptome variance explained by each CC ---
    ## For each CC k, sum of squared correlations between each gene and
    ## the k-th sample weight vector, divided by number of genes.
    ## scaled_Y is sample x gene, v is sample x n_cc
    ## Gene projections = t(scaled_Y) %*% v -> gene x n_cc
    trans_var <- numeric(n_cc)
    n_genes_y <- ncol(scaled_Y)
    for (k in seq_len(n_cc)) {
        ss <- 0
        for (j in seq_len(n_genes_y)) {
            ss <- ss + cor(scaled_Y[, j], v[, k])^2
        }
        trans_var[k] <- ss / n_genes_y
    }

    ## --- TF variance explained by each CC ---
    tf_var <- numeric(n_cc)
    n_tfs <- ncol(scaled_X)  # scaled_X is TF x gene, but we need TF-wise
    ## Actually scaled_X is TF x gene; u is TF x n_cc
    ## For each CC k, sum of squared cor between each TF and k-th TF weight
    for (k in seq_len(n_cc)) {
        ss <- 0
        for (j in seq_len(ncol(scaled_X))) {
            ss <- ss + cor(scaled_X[, j], u[, k])^2
        }
        tf_var[k] <- ss / ncol(scaled_X)
    }

    ## --- Singular value proportions ---
    sv_prop <- d[seq_len(n_cc)] / sum(d)

    ## --- Squared Canonical Fraction (SCF) ---
    scf <- d[seq_len(n_cc)]^2 / sum(d^2)

    ## --- Canonical correlations ---
    ## cor(X' %*% u_k, Y' %*% v_k) for each CC k
    ## X_proj = u' %*% scaled_X  -> n_cc x gene
    ## Y_proj = scaled_Y %*% v   -> but we need gene-level:
    ##   scaled_Y is sample x gene; v is sample x n_cc
    ##   Y_proj_gene = t(scaled_Y) %*% v -> gene x n_cc (gene-level projection)
    ## X is TF x gene; u is TF x n_cc
    ##   X_proj = t(u) %*% scaled_X -> n_cc x gene
    X_proj <- t(u) %*% scaled_X        # n_cc x gene
    Y_proj <- t(scaled_Y) %*% v        # gene x n_cc

    can_cor <- numeric(n_cc)
    for (k in seq_len(n_cc)) {
        can_cor[k] <- cor(X_proj[k, ], Y_proj[, k])
    }

    ## Build DataFrame
    param_mat <- rbind(
        transcriptome_variance = trans_var,
        tf_variance            = tf_var,
        sv_proportion          = sv_prop,
        scf                    = scf,
        canonical_correlation  = can_cor
    )
    colnames(param_mat) <- cc_names
    S4Vectors::DataFrame(param_mat)
}

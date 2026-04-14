test_that("filter_low_expression works with basic matrix", {
    set.seed(1)
    ## Create a matrix with some low-expression genes
    counts <- matrix(rnbinom(2000, mu = 50, size = 2),
                     nrow = 200, ncol = 10)
    rownames(counts) <- paste0("gene", seq_len(200))
    ## Add some zero genes
    counts[1:10, ] <- 0

    result <- filter_low_expression(counts, verbose = FALSE)

    expect_true(is.matrix(result))
    expect_true(nrow(result) < nrow(counts))
    expect_true(nrow(result) > 0)
    ## Zero genes should be removed
    expect_false(any(paste0("gene", 1:10) %in% rownames(result)))
    ## Column count unchanged
    expect_equal(ncol(result), ncol(counts))
})

test_that("filter_low_expression rejects bad input", {
    expect_error(filter_low_expression("not a matrix"),
                 "numeric matrix")
    m <- matrix(1:10, nrow = 5)
    expect_error(filter_low_expression(m),
                 "row names")
})

test_that("prepare_matrices intersects genes and normalises", {
    set.seed(2)
    ngenes <- 100
    counts <- matrix(rnbinom(ngenes * 10, mu = 50, size = 2),
                     nrow = ngenes, ncol = 10)
    rownames(counts) <- paste0("gene", seq_len(ngenes))
    colnames(counts) <- paste0("sample", seq_len(10))

    ## TF matrix with partial overlap
    tf_genes <- paste0("gene", 51:150)  # 50 overlap with counts
    tf_mat <- matrix(sample(0:1, 100 * 5, replace = TRUE),
                     nrow = 100, ncol = 5)
    rownames(tf_mat) <- tf_genes
    colnames(tf_mat) <- paste0("TF", seq_len(5))

    result <- prepare_matrices(counts, tf_mat, verbose = FALSE)

    expect_type(result, "list")
    expect_named(result, c("tf", "expr"))
    expect_equal(nrow(result$tf), nrow(result$expr))
    expect_equal(nrow(result$tf), 50)  # 50 common genes
    expect_equal(rownames(result$tf), rownames(result$expr))
})

test_that("run_cca produces valid CATaNResult", {
    set.seed(42)
    ngenes <- 200
    ntf <- 20
    nsamp <- 30

    tf_mat <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
                     nrow = ngenes, ncol = ntf)
    rownames(tf_mat) <- paste0("gene", seq_len(ngenes))
    colnames(tf_mat) <- paste0("TF", seq_len(ntf))

    expr_mat <- matrix(rnorm(ngenes * nsamp, mean = 5, sd = 2),
                       nrow = ngenes, ncol = nsamp)
    rownames(expr_mat) <- paste0("gene", seq_len(ngenes))
    colnames(expr_mat) <- paste0("sample", seq_len(nsamp))

    res <- run_cca(tf_mat, expr_mat, n_cc = 5L, n_hvg = 100L, verbose = FALSE)

    ## Check class
    expect_s4_class(res, "CATaNResult")

    ## Check slots
    expect_equal(res@n_cc, 5L)
    expect_equal(ncol(res@tf_weights), 5)
    expect_equal(ncol(res@sample_weights), 5)
    expect_true(length(res@hvg) > 0)

    ## Check SVD components
    expect_named(res@svd_result, c("u", "d", "v"))
    expect_true(all(res@svd_result$d >= 0))  # singular values non-negative

    ## Check parameters
    params <- as.data.frame(res@parameters)
    expect_true("canonical_correlation" %in% rownames(params))
    expect_true("scf" %in% rownames(params))
    expect_true("transcriptome_variance" %in% rownames(params))
    expect_true("tf_variance" %in% rownames(params))

    ## SCF should sum to <= 1
    scf_sum <- sum(as.numeric(params["scf", ]))
    expect_true(scf_sum <= 1)
    expect_true(scf_sum > 0)

    ## Canonical correlations should be between -1 and 1
    cancor <- as.numeric(params["canonical_correlation", ])
    expect_true(all(abs(cancor) <= 1))

    ## Validity check
    expect_true(validObject(res))
})

test_that("run_cca handles n_cc larger than possible", {
    set.seed(10)
    ngenes <- 50
    ntf <- 3
    nsamp <- 5

    tf_mat <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
                     nrow = ngenes, ncol = ntf)
    rownames(tf_mat) <- paste0("gene", seq_len(ngenes))
    colnames(tf_mat) <- paste0("TF", seq_len(ntf))

    expr_mat <- matrix(rnorm(ngenes * nsamp), nrow = ngenes, ncol = nsamp)
    rownames(expr_mat) <- paste0("gene", seq_len(ngenes))
    colnames(expr_mat) <- paste0("sample", seq_len(nsamp))

    ## n_cc = 3 should work (min of TF and sample dimensions)
    res <- run_cca(tf_mat, expr_mat, n_cc = 3L, n_hvg = 30L, verbose = FALSE)
    expect_s4_class(res, "CATaNResult")
})

test_that("run_cca fails gracefully with no common genes", {
    tf_mat <- matrix(1, nrow = 5, ncol = 3)
    rownames(tf_mat) <- paste0("geneA", 1:5)
    colnames(tf_mat) <- paste0("TF", 1:3)

    expr_mat <- matrix(1, nrow = 5, ncol = 3)
    rownames(expr_mat) <- paste0("geneB", 1:5)
    colnames(expr_mat) <- paste0("s", 1:3)

    expect_error(run_cca(tf_mat, expr_mat, verbose = FALSE),
                 "No common genes")
})

test_that("bidirectional scaling produces expected dimensions", {
    mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(mat) <- paste0("gene", 1:10)
    colnames(mat) <- paste0("feat", 1:10)

    result <- CATaN:::.bidirectional_scale(mat)

    ## After bidirectional scaling: result is feature x gene
    expect_equal(nrow(result), ncol(mat))  # features
    expect_equal(ncol(result), nrow(mat))  # genes
})

test_that("show methods work without error", {
    res <- new("CATaNResult")
    expect_output(show(res), "CATaNResult")

    annot <- new("CATaNAnnotation")
    expect_output(show(annot), "CATaNAnnotation")
})

test_that("accessors work on CATaNResult", {
    set.seed(42)
    ngenes <- 100
    ntf <- 10
    nsamp <- 15

    tf_mat <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
                     nrow = ngenes, ncol = ntf)
    rownames(tf_mat) <- paste0("gene", seq_len(ngenes))
    colnames(tf_mat) <- paste0("TF", seq_len(ntf))

    expr_mat <- matrix(rnorm(ngenes * nsamp, mean = 5),
                       nrow = ngenes, ncol = nsamp)
    rownames(expr_mat) <- paste0("gene", seq_len(ngenes))
    colnames(expr_mat) <- paste0("sample", seq_len(nsamp))

    res <- run_cca(tf_mat, expr_mat, n_cc = 3L, n_hvg = 50L, verbose = FALSE)

    expect_s4_class(tfWeights(res), "DataFrame")
    expect_s4_class(sampleWeights(res), "DataFrame")
    expect_s4_class(ccParameters(res), "DataFrame")
    expect_equal(ncol(tfWeights(res)), 3)
    expect_equal(ncol(sampleWeights(res)), 3)
})

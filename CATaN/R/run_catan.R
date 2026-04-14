#' Run the full CATaN pipeline
#'
#' End-to-end wrapper that runs CCA, aligns CC scores to SNPs, extracts
#' top/bottom SNPs, and optionally exports sLDSC-ready files.
#'
#' @param counts Numeric matrix of raw counts (gene x sample) or a
#'   SummarizedExperiment.
#' @param tf_matrix Numeric matrix of TF-gene binding (gene x TF).
#' @param peak_dir Character. Path to TF peak BED directory. If NULL,
#'   \code{\link{download_peak_beds}} is called automatically.
#' @param snp_gr GRanges of target SNP positions.
#' @param sample_metadata data.frame of sample annotations (optional).
#' @param n_cc Integer. Number of CCs (default 10L).
#' @param n_hvg Integer. Number of HVGs (default 10000L).
#' @param percentile Numeric. Top/bottom fraction (default 0.1).
#' @param population Character (default "custom").
#' @param output_dir Character. If non-NULL, export sLDSC files here.
#' @param BPPARAM BiocParallelParam (default bpparam()).
#' @param verbose Logical (default TRUE).
#'
#' @return A \code{\linkS4class{CATaNAnnotation}} object.
#'
#' @examples
#' set.seed(1)
#' ngenes <- 100; ntf <- 10; nsamp <- 15
#' tf <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
#'     ngenes, ntf)
#' rownames(tf) <- paste0("g", seq_len(ngenes))
#' colnames(tf) <- paste0("TF", seq_len(ntf))
#' cts <- matrix(rnbinom(ngenes * nsamp, mu = 50, size = 2),
#'     ngenes, nsamp)
#' rownames(cts) <- rownames(tf)
#' colnames(cts) <- paste0("s", seq_len(nsamp))
#'
#' # Create mock peak BED files
#' peak_dir <- tempfile("peaks")
#' dir.create(peak_dir)
#' for (t in colnames(tf)) {
#'     bed <- data.frame("chr1", c(100, 500), c(200, 600))
#'     write.table(bed,
#'         file.path(peak_dir, paste0(t, ".hg19.bed")),
#'         sep = "\t", row.names = FALSE,
#'         col.names = FALSE, quote = FALSE)
#' }
#'
#' # Create mock SNPs
#' library(GenomicRanges)
#' snps <- GRanges("chr1",
#'     IRanges::IRanges(c(150, 550, 9999), width = 1))
#'
#' result <- run_catan(cts, tf, peak_dir = peak_dir,
#'     snp_gr = snps, n_cc = 2L, n_hvg = 50L,
#'     BPPARAM = BiocParallel::SerialParam(),
#'     verbose = FALSE)
#'
#' @export
run_catan <- function(counts,
                      tf_matrix,
                      peak_dir        = NULL,
                      snp_gr,
                      sample_metadata = NULL,
                      n_cc            = 10L,
                      n_hvg           = 10000L,
                      percentile      = 0.1,
                      population      = "custom",
                      output_dir      = NULL,
                      BPPARAM         = BiocParallel::bpparam(),
                      verbose         = TRUE) {

    ## Step 1: Preprocess
    if (verbose) message("=== Step 1: Preprocessing ===")
    filtered <- filter_low_expression(counts, verbose = verbose)
    mats     <- prepare_matrices(filtered, tf_matrix, verbose = verbose)

    ## Step 2: CCA
    if (verbose) message("=== Step 2: CCA ===")
    cca_res <- run_cca(
        tf_matrix       = mats$tf,
        transcriptome   = mats$expr,
        n_cc            = n_cc,
        n_hvg           = n_hvg,
        sample_metadata = sample_metadata,
        verbose         = verbose
    )

    ## Step 3: Align to SNPs
    if (verbose) message("=== Step 3: SNP annotation ===")
    if (is.null(peak_dir)) {
        peak_dir <- download_peak_beds(
            dest_dir = tempdir(),
            verbose  = verbose
        )
    }
    annot <- align_cc_to_snps(
        catan_result = cca_res,
        peak_dir     = peak_dir,
        snp_gr       = snp_gr,
        population   = population,
        BPPARAM      = BPPARAM,
        verbose      = verbose
    )

    ## Step 4: Extract top/bottom
    annot <- extract_top_bottom_snps(annot, percentile = percentile)

    ## Step 5: Export
    if (!is.null(output_dir)) {
        if (verbose) message("=== Step 4: Exporting sLDSC files ===")
        export_for_sldsc(annot, output_dir = output_dir)
    }

    annot
}

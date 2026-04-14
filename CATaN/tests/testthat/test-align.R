## ---- Helper: create mock data for Phase 2 tests ----

.make_test_cca_result <- function() {
    set.seed(42)
    ngenes <- 100; ntf <- 5; nsamp <- 10

    tf_mat <- matrix(sample(0:1, ngenes * ntf, replace = TRUE),
                     nrow = ngenes, ncol = ntf)
    rownames(tf_mat) <- paste0("gene", seq_len(ngenes))
    colnames(tf_mat) <- paste0("TF", seq_len(ntf))

    expr_mat <- matrix(rnorm(ngenes * nsamp, mean = 5, sd = 2),
                       nrow = ngenes, ncol = nsamp)
    rownames(expr_mat) <- paste0("gene", seq_len(ngenes))
    colnames(expr_mat) <- paste0("sample", seq_len(nsamp))

    run_cca(tf_mat, expr_mat, n_cc = 3L, n_hvg = 50L, verbose = FALSE)
}

.make_test_peaks <- function(tmp_dir) {
    ## Create simple BED files for 5 TFs
    for (i in 1:5) {
        tf_name <- paste0("TF", i)
        bed_content <- data.frame(
            chr   = rep("chr1", 10),
            start = seq(100, 1000, by = 100) + (i - 1) * 10,
            end   = seq(200, 1100, by = 100) + (i - 1) * 10
        )
        ## Add some chr2 peaks
        bed_chr2 <- data.frame(
            chr   = rep("chr2", 5),
            start = seq(500, 900, by = 100),
            end   = seq(600, 1000, by = 100)
        )
        bed_all <- rbind(bed_content, bed_chr2)
        write.table(bed_all,
                    file.path(tmp_dir, paste0(tf_name, ".hg19.bed")),
                    sep = "\t", row.names = FALSE, col.names = FALSE,
                    quote = FALSE)
    }
}

.make_test_snps <- function() {
    GenomicRanges::GRanges(
        seqnames = c(rep("chr1", 20), rep("chr2", 10)),
        ranges   = IRanges::IRanges(
            start = c(seq(150, 1050, length.out = 20),
                      seq(550, 950, length.out = 10)),
            width = 1
        )
    )
}


## ---- Tests for align_cc_to_snps ----

test_that("align_cc_to_snps works with test data", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    .make_test_peaks(tmp_dir)
    snp_gr <- .make_test_snps()

    annot <- align_cc_to_snps(
        catan_result = cca_res,
        peak_dir     = tmp_dir,
        snp_gr       = snp_gr,
        chromosomes  = c("chr1", "chr2"),
        BPPARAM      = BiocParallel::SerialParam(),
        verbose      = FALSE
    )

    expect_s4_class(annot, "CATaNAnnotation")
    expect_equal(length(annot@snp_scores), 3)  # 3 CCs

    ## Each CC should have all 30 SNPs
    for (cc in names(annot@snp_scores)) {
        expect_equal(length(annot@snp_scores[[cc]]), 30)
    }

    ## Some SNPs should have non-zero scores
    scores_cc1 <- GenomicRanges::mcols(annot@snp_scores[[1]])$score
    expect_true(any(scores_cc1 != 0))
})

test_that("SNPs outside peaks receive score zero", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    .make_test_peaks(tmp_dir)

    ## Create SNPs that are far outside any peak region
    snp_gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1"),
        ranges   = IRanges::IRanges(start = c(150, 50000, 99999), width = 1)
    )

    annot <- align_cc_to_snps(
        cca_res, tmp_dir, snp_gr,
        chromosomes = "chr1",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    scores <- GenomicRanges::mcols(annot@snp_scores[[1]])$score
    ## SNP at 150 should overlap peaks, but 50000 and 99999 should not
    expect_true(scores[1] != 0)  # inside peak
    expect_equal(scores[2], 0)   # outside peak
    expect_equal(scores[3], 0)   # outside peak
})

test_that("align_cc_to_snps errors on bad input", {
    expect_error(align_cc_to_snps("not_a_result", ".", GenomicRanges::GRanges()),
                 "CATaNResult")

    cca_res <- .make_test_cca_result()
    expect_error(align_cc_to_snps(cca_res, "/nonexistent/path",
                                  GenomicRanges::GRanges()),
                 "does not exist")
})

test_that("align_cc_to_snps handles missing BED files gracefully", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_partial_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    ## Only create 2 of 5 TF BED files
    for (i in 1:2) {
        bed <- data.frame(chr = "chr1", start = i * 100, end = i * 100 + 50)
        write.table(bed,
                    file.path(tmp_dir, paste0("TF", i, ".hg19.bed")),
                    sep = "\t", row.names = FALSE, col.names = FALSE,
                    quote = FALSE)
    }

    snp_gr <- .make_test_snps()

    ## Should work with warning about missing TFs
    annot <- align_cc_to_snps(
        cca_res, tmp_dir, snp_gr,
        chromosomes = c("chr1", "chr2"),
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    expect_s4_class(annot, "CATaNAnnotation")
})


## ---- Tests for extract_top_bottom_snps ----

test_that("extract_top_bottom_snps works correctly", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    .make_test_peaks(tmp_dir)
    snp_gr <- .make_test_snps()

    annot <- align_cc_to_snps(
        cca_res, tmp_dir, snp_gr,
        chromosomes = c("chr1", "chr2"),
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    result <- extract_top_bottom_snps(annot, percentile = 0.1)

    expect_s4_class(result, "CATaNAnnotation")
    expect_equal(length(result@top_snps), 3)
    expect_equal(length(result@bottom_snps), 3)
    expect_equal(result@percentile, 0.1)

    ## Check counts: 30 SNPs, 10% = 3 per tail
    for (cc in names(result@top_snps)) {
        expect_equal(length(result@top_snps[[cc]]), 3)
        expect_equal(length(result@bottom_snps[[cc]]), 3)
    }

    ## Top scores should be >= bottom scores
    for (cc in names(result@top_snps)) {
        top_scores <- GenomicRanges::mcols(result@top_snps[[cc]])$score
        bot_scores <- GenomicRanges::mcols(result@bottom_snps[[cc]])$score
        expect_true(min(top_scores) >= max(bot_scores))
    }
})

test_that("extract_top_bottom_snps rejects bad input", {
    expect_error(extract_top_bottom_snps("not_annotation"),
                 "CATaNAnnotation")

    annot <- new("CATaNAnnotation")
    expect_error(extract_top_bottom_snps(annot),
                 "snp_scores is empty")
})

test_that("extract_top_bottom_snps with different percentiles", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))

    .make_test_peaks(tmp_dir)
    snp_gr <- .make_test_snps()

    annot <- align_cc_to_snps(
        cca_res, tmp_dir, snp_gr,
        chromosomes = c("chr1", "chr2"),
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    ## 20% = 6 SNPs per tail
    result_20 <- extract_top_bottom_snps(annot, percentile = 0.2)
    for (cc in names(result_20@top_snps)) {
        expect_equal(length(result_20@top_snps[[cc]]), 6)
    }
})


## ---- Tests for export_for_sldsc ----

test_that("export_for_sldsc writes correct files", {
    cca_res <- .make_test_cca_result()
    tmp_dir <- tempfile("peaks_")
    dir.create(tmp_dir)
    .make_test_peaks(tmp_dir)
    snp_gr <- .make_test_snps()

    annot <- align_cc_to_snps(
        cca_res, tmp_dir, snp_gr,
        chromosomes = c("chr1", "chr2"),
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    annot <- extract_top_bottom_snps(annot, percentile = 0.1)

    out_dir <- tempfile("sldsc_out_")
    on.exit(unlink(c(tmp_dir, out_dir), recursive = TRUE))

    files <- export_for_sldsc(annot, out_dir, prefix = "test")

    ## Should create 6 files (3 CCs x 2 top/bottom)
    expect_equal(length(files), 6)
    expect_true(all(file.exists(files)))

    ## Check file naming
    expect_true(any(grepl("top\\.bed\\.gz$", files)))
    expect_true(any(grepl("bottom\\.bed\\.gz$", files)))

    ## Check content of one file
    top_file <- files[grep("CC1_rotation_top", files)]
    bed <- read.table(gzfile(top_file), sep = "\t", header = FALSE)
    expect_equal(ncol(bed), 4)  # chr, start, end, score
    expect_equal(nrow(bed), 3)  # 10% of 30 = 3

    ## Verify 0-based coordinates (start < end)
    expect_true(all(bed[, 2] < bed[, 3]))
})

test_that("export_for_sldsc errors without top/bottom", {
    annot <- new("CATaNAnnotation")
    expect_error(export_for_sldsc(annot, tempdir()),
                 "top_snps.*empty")
})

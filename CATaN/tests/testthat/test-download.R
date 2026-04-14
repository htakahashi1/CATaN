test_that("download_peak_beds errors when URL not configured", {
    expect_error(download_peak_beds(tempdir()),
                 "Zenodo URL has not been configured")
})

test_that("download_tf_matrix errors when URL not configured", {
    expect_error(download_tf_matrix(tempdir()),
                 "Zenodo URL has not been configured")
})

test_that("download_peak_beds works with custom URL", {
    ## Create a mock directory with fake BED files directly
    ## and test the extraction verification logic
    tmp_src <- tempfile("mock_beds_")
    dir.create(tmp_src)
    on.exit(unlink(tmp_src, recursive = TRUE))

    ## Create fake BED files
    for (tf in c("TF1", "TF2")) {
        bed <- data.frame(chr = "chr1", start = 100, end = 200)
        write.table(bed,
                    file.path(tmp_src, paste0(tf, ".hg19.bed")),
                    sep = "\t", row.names = FALSE, col.names = FALSE,
                    quote = FALSE)
    }

    ## Create tar.gz using setwd for relative paths
    tar_path <- tempfile(fileext = ".tar.gz")
    old_wd <- setwd(tmp_src)
    on.exit(setwd(old_wd), add = TRUE, after = FALSE)
    utils::tar(tar_path, files = list.files("."),
               compression = "gzip", tar = "internal")
    setwd(old_wd)

    ## Test: use local path
    dest <- tempfile("peaks_dest_")
    result <- download_peak_beds(dest, url = tar_path, cache = FALSE,
                                 verbose = FALSE)

    expect_true(dir.exists(result))
    bed_files <- list.files(result, pattern = "\\.hg19\\.bed$")
    expect_equal(length(bed_files), 2)
})

test_that("catan_data_status returns correct structure", {
    status <- catan_data_status()
    expect_true(is.data.frame(status))
    expect_true("name" %in% colnames(status))
    expect_true("cached" %in% colnames(status))
    expect_equal(nrow(status), 2)
})

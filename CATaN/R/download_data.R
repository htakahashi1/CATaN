## Download CATaN reference data from Zenodo

## ---- Default URLs (update after Zenodo upload) ----


.CATAN_ZENODO_PEAK_URL <- paste0(
    "https://zenodo.org/records/19507630/files/",
    "TF_ChIPseq_bed_hg19.tar.gz"
)
.CATAN_ZENODO_MATRIX_URL <- paste0(
    "https://zenodo.org/records/19507630/files/",
    "TF_GRN_matrix.txt"
)


#' Download TF peak BED files from Zenodo
#'
#' Downloads the TF ChIP-seq peak BED files (hg19) used in the CATaN paper.
#' The archive is cached via BiocFileCache and extracted to the specified
#' directory.
#'
#' @param dest_dir Character. Directory to extract BED files into. Created
#'   if it does not exist.
#' @param url Character. URL of the tar.gz archive. If NULL (default), uses
#'   the default Zenodo URL.
#' @param cache Logical. Whether to use BiocFileCache for caching
#'   (default TRUE).
#' @param verbose Logical (default TRUE).
#'
#' @return Character. Path to the directory containing the extracted BED files.
#'
#' @details
#' On first call, the archive (~1 GB) is downloaded and cached. Subsequent
#' calls skip the download and use the cached copy. The archive is extracted
#' to \code{dest_dir}, which can then be passed to
#' \code{\link{align_cc_to_snps}} as the \code{peak_dir} argument.
#'
#' @examples
#' \donttest{
#' peak_dir <- download_peak_beds("~/catan_data/peaks")
#' }
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath bfcnew
#' @importFrom utils download.file untar
#' @export
download_peak_beds <- function(dest_dir,
                               url     = NULL,
                               cache   = TRUE,
                               verbose = TRUE) {

    if (is.null(url)) {
        url <- .CATAN_ZENODO_PEAK_URL
    }


    if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
    }

    ## Download (with cache)
    archive_path <- .cached_download(url, "TF_ChIPseq_bed_hg19.tar.gz",
                                     cache = cache, verbose = verbose)

    ## Extract
    if (verbose) message("Extracting BED files to ", dest_dir, "...")
    utils::untar(archive_path, exdir = dest_dir)

    ## Verify
    bed_files <- list.files(dest_dir, pattern = "\\.hg19\\.bed$")
    if (length(bed_files) == 0) {
        warning("No .hg19.bed files found after extraction. ",
                "Check archive contents.")
    } else if (verbose) {
        message(sprintf("  %d TF peak BED files available", length(bed_files)))
    }

    dest_dir
}


#' Download TF-gene binding matrix from Zenodo
#'
#' Downloads the TF-gene binding association matrix used in the CATaN paper.
#'
#' @param dest_dir Character. Directory to save the file. Created if it does
#'   not exist.
#' @param url Character. URL of the gzipped matrix file. If NULL (default),
#'   uses the default Zenodo URL.
#' @param cache Logical. Whether to use BiocFileCache (default TRUE).
#' @param verbose Logical (default TRUE).
#'
#' @return Character. Path to the downloaded (uncompressed) matrix file.
#'
#' @examples
#' \donttest{
#' matrix_file <- download_tf_matrix("~/catan_data")
#' tf_mat <- read.table(matrix_file, header = TRUE, row.names = 1)
#' }
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath
#' @export
download_tf_matrix <- function(dest_dir,
                               url     = NULL,
                               cache   = TRUE,
                               verbose = TRUE) {

    if (is.null(url)) {
        url <- .CATAN_ZENODO_MATRIX_URL
    }



    if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
    }

    ## Download (with cache)
    gz_path <- .cached_download(url, "TF_GRN_matrix.txt",
                                cache = cache, verbose = verbose)

    ## Copy to dest_dir
    out_path <- file.path(dest_dir, "TF_GRN_matrix.txt")
    if (!file.exists(out_path)) {
      if (verbose) message("Copying to ", out_path, "...")
      file.copy(gz_path, out_path)
    } else if (verbose) {
      message("TF matrix already exists at ", out_path)
    }

    out_path
}


#' Internal: download file with BiocFileCache
#' @noRd
.cached_download <- function(url, fname, cache = TRUE, verbose = TRUE) {

    if (cache) {
        bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

        ## Check if already cached
        res <- BiocFileCache::bfcquery(bfc, fname, exact = TRUE)

        if (nrow(res) > 0) {
            cached_path <- BiocFileCache::bfcrpath(bfc, rids = res$rid[1])
            if (verbose) message("Using cached file: ", cached_path)
            return(cached_path)
        }

        ## Download and cache
        if (verbose) message("Downloading ", fname, " from Zenodo...")
        cached_path <- BiocFileCache::bfcadd(bfc, rname = fname, fpath = url)
        if (verbose) message("  Cached at: ", cached_path)
        return(cached_path)

    } else {
        ## Direct download to temp
        tmp_path <- file.path(tempdir(), fname)
        if (!file.exists(tmp_path)) {
            if (verbose) message("Downloading ", fname, "...")
            if (file.exists(url)) {
                ## Local file path
                file.copy(url, tmp_path)
            } else {
                utils::download.file(
                    url, tmp_path,
                    mode = "wb", quiet = !verbose)
            }
        }
        return(tmp_path)
    }
}


#' List available CATaN reference data
#'
#' Shows the status of locally cached CATaN data files.
#'
#' @return A data.frame with columns: name, cached, path.
#'
#' @examples
#' status <- catan_data_status()
#' status
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery
#' @export
catan_data_status <- function() {
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

    files <- c(
        peaks  = "TF_ChIPseq_bed_hg19.tar.gz",
        matrix = "TF_GRN_matrix.txt"
    )

    status <- data.frame(
        name   = names(files),
        file   = unname(files),
        cached = FALSE,
        path   = NA_character_,
        stringsAsFactors = FALSE
    )

    for (i in seq_along(files)) {
        res <- BiocFileCache::bfcquery(bfc, files[i], exact = TRUE)
        if (nrow(res) > 0) {
            status$cached[i] <- TRUE
            status$path[i]   <- BiocFileCache::bfcrpath(bfc, rids = res$rid[1])
        }
    }

    status
}

#' Export CC annotations as S-LDSC-ready BED files
#'
#' Writes top and bottom SNP sets for each CC as gzipped BED files
#' suitable for use as annotations in stratified LD Score Regression.
#' Output files are tab-delimited BED4 format (chr, start, end, score)
#' sorted by coordinate.
#'
#' @param annotation A \code{\linkS4class{CATaNAnnotation}} object with
#'   \code{top_snps} and \code{bottom_snps} populated (output of
#'   \code{\link{extract_top_bottom_snps}}).
#' @param output_dir Character. Directory to write output files. Created
#'   if it does not exist.
#' @param prefix Character. File name prefix (default \code{"CATaN"}).
#'
#' @return Invisibly returns a character vector of output
#'   file paths.
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges::IRanges(seq(1, 100), width = 1))
#' mcols(gr)$score <- rnorm(100)
#' annot <- new("CATaNAnnotation",
#'     snp_scores = list(CC1_rotation = gr),
#'     percentile = 0.1, population = "test")
#' annot <- extract_top_bottom_snps(annot)
#' out <- export_for_sldsc(annot, tempdir(), prefix = "ex")
#'
#' @details
#' For each CC, two files are created:
#' \itemize{
#'   \item \code{{prefix}_{CC}_top.bed.gz} — SNPs in the top percentile
#'   \item \code{{prefix}_{CC}_bottom.bed.gz} — SNPs in the bottom percentile
#' }
#'
#' The BED format uses 0-based half-open coordinates, matching the output
#' of the original bedtools-based pipeline.
#'
#' @importFrom GenomicRanges mcols seqnames start end
#' @importFrom methods is
#' @importFrom utils write.table
#' @export
export_for_sldsc <- function(annotation, output_dir, prefix = "CATaN") {

    ## ---- Input validation ----
    if (!is(annotation, "CATaNAnnotation")) {
        stop("annotation must be a CATaNAnnotation object")
    }
    if (length(annotation@top_snps) == 0 ||
        length(annotation@bottom_snps) == 0) {
        stop("top_snps and/or bottom_snps are empty. ",
             "Run extract_top_bottom_snps() first.")
    }

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    output_files <- character()

    for (cc_name in names(annotation@top_snps)) {

        ## ---- Top SNPs ----
        top_gr <- annotation@top_snps[[cc_name]]
        top_file <- file.path(output_dir,
                              paste0(prefix, "_", cc_name, "_top.bed.gz"))
        .write_bed_gz(top_gr, top_file)
        output_files <- c(output_files, top_file)

        ## ---- Bottom SNPs ----
        bot_gr <- annotation@bottom_snps[[cc_name]]
        bot_file <- file.path(output_dir,
                              paste0(prefix, "_", cc_name, "_bottom.bed.gz"))
        .write_bed_gz(bot_gr, bot_file)
        output_files <- c(output_files, bot_file)
    }

    message(sprintf("Exported %d files to %s",
                    length(output_files), output_dir))
    invisible(output_files)
}


#' Write a GRanges to a gzipped BED4 file
#' @noRd
.write_bed_gz <- function(gr, filepath) {
    if (length(gr) == 0) {
        ## Write empty file
        con <- gzfile(filepath, "wb")
        close(con)
        return(invisible(filepath))
    }

    ## Build BED4 data.frame (0-based start, 1-based end)
    bed_df <- data.frame(
        chr   = as.character(GenomicRanges::seqnames(gr)),
        start = GenomicRanges::start(gr) - 1L,
        end   = GenomicRanges::end(gr),
        score = GenomicRanges::mcols(gr)$score,
        stringsAsFactors = FALSE
    )

    ## Sort by coordinate
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]

    ## Write gzipped
    con <- gzfile(filepath, "wt")
    utils::write.table(bed_df, con, sep = "\t", row.names = FALSE,
                       col.names = FALSE, quote = FALSE)
    close(con)

    invisible(filepath)
}

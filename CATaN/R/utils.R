## Internal utility functions for CATaN


#' Read a PLINK .bim file as GRanges
#'
#' Converts a PLINK .bim file to a GRanges object with hg19 coordinates.
#' This is a convenience function for users who have .bim files from
#' 1000 Genomes Phase 3.
#'
#' @param bim_file Path to a .bim file.
#' @return A GRanges object with SNP positions (1-based, width 1).
#'
#' @examples
#' # Create a minimal mock .bim file
#' tmp <- tempfile(fileext = ".bim")
#' bim <- data.frame(1, c("rs1", "rs2"), 0,
#'     c(1000, 2000), "A", "G")
#' write.table(bim, tmp, sep = "\t",
#'     row.names = FALSE, col.names = FALSE, quote = FALSE)
#' gr <- bim_to_granges(tmp)
#' gr
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom utils read.table
#' @export
bim_to_granges <- function(bim_file) {
    if (!file.exists(bim_file)) {
        stop("File not found: ", bim_file)
    }

    bim <- utils::read.table(bim_file, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)
    colnames(bim) <- c("chr", "id", "cm", "pos", "a1", "a2")

    gr <- GenomicRanges::GRanges(
        seqnames = bim$chr,
        ranges   = IRanges::IRanges(start = bim$pos, end = bim$pos),
        snp_id   = bim$id
    )
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    GenomeInfoDb::seqlevels(gr) <- paste0("chr", c(1:22))
    sort(gr)
}


#' Validate genome build consistency
#' @noRd
.check_genome <- function(..., expected = "hg19") {
    grs <- list(...)
    for (i in seq_along(grs)) {
        if (is(grs[[i]], "GRanges")) {
            g <- unique(GenomeInfoDb::genome(grs[[i]]))
            g <- g[!is.na(g)]
            if (length(g) > 0 && !all(g == expected)) {
                warning(sprintf(
                    "GRanges %d has genome '%s', expected '%s'",
                                i, paste(g, collapse = ","), expected))
            }
        }
    }
    invisible(TRUE)
}

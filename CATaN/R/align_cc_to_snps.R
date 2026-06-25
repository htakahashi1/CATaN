#' Align CC scores to SNP positions
#'
#' Maps TF sample loadings from a CATaNResult to genomic SNP positions using
#' TF ChIP-seq peak BED files. For each CC, every TF's peak regions are
#' assigned that TF's sample loadings, overlapped with the target SNP set, and
#' aggregated by summing scores at each SNP. SNPs not overlapping any peak
#' receive a score of zero.
#'
#' This replicates the logic of the original bedtools-based shell script
#' using GenomicRanges for full R compatibility.
#'
#' @param catan_result A \code{\linkS4class{CATaNResult}} object.
#' @param peak_dir Character. Path to directory containing TF peak BED files
#'   (one per TF, named \code{{TF_name}.hg19.bed}).
#' @param snp_gr A \code{GRanges} object of target SNP positions. Can be
#'   created from a .bim file using \code{\link{bim_to_granges}}.
#' @param chromosomes Character vector of chromosome names to process
#'   (default \code{paste0("chr", 1:22)}).
#' @param population Character. Label for the reference population
#'   (default \code{"custom"}).
#' @param BPPARAM A \code{BiocParallelParam} object for parallel processing
#'   (default \code{BiocParallel::bpparam()}).
#' @param verbose Logical (default TRUE).
#'
#' @return A \code{\linkS4class{CATaNAnnotation}} object with the
#'   \code{snp_scores} slot populated.
#'
#' @examples
#' set.seed(1)
#' tf <- matrix(sample(0:1, 200, replace = TRUE), 20, 10)
#' rownames(tf) <- paste0("g", seq_len(20))
#' colnames(tf) <- paste0("TF", seq_len(10))
#' ex <- matrix(rnorm(100, 5), 20, 5)
#' rownames(ex) <- rownames(tf)
#' colnames(ex) <- paste0("s", seq_len(5))
#' res <- run_cca(tf, ex, n_cc = 2L, n_hvg = 15L,
#'                verbose = FALSE)
#'
#' # Create mock peak BED files
#' peak_dir <- tempfile("peaks")
#' dir.create(peak_dir)
#' for (tf_name in colnames(tf)) {
#'     bed <- data.frame("chr1", c(100, 500), c(200, 600))
#'     write.table(bed,
#'         file.path(peak_dir, paste0(tf_name, ".hg19.bed")),
#'         sep = "\t", row.names = FALSE,
#'         col.names = FALSE, quote = FALSE)
#' }
#'
#' # Create mock SNPs
#' library(GenomicRanges)
#' snps <- GRanges("chr1",
#'     IRanges::IRanges(c(150, 550, 9999), width = 1))
#' annot <- align_cc_to_snps(res, peak_dir, snps,
#'     chromosomes = "chr1",
#'     BPPARAM = BiocParallel::SerialParam(),
#'     verbose = FALSE)
#'
#' @details
#' The algorithm proceeds per CC, per chromosome:
#' \enumerate{
#'   \item Load each TF's peak BED file and assign that TF sample loading as
#'     the score.
#'   \item Combine all TF peaks and find overlaps with the target SNP set
#'     (equivalent to \code{bedtools intersect}).
#'   \item For SNPs overlapping multiple TF peaks, sum the scores
#'     (equivalent to \code{bedtools merge -d -1 -c 4 -o sum}).
#'   \item Assign score 0 to SNPs with no peak overlap
#'     (equivalent to \code{bedtools intersect -v} + zero fill).
#'   \item Concatenate across chromosomes.
#' }
#'
#' Chromosomes are processed in parallel via \code{BiocParallel}.
#'
#' @importFrom GenomicRanges GRanges mcols mcols<- findOverlaps seqnames
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle keepSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom rtracklayer import.bed
#' @importFrom methods new is
#' @export
align_cc_to_snps <- function(catan_result,
                             peak_dir,
                             snp_gr,
                             chromosomes = paste0("chr", 1:22),
                             population  = "custom",
                             BPPARAM     = BiocParallel::bpparam(),
                             verbose     = TRUE) {

    ## ---- Input validation ----
    if (!is(catan_result, "CATaNResult")) {
        stop("catan_result must be a CATaNResult object")
    }
    if (!is(snp_gr, "GRanges")) {
        stop("snp_gr must be a GRanges object")
    }
    if (!dir.exists(peak_dir)) {
        stop("peak_dir does not exist: ", peak_dir)
    }

    tf_sample_loading <- as.data.frame(catan_result@tf_sample_loading)
    tf_names   <- rownames(tf_sample_loading)
    cc_names   <- colnames(tf_sample_loading)
    n_cc       <- length(cc_names)

    ## ---- Check that BED files exist ----
    bed_files <- file.path(peak_dir, paste0(tf_names, ".hg19.bed"))
    names(bed_files) <- tf_names
    exists_mask <- file.exists(bed_files)

    if (!any(exists_mask)) {
        stop("No TF peak BED files found in ", peak_dir,
             "\nExpected files like: ", basename(bed_files[1]))
    }
    if (!all(exists_mask)) {
        missing_n <- sum(!exists_mask)
        if (verbose) {
            message(sprintf(
                "Warning: %d/%d TF BED files not found, skipping",
                missing_n, length(tf_names)))
        }
        tf_names  <- tf_names[exists_mask]
        bed_files <- bed_files[exists_mask]
    }

    ## ---- Load all TF peaks once ----
    if (verbose) message("Loading TF peak BED files...")
    peak_list <- lapply(seq_along(tf_names), function(i) {
        tf <- tf_names[i]
        gr <- tryCatch(
            rtracklayer::import.bed(bed_files[i]),
            error = function(e) {
                ## Fallback: read as plain BED3
                bed <- utils::read.table(bed_files[i], header = FALSE,
                                         sep = "\t", stringsAsFactors = FALSE)
                GenomicRanges::GRanges(
                    seqnames = bed[, 1],
                    ranges   = IRanges::IRanges(start = bed[, 2] + 1L,
                                                end   = bed[, 3])
                )
            }
        )
        GenomicRanges::mcols(gr)$tf_name <- tf
        gr
    })
    names(peak_list) <- tf_names

    if (verbose) {
        total_peaks <- sum(vapply(peak_list, length, integer(1)))
        message(sprintf("  Loaded %d TFs, %d total peaks",
                        length(tf_names), total_peaks))
    }

    ## ---- Process each CC ----
    snp_scores_list <- list()

    for (cc_idx in seq_len(n_cc)) {
        cc_name <- cc_names[cc_idx]
        if (verbose) message(sprintf("Processing %s (%d/%d)...",
                                     cc_name, cc_idx, n_cc))

        ## Get scores for this CC
        scores <- tf_sample_loading[tf_names, cc_idx]
        names(scores) <- tf_names

        ## Process chromosomes in parallel
        chr_results <- BiocParallel::bplapply(chromosomes, function(chr) {
            .align_one_chromosome(
                chr       = chr,
                snp_gr    = snp_gr,
                peak_list = peak_list,
                scores    = scores,
                tf_names  = tf_names
            )
        }, BPPARAM = BPPARAM)

        ## Combine all chromosomes
        non_empty <- vapply(chr_results, length, integer(1)) > 0
        if (any(non_empty)) {
            all_snp_gr <- do.call(c, chr_results[non_empty])
            all_snp_gr <- sort(all_snp_gr)
        } else {
            all_snp_gr <- GenomicRanges::GRanges()
        }
        snp_scores_list[[cc_name]] <- all_snp_gr

        if (verbose) {
            n_scored <- sum(GenomicRanges::mcols(all_snp_gr)$score != 0)
            message(sprintf("  %d SNPs total, %d with non-zero score",
                            length(all_snp_gr), n_scored))
        }
    }

    ## ---- Build CATaNAnnotation ----
    new("CATaNAnnotation",
        snp_scores   = snp_scores_list,
        top_snps     = list(),
        bottom_snps  = list(),
        percentile   = 0.1,
        population   = population,
        catan_result = catan_result)
}


#' Process one chromosome: align CC scores to SNPs
#' @noRd
.align_one_chromosome <- function(chr, snp_gr, peak_list, scores, tf_names) {

    ## Subset SNPs to this chromosome
    chr_snps <- snp_gr[GenomicRanges::seqnames(snp_gr) == chr]
    if (length(chr_snps) == 0) {
        return(GenomicRanges::GRanges())
    }

    ## ---- Step 1: Collect all TF peaks on this chr with scores ----
    scored_peaks_list <- lapply(tf_names, function(tf) {
        tf_peaks <- peak_list[[tf]]
        chr_peaks <- tf_peaks[GenomicRanges::seqnames(tf_peaks) == chr]
        if (length(chr_peaks) == 0) return(NULL)
        GenomicRanges::mcols(chr_peaks)$score <- scores[tf]
        chr_peaks
    })
    scored_peaks_list <- scored_peaks_list[
        !vapply(scored_peaks_list, is.null, logical(1))
    ]

    if (length(scored_peaks_list) == 0) {
        ## No peaks on this chromosome -> all SNPs get score 0
        GenomicRanges::mcols(chr_snps)$score <- 0
        return(chr_snps)
    }

    all_peaks <- do.call(c, scored_peaks_list)
    all_peaks <- sort(all_peaks)

    ## ---- Step 2: Intersect peaks with SNPs ----
    ## equivalent to: bedtools intersect -a peaks -b snps
    hits <- GenomicRanges::findOverlaps(chr_snps, all_peaks)

    ## ---- Step 3: Sum scores per SNP ----
    ## equivalent to: bedtools merge -d -1 -c 4 -o sum
    all_scores <- rep(0, length(chr_snps))

    if (length(hits) > 0) {
        snp_indices  <- S4Vectors::queryHits(hits)
        peak_indices <- S4Vectors::subjectHits(hits)
        peak_scores  <- GenomicRanges::mcols(all_peaks)$score[peak_indices]

        ## Aggregate: sum of scores per SNP
        score_sums <- tapply(peak_scores, snp_indices, sum)
        scored_snp_idx <- as.integer(names(score_sums))
        all_scores[scored_snp_idx] <- as.numeric(score_sums)
    }

    ## ---- Step 4: Assign scores (zero-fill is implicit) ----
    GenomicRanges::mcols(chr_snps)$score <- all_scores
    chr_snps
}

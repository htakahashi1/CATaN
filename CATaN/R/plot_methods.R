#' Plot CC scatter plots from CATaNResult
#'
#' @param object A \code{\linkS4class{CATaNResult}} object.
#' @param cc_x Integer. CC number for x-axis (default 1).
#' @param cc_y Integer. CC number for y-axis (default 2).
#' @param space Character. Either \code{"tf"} or \code{"sample"} to plot
#'   TF weights or sample weights respectively (default \code{"tf"}).
#' @param color_by Character. Column name in sample_metadata for colouring
#'   points (only used when \code{space = "sample"}).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object (if ggplot2 is available) or
#'   base R plot (invisible NULL).
#'
#' @examples
#' set.seed(42)
#' tf <- matrix(sample(0:1, 500, replace = TRUE), 50, 10)
#' rownames(tf) <- paste0("g", seq_len(50))
#' colnames(tf) <- paste0("TF", seq_len(10))
#' ex <- matrix(rnorm(250, 5), 50, 5)
#' rownames(ex) <- rownames(tf)
#' colnames(ex) <- paste0("s", seq_len(5))
#' res <- run_cca(tf, ex, n_cc = 3L, n_hvg = 30L,
#'                verbose = FALSE)
#' plotCC(res, cc_x = 1, cc_y = 2, space = "tf")
#'
#' @export
setGeneric("plotCC",
    function(object, cc_x = 1L, cc_y = 2L, space = "tf", color_by = NULL, ...)
        standardGeneric("plotCC"))

#' @rdname plotCC
#' @export
setMethod("plotCC", "CATaNResult",
    function(object, cc_x = 1L, cc_y = 2L, space = "tf", color_by = NULL, ...) {
        space <- match.arg(space, c("tf", "sample"))

        if (space == "tf") {
            df <- as.data.frame(object@tf_sample_loading)
        } else {
            df <- as.data.frame(object@transcriptome_sample_loading)
        }

        x_col <- paste0("CC", cc_x, "_rotation")
        y_col <- paste0("CC", cc_y, "_rotation")

        if (!x_col %in% colnames(df) || !y_col %in% colnames(df)) {
            stop(sprintf("CC%d or CC%d not found in loading", cc_x, cc_y))
        }

        ## SCF for axis labels
        d <- object@svd_result$d
        scf_x <- round(d[cc_x]^2 / sum(d^2) * 100, 1)
        scf_y <- round(d[cc_y]^2 / sum(d^2) * 100, 1)

        xlab <- sprintf("CC%d (SCF=%.1f%%)", cc_x, scf_x)
        ylab <- sprintf("CC%d (SCF=%.1f%%)", cc_y, scf_y)

        if (requireNamespace("ggplot2", quietly = TRUE)) {
            p <- ggplot2::ggplot(df, ggplot2::aes(
                    x = .data[[x_col]], y = .data[[y_col]])) +
                ggplot2::geom_point(shape = 21, size = 3, alpha = 0.7) +
                ggplot2::labs(x = xlab, y = ylab) +
                ggplot2::theme_bw(base_size = 14)

            if (!is.null(color_by) && space == "sample") {
                meta <- as.data.frame(object@sample_metadata)
                if (color_by %in% colnames(meta)) {
                    df[[color_by]] <- meta[rownames(df), color_by]
                    p <- ggplot2::ggplot(df, ggplot2::aes(
                            x = .data[[x_col]], y = .data[[y_col]],
                            fill = .data[[color_by]])) +
                        ggplot2::geom_point(shape = 21, size = 3, alpha = 0.7) +
                        ggplot2::labs(x = xlab, y = ylab) +
                        ggplot2::theme_bw(base_size = 14)
                }
            }
            return(p)
        } else {
            plot(df[[x_col]], df[[y_col]],
                 xlab = xlab, ylab = ylab,
                 pch = 21, cex = 1.2, ...)
        }
    }
)

#' Functions to indent the plot to include the entire labels
#'
#' @param plt a ggplot object
#' @param xInd a scalar or a vector of length 2,
#' specifying the indentation left and right of the plot to allow for the labels
#'  to be printed entirely
#' @param yInd a a scalar or a vector of length 2,
#' specifying the indentation top and bottom of the plot
#'  to allow for the labels to be printed entirely
#' @return a ggplot object, squared
indentPlot <- function(plt, xInd = 0, yInd = 0) {
    return(plt +
            expand_limits(
            x = ggplot_build(plt)$layout$panel_params[[1]]$x.range +
        if (length(xInd) == 1) xInd * c(-1,
            1) else xInd,
        y = ggplot_build(plt)$layout$panel_params[[1]]$y.range +
        if (length(yInd) == 1) yInd * c(-1,
            1) else yInd))
}

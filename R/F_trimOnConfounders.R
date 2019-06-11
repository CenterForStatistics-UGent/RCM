#' Trim based on confounders to avoid taxa with only zero counts
#'
#' @param confounders a nxt confounder matrix
#' @param X the nxp data matrix
#' @param prevCutOff a scalar between 0 and 1, the prevalence cut off
#' @param minFraction a scalar between 0 and 1,
#'  each taxon's total abundance should equal at least the number of samples n
#'   times minFraction, otherwise it is trimmed
#' @param n the number of samples
#'
#' Should be called prior to fitting the independence model
#'
#' @return A trimmed data matrix nxp'
trimOnConfounders = function(confounders,
    X, prevCutOff, minFraction, n) {
    trimmingID = apply(X, 2, function(x) {
        # Over taxa Over confounding variables
        any(apply(confounders, 2, function(conf) {
            tapply(X = x, INDEX = conf, FUN = function(y) {
                mean(!(y %in% c(0, NA))) <= prevCutOff |
                sum(y, na.rm = TRUE) < (n * minFraction)
            })  #Any all-zero subgroup?
        }))
    })

    if (sum(!trimmingID) <= 1) {
        stop("All taxa would be trimmed,
        please provide a covariate with less levels,
        or reduce the prevalence cut-off! \n")
    }

    X[, !trimmingID]  #Return trimmed X
}

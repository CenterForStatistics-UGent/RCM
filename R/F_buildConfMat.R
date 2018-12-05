#' A function to build the confounder matrices
#'
#' @param x a matrix, data frame or character string
#' @param ... further arguments passed on to other methods
#'
#' For the preliminary trimming, we do not include an intercept,
#' but we do include all the levels of the factors using contrasts=FALSE:
#'  we want to do the trimming in every subgroup, so no hidden reference levels
#'   For the filtering we just use a model with an intercept and
#'    treatment coding, here the interest is only in adjusting the offset
#'
#' @return a list with components
#' \item{confModelMatTrim}{A confounder matrix without intercept, with all
#'  levels of factors present. This will be used to trim out taxa that have
#'   zero abundances in any subgroup defined by confounders}
#' \item{confModelMat}{A confounder matrix with intercept,
#' and with reference levels for factors absent.
#' This will be used to fit the model to modify the independence model,
#' and may include continuous variables}
buildConfMat = function(x, ...) {
    UseMethod("buildConfMat", x)
}

#' buildConfMat.data.frame
#'
#' @param confounders a data frame of confounders
#' @param n the number of rows of the count matrix
#'
#' @return see buidConfMat
buildConfMat.data.frame = function(confounders, n) {
    if (n != NROW(confounders)) {
        # Check dimensions
        stop("Data and confounder matrix do not have the same number
        of samples! \n")
    }
    if (anyNA(confounders)) {
        stop("Confounders contain missing values!\n")
    }
    # No intercept or continuous variables for preliminary
    # trimming
    confModelMatTrim = model.matrix(object = as.formula(paste("~",
        paste(names(confounders)[vapply(FUN.VALUE = TRUE, confounders,
            is.factor)], collapse = "+"), "-1")), data = confounders,
        contrasts.arg = lapply(confounders[vapply(FUN.VALUE = TRUE,
            confounders, is.factor)], contrasts, contrasts = FALSE))
    # With intercept for filtering
    confModelMat = model.matrix(object = as.formula(paste("~",
        paste(names(confounders), collapse = "+"))), data = confounders,
        contrasts.arg = lapply(confounders[vapply(FUN.VALUE = TRUE,
            confounders, is.factor)], contrasts, contrasts = TRUE))
    list(confModelMatTrim = confModelMatTrim, confModelMat = confModelMat)
}
#' buildConfMat.character
#' @param confounders a numeric matrix of confounders
#' @param physeq a physeq object with a sample_data slot
#'
#' @return see buidConfMat.numeric
buildConfMat.character = function(confounders, physeq) {
    if (!is(physeq, "phyloseq")) {
        stop("Providing confounders through variable names is only allowed
        if phyloseq object is provided! \n")
    }
    confounders = data.frame(get_variable(physeq, confounders))
    # The dataframe with the confounders
    buildConfMat(confounders, n = nsamples(physeq))
}
buildConfMat.default = function(...) {
    stop("Please provide the confounders either as
    dataframe, or character string! \n")
}

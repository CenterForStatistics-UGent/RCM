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
#' @importFrom stats model.matrix
#' @rdname buildConfMat
#' @export
setGeneric("buildConfMat", function(x, ...) standardGeneric("buildConfMat"))

#' buildConfMat.data.frame
#'
#' @param x a data frame of confounders
#' @param n the number of rows of the count matrix
#'
#' @return see buidConfMat
setMethod("buildConfMat", "data.frame",  function(x, n) {
    if (n != NROW(x)) {
        # Check dimensions
        stop("Data and confounder matrix do not have the same number
        of samples! \n")
    }
    if (anyNA(x)) {
        stop("Confounders contain missing values!\n")
    }
    #Check alias structure
    checkAlias(x, names(x))
    # No intercept or continuous variables for preliminary
    # trimming
    confModelMatTrim = model.matrix(object = as.formula(paste("~",
        paste(names(x)[vapply(FUN.VALUE = TRUE, x,
            is.factor)], collapse = "+"), "-1")), data = x,
        contrasts.arg = lapply(x[vapply(FUN.VALUE = TRUE,
            x, is.factor)], contrasts, contrasts = FALSE))
    # With intercept for filtering
    confModelMat = model.matrix(object = as.formula(paste("~",
        paste(names(x), collapse = "+"))), data = x,
        contrasts.arg = lapply(x[vapply(FUN.VALUE = TRUE,
            x, is.factor)], contrasts, contrasts = TRUE))
    list(confModelMatTrim = confModelMatTrim, confModelMat = confModelMat)
})
#' buildConfMat.character
#' @param x a numeric matrix of confounders
#' @param physeq a physeq object with a sample_data slot
#'
#' @return see buidConfMat.numeric
setMethod("buildConfMat", "character",  function(x, physeq) {
    if (!is(physeq, "phyloseq")) {
        stop("Providing confounders through variable names is only allowed
        if phyloseq object is provided! \n")
    }
    confounders = as(sample_data(physeq),"data.frame")[,make.names(x),
                                           drop = FALSE]
    # The dataframe with the confounders
    buildConfMat(confounders, n = nsamples(physeq))
})

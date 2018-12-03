#' A function to build the design matrix
#'
#' @param sampleScore a vector of environmental scores
#' @param responseFun A character string, indicating the shape
#' of the response function
#'
#' For dynamic response function estimation, the same desing matrix
#' as for the quadratic one is returned.
#' Will throw an error when an unknown repsonse function is provided
#'
#' @return A design matrix of dimension n-by-f
buildDesign = function(sampleScore, responseFun) {
    # With intercept
    design = switch(responseFun, linear = model.matrix(~sampleScore), 
        quadratic = model.matrix(~sampleScore + 
            I(sampleScore^2)), dynamic = model.matrix(~sampleScore + 
            I(sampleScore^2)), stop("Unknown response function"))
}

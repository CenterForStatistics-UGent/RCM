#' A function to construct a model matrix of a certain degree
#' @param y the variable
#' @param degree the degree
#'
#' @return A model matrix with degree+1 columns and as many rows as lenght(y)
getModelMat = function(y, degree) {
    model.matrix(formula(paste("~", paste(paste("I(y^", 
        seq_len(degree), ")"), collapse = "+"))))
}

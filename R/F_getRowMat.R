#' A function that returns a matrix of scores resulting from a response function and an environmental gradient.
#'
#' @param sampleScore a vector of length n with sample scores
#' @param responseFun a character string, the type of response function, either "linear" or "quadratic"
#' @param NB_params a v-by-p matrix of parameters of theresponse function
#'
#' Multiplying the old offset with the exponent matrix times the importance parameter obtains the new one based on lower dimension
#'
#' @return a n-by-p matrix of scores
getRowMat = function(sampleScore, responseFun, NB_params){
getModelMat(sampleScore, responseFun) %*% NB_params
}
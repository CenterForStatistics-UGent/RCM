#' A function that returns the model matrix, given sample score vector and response function.
#'
#' @param sampleScore a vector of length n with sample scores
#' @param responseFun a character string, the type of response function, either "linear" or "quadratic"
#'
#' @return the design matrix
getModelMat = function(sampleScore,responseFun){
  switch(responseFun,
       linear = model.matrix(~sampleScore),
       quadratic = model.matrix(~sampleScore + I(sampleScore^2) ),
       stop("Response function unknown! \n"))
}
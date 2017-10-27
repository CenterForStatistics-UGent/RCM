#' A function to perform correspondence analysis, regular or canonical, with the possibility of conditioning.
#'
#' @param X a n-by-p count matrix
#' @param Y an n-by-d covariate matrix
#' @param Z an n-by-z matrix of conditioning variables
#'
#' @return an svd() object in the constrained case, and an cca() object in the constrained one
#' The unconstrained approach is the approximation to the Pearson's chi-squared, as defined by Gower _et al._ in their book "Understanding Biplots". The constrained approach is the one implemented in vegan, based on a weighted regression of the departure matrix on the environmental scores. We cannot reproduce its results by the code given by ter Braak based on eigenvlaue or singular value decompositions.
caSVD = function(X, Y = NULL, Z = NULL){
  # @param X: the nxp count matrix
  # @param Y(optional): a nxd matrix of covariates, or a vector of variable names if X is a phyloseq object. If null regular CA is performed
  # @param Z(optional): a nxz matrix of variables to condition on, or a vector of variable names if X is a phyloseq object.

  # @return: the singular value decomposition of the matrix of pearson residuals
  if(class(X)=="phyloseq"){
    if(!is.null(Y)){Y = model.matrix(data = data.frame(sample_data(object = X))[,Y], ~.)}
    if(!is.null(Z)){Z = model.matrix(data = sample_data(object = X)[,Z], ~.)}
    X = if (taxa_are_rows(X)) t(otu_table(X)@.Data) else otu_table(X)@.Data
  }
  C = colSums(X)
  R = rowSums(X)
  E = outer(R,C)/sum(X) #Expected counts under independence
  if(is.null(Y) & is.null(Z)){ #Unconstrained analysis
    Goal = diag(1/sqrt(R)) %*% (X-E) %*% diag(1/sqrt(C))
    dimnames(Goal) = dimnames(X)
    svd(Goal)
  } else { #Constrained analysis
    vegan:::cca(X = X, Y=Y, Z=Z)
  }
}
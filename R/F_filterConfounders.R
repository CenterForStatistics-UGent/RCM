#'Filters out the effect of known confounders. This is done by fitting interactions of every taxon with the levels of the confounders. It returns a modified offset matrix for the remainder of the fitting procedure.
#'
#' @param muMarg a nxp matrix, the current offset
#' @param confounders a nxt confounder matrix
#' @param X the nxp data matrix
#' @param thetas a vector of length p with the current dispersion estimates
#' @param p an integer, the number of columns of X
#' @param n an integer, the number of rows of X
#' @param tol a scalar, the convergence tolerance
#' @param maxIt maximum number of iterations
#'
#' Fits the negative binomial mean parameters and overdispersion parameters iteratively. Convergence is determined based on the L2-norm of the absolute change of mean parameters
#'
#' @return a list with components:
#' \item{thetas}{new theta estimates}
#' \item{muMarg}{a new offset matrix}
#' \item{NB_params}{The estimated parameters of the interaction terms}

filterConfounders = function(muMarg, confMat, X, thetas, p, n, nleqslv.control, trended.dispersion, tol = 1e-3, maxIt = 20){
  NB_params = matrix(0, ncol(confMat), p)

  iter = 1
  while((iter ==1) || ((iter <= maxIt) && (!convergence))){

    NB_params_old = NB_params

    NB_params = sapply(seq_len(p), function(i){
      nleq = try(nleqslv(NB_params[,i] , reg = confMat,  fn = dNBllcol_constr, theta = thetas[i], muMarg = muMarg[,i], X = X[,i], control = nleqslv.control, jac = JacCol_constr, psi = 1)$x)  # , ... Fit the taxon-by taxon NB with given overdispersion parameters and return predictions
      if(class(nleq)=="try-error"| anyNA(nleq)| any(is.infinite(nleq))){
        nleq = nleqslv(NB_params[,i] , reg = confMat,  fn = dNBllcol_constr, theta = thetas[i], muMarg = muMarg[,i], X = X[,i], control = nleqslv.control, psi = 1)$x # If fails try with numeric jacobian
      }
      return(nleq)
    }) #Estimate response functions

    thetas = estDisp(X = X, cMat = matrix(0, ncol=p), rMat = matrix(0, nrow=n), psis = 0, muMarg = muMarg * exp(confMat %*% NB_params), trended.dispersion = trended.dispersion) #Estimate overdispersion

    iter = iter +1

    convergence = sqrt(mean((NB_params_old-NB_params)^2)) < tol #Check for convergence, L2-norm

  }

  list(muMarg = muMarg * exp(confMat %*% NB_params), thetas = thetas, NB_params = NB_params)
}
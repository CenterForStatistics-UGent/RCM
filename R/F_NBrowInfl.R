#' The influence function for the row scores
#'
#' @param rcm an rcm object
#' @param Dim the required dimension
#'
#' @return A list with components
#' \item{score}{a matrix with components of the score function}
#' \item{InvJac}{A square matrix of dimension n with the components of the Jacobian related to the row scores}
NBrowInfl = function(rcm, Dim){
  reg = rcm$psis[Dim] * rcm$cMat[Dim,]
  mu = extractE(rcm, seq_len(Dim)) #Take also lower dimensions into account here
  thetaMat = extractDisp(rcm, mu, Dim)
  lambdaRow = rcm$lambdaRow[seq_k(Dim)]
  rMatK = rcm$rMat[,seq_len(Dim-1),drop=FALSE]
  tmp = if(k>1) lambdaRow[-(1:2)] %*% cMatK else 0

  score= reg*(X-mu)/(1+mu/thetaMat) + c(rcm$rowWeights*(lambdaRow[1] + lambdaRow[2]*2*rcm$rMat[,Dim] + tmp))

  JacobianInv = solve(NBjacobianRow(beta = c(rcm$rMat[,Dim], lambdaRow), X = X, reg= reg, thetas = thetaMat, muMarg = muMarg, k = k, p = p, n=n, rowWeights = rowWeights , nLambda = nLambda, rMatK = rMatK)) #Inverse Jacobian

  #After a long thought: The X's do not affect the estimation of the lambda parameters!
  #Matrix becomes too large: return score and inverse jacobian
  return(list(score=score, InvJac = JacobianInv[1:n,1:n]))
}
#' A function that the components of the influence functions of the environmental gradient.
#'
#' @param rcm an rcm object
#' @param Dim the required dimension
#'
#' @return An n-by-p-by-d array with the influence of every observation on every alpha parameter
NBalphaInfl = function(rcm, Dim){
  if(length(Dim)>1) {stop("Influence of only one dimension at the time can be extratced! \n")}
  #Extract the parameters
  alpha = rcm$alpha[, Dim]
  centMat = buildCentMat(rcm)
  nLambda1s = NROW(centMat)
  lambdas = rcm$lambdasAlpha[seq_k(Dim, nLambda1s)]
  lambda1s = lambdas[seq_len(nLambda1s)] #Multiple centering requirements now!
  lambda2 = lambdas[nLambda1s+1]
  lambda3 = if(Dim==1) {0} else {lambdas[-seq_len(nLambda1s)]}
  responseFun = rcm$responseFun
  CC = rcm$covariates
  X = rcm$X
  p = ncol(X)
  n = nrow(X)
  d = ncol(CC)
  envGradEst = rcm$envGradEst
  thetaMat = matrix(rcm$thetas, byrow = TRUE, n, p)
  NB_params = rcm$NB_params[,,Dim]
  NB_params_noLab = rcm$NB_params_noLab[,Dim]
  psi = rcm$psis[Dim]

  sampleScore = CC %*% alpha #A linear combination of the environmental variables yields the sampleScore
  mu = extractE(rcm, seq_len(Dim))
  muMarg = extractE(rcm, seq_len(Dim-1))
  tmp = (X-mu)/(1+mu/thetaMat)
  tmp2 = rowMultiply(tmp, NB_params[2,])

  if(envGradEst == "LR"){
    mu0 = muMarg * c(exp(buildDesign(sampleScore, responseFun) %*% NB_params_noLab*psi))
    tmp0 = (X-mu0)/(1+mu0/thetaMat)
  }
  # score = array(0, dim = c(n,p,d + nLambda1s + Dim +1))
  score = switch(responseFun, #A n-by-p-by-d array
               "linear" = if(envGradEst == "LR"){psi * (vapply(seq_len(d), FUN.VALUE = tmp, function(i){tmp2*CC[,i]}) - NB_params_noLab[2]*vapply(seq_len(d), FUN.VALUE = tmp0, function(i){tmp0*CC[,i]}))
                 } else {
                   psi * (vapply(seq_len(d), FUN.VALUE = tmp, function(i){tmp2*CC[,i]}))
                   },
               "quadratic" = if(envGradEst == "LR"){psi * (
                 vapply(seq_len(d), FUN.VALUE = tmp, function(i){tmp2*CC[,i]}) +
                   2*vapply(seq_len(d), FUN.VALUE = tmp, function(i){rowMultiply(tmp, NB_params[3,])*CC[,i]*c(sampleScore)}) -
                   NB_params_noLab[2]*vapply(seq_len(d), FUN.VALUE = tmp0, function(i){tmp0*CC[,i]}) -
                   2*NB_params_noLab[3]*vapply(seq_len(d), FUN.VALUE = tmp, function(i){tmp0*CC[,i]*c(sampleScore)}))
                 } else {
                     psi * (vapply(seq_len(d), FUN.VALUE = tmp, function(i){tmp2*CC[,i]}) -
                            NB_params_noLab[2]*vapply(seq_len(d), FUN.VALUE = tmp0, function(i){tmp0*CC[,i]}))
                   },
               stop("Unknown response function provided! \n")) + #Restrictions do not depend on response function
    rep(c(lambda1s %*% centMat) + lambda2 * 2 * alpha + if(Dim>1) rowSums(rcm$alpha[, seq_len(Dim-1), drop = FALSE] %*% lambda3) else 0, each = n*p)

  JacobianInv = -solve(LR_nb_Jac(Alpha = c(alpha, lambda1s, lambda2, lambda3), X = X, CC = CC, responseFun = responseFun, psi = psi, NB_params = NB_params, NB_params_noLab = NB_params_noLab, d = d, k = Dim, centMat = centMat, nLambda = nLambda1s + Dim, nLambda1s = nLambda1s, thetaMat = thetaMat, muMarg = extractE(rcm, seq_len(Dim-1)), n = n, ncols = p, envGradEst = envGradEst, alphaK = rcm$alpha[, seq_len(Dim-1), drop = FALSE], preFabMat = 1+X/thetaMat))[seq_len(d), seq_len(d)] #Return only alpha indices, don't forget the minus sign!
rownames(JacobianInv) = colnames(JacobianInv) = names(alpha)
  tensor(score, JacobianInv, 3,1)
}
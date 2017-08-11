### Zero-inflated Poisson

```{r Zero-inflated poisson 1B1, purl=TRUE, echo = FALSE}
#----------------------------#
#expit
expit=function(x){
  tmp = exp(x)/(1+exp(x))
  tmp[is.na(tmp)]=1 #Adjust for overflow
  tmp}
#-----------------------------------#

## A function to perform the M step: maximize the likelihoods. This will again be an iterative process, estimating the parameters step by step. estimation of poisson and zero-inflated part can occur independently, which opens up opportunities for parallelization.

MstepZIP = function(Z, X, rMat, cMat, tMat, vMat,  k, muMarg,  zeroMarg, psis, chis, lambdaCol, lambdaRow, lambdaColZero, lambdaRowZero, nLambda, rowWeights, colWeights, rMatK, cMatK, vMatK, tMatK, twoCores=TRUE, tol=1e-3, psiTol = 1e-4, chiTol = psiTol, convNorm = 2, maxItMean=5, maxItZeroes= 20, n=n, p=p, global=global, nleqslv.control= nleqslv.control){

  #Optimization of the mean and zero-inflated components are independent (see Lambert 1992), so fork here
  resList = mclapply(mc.cores=twoCores + 1, c(meanEstZIP, ZIestZIP), function(fun){
    fun(X=X, rMat=rMat, cMat=cMat, tMat=tMat, chis=chis, vMat=vMat, zeroMarg = zeroMarg, lambdaCol=lambdaCol, lambdaRow=lambdaRow, lambdaRowZero=lambdaRowZero, lambdaColZero=lambdaColZero, psiTol=psiTol, chiTol=chiTol, tol=tol, convNorm = convNorm, nleqslv.control = nleqslv.control, global=global, nLambda=nLambda, k=k, Z=Z, muMarg=muMarg,n=n, p=p, psis=psis, maxItMean = maxItMean, maxItZeroes = maxItZeroes, rowWeights=rowWeights, colWeights=colWeights, rMatK = rMatK, cMatK = cMatK, tMatK = tMatK, vMatK = vMatK)
  })

  return(unlist(resList, recursive=FALSE))
}
#--------------------------------------#

# A function to estimate the mean component of the ZIP model by 1B1
meanEstZIP = function(X, rMat, cMat, Z, muMarg, k, global, nleqslv.control, tol, psiTol,lambdaCol, lambdaRow, convNorm,  nLambda, n, p, psis, rowWeights, colWeights, maxItMean = 10, maxItZeroes = 10, rMatK, cMatK,...){
  #Mean component

  iter = 1
  while((iter==1 || !converged) && iter<maxItMean){

    cat("Inner iteration(mean)", iter, "\n")

    psiOld = psis
    rMatOld = rMat
    cMatOld = cMat

    ## Row scores
    cat("Estimating row scores mean \n")
    regRow = cMat*psis
    rMatSol = try(nleqslv(fn = dZipMeanRmat, x = c(rMat, lambdaRow),X=X, reg =regRow, muMarg=muMarg, n=n, k=k, global=global, control = nleqslv.control, jac=ZipJacobianRmat, Z=Z, nLambda=nLambda, rowWeights=rowWeights, rMatK = rMatK)$x, silent=TRUE)

    if(class(rMatSol)!="try-error"){
      rMat = matrix(rMatSol[1:n], byrow=FALSE, ncol=1, nrow=n)
      lambdaRow = rMatSol[(n+1):length(rMatSol)]
    }

    #Normalize (speeds up algorithm if previous step had not converged)
    rMat =  rMat - sum(rMat * rowWeights)/sum(rowWeights)

    rMat = rMat/sqrt(sum(rowWeights * rMat^2))

    ## Column scores
    cat("Estimating column scores mean \n")
    regCol = rowMultiply(rMat,psis)
    cMatSol = try(nleqslv(fn = dZipMeanCmat, x = c(t(cMat), lambdaCol), X=X, reg=regCol, muMarg=muMarg, p=p, k=k, global=global, control = nleqslv.control, jac=ZipJacobianCmat, Z=Z, nLambda=nLambda, colWeights=colWeights, cMatK = cMatK)$x, silent=TRUE)
    if(class(cMatSol)!="try-error"){
      cMat = matrix(cMatSol[1:p], byrow=TRUE, nrow=1, ncol=p)
      lambdaCol = cMatSol[(p+1):length(cMatSol)]
    }

    cMat = cMat - sum(cMat * colWeights)/sum(colWeights)
    cMat = cMat/sqrt(sum(colWeights * cMat^2))

    ## Psis
    cat("\n Estimating psis (k =",k,") \n")

    regPsi =  rMat %*% cMat

    psisSol = try(sort(abs(nleqslv(fn = dZipMeanPsi, x = psis, X=X, reg=regPsi, Z=Z, muMarg=muMarg, global=global, control = nleqslv.control, jac=ZipJacobianPsi)$x), decreasing=TRUE), silent=TRUE)
    if(class(psisSol)!="try-error") psis=psisSol

    converged = all(abs(psiOld-psis) < psiTol) &&  (sum(abs(1-rMat/rMatOld)^convNorm))^(1/convNorm) < tol &&  (sum(abs(1-cMat/cMatOld)^convNorm))^(1/convNorm) < tol
    iter = iter +1
  }

  return(list(cMat=cMat, rMat=rMat, iterMean = iter, psis=psis, convergedMean=converged, lambdaCol = lambdaCol, lambdaRow=lambdaRow))
}
#--------------------------------------#

# A function to estimate the zero inflated component of the ZIP model

ZIestZIP = function(X, Z, muMarg, k, global, nleqslv.control, tol, chiTol, tMat, vMat, chis, zeroMarg, lambdaColZero, lambdaRowZero, convNorm, n, p, nLambda, psis, rowWeightsZeroNum, colWeightsZeroNum, rMatK, cMatK, tMatK, vMatK, maxItMean = 10, maxItZeroes = 10, ...){

  iter = 1
  while((iter==1 || !converged) && iter < maxItZeroes){
    chiOld = chis
    tMatOld = tMat
    vMatOld = vMat

    cat("Inner iteration(zeroes)", iter, "\n")

    ## Row scoers
    cat("Estimating row scores zeroes \n")
    regRowZero = vMat*chis
    tMatSol = try(nleqslv(fn = dZipMeanTmat, x = c(tMat, lambdaRowZero),  n=n,k=k, reg=regRowZero, global=global, control = nleqslv.control, zeroMarg = zeroMarg, jac=ZipJacobianTmat, Z=Z, nLambda=nLambda, rowWeights = rowWeightsZeroNum, tMatK = tMatK)$x, silent=TRUE)
    if(!inherits(tMatSol,"try-error")){
      tMat = matrix(tMatSol[1:n], byrow=FALSE, ncol=1, nrow=n)
      lambdaRowZero = tMatSol[(n+1):(n+nLambda)]
    }

    ## Column scores
    cat("Estimating column scores zeroes \n")
    regColZero = tMat*chis
    vMatSol = try(nleqslv(fn = dZipMeanVmat, x = c(t(vMat), lambdaColZero), reg=regColZero, p=p,k=k, global=global, control = nleqslv.control, zeroMarg = zeroMarg, jac=ZipJacobianVmat, Z=Z, nLambda=nLambda, colWeights=colWeightsZeroNum, vMatK = vMatK)$x, silent=TRUE)
    if(!inherits(vMatSol,"try-error")){
      vMat = matrix(vMatSol[1:p], byrow=TRUE, nrow=1, ncol=p)
      lambdaColZero = vMatSol[(p+1):(p+nLambda)]
    }

    # Chis
    cat("Estimating chis (zeroes) \n")
    regChis =  tMat %*% vMat

    chisSol = try(sort(abs(nleqslv(fn = dZipMeanChi, x = chis, reg=regChis, Z=Z, global=global, control = nleqslv.control, zeroMarg = zeroMarg, jac=ZipJacobianChi)$x), decreasing=TRUE), silent=TRUE)
    if(!inherits(chisSol,"try-error")){
      chis=chisSol
    }

    converged = all ((chiOld-chis) < chiTol) &&  (sum(abs(1-tMat/tMatOld)^convNorm))^(1/convNorm) < tol &&  (sum(abs(1-vMat/vMatOld)^convNorm))^(1/convNorm) < tol
    iter = iter +1
  }
  return(list(vMat=vMat, tMat=tMat, iterZI = iter, chis=chis, convergedZI=converged, lambdaColZero=lambdaColZero, lambdaRowZero=lambdaRowZero))
}

###Estimate the offsets
#--------------------------------------#
dZipMeanLibsizes = function(beta, X, Z, reg){
  # @param beta: a vector of logged library size estimates
  # @param y: the nxp data matrix
  # @param reg: the current logged abundance estimates
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(outer(beta, reg, "+"))
  rowSums((1 - Z)*(X - mu))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianLibsizes = function(beta, X, reg, Z){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(outer(beta, reg, "+"))
  diag(rowSums(mu*(Z-1)))
}
#--------------------------------------#
dZipMeanAbunds = function(beta, X, Z, reg){
  # @param beta: a vector of logged library size estimates
  # @param y: the nxp data matrix
  # @param reg: the current logged abundance estimates
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(outer(reg,beta, "+"))
  colSums((1 - Z)*(X - mu))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianAbunds = function(beta, X, reg, Z){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(outer(reg,beta, "+"))
  diag(colSums(mu*(Z-1)))
}

#All matrices X are considered to be nxp, i.e. samples are rows and taxa are columns

#--------------------------------------#
dZipMeanPsi = function(beta, X, muMarg, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(reg* beta) * muMarg
  sum((1 - Z)*(X - mu)*reg)
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianPsi = function(beta, X, reg, muMarg,Z){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(reg* beta) * muMarg
  sum(mu*reg^2*(Z-1))
}

#--------------------------------------#
dZipMeanRmat = function(beta, X, reg, muMarg, n, k, Z, nLambda, rowWeights, rMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg) * muMarg

  lambda1 = beta[n+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[n+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(n+3):length(beta)]}

  score = tcrossprod(reg ,(1-Z)*(X-mu)) + c(rowWeights*(lambda1 + lambda2*2*rMat + rMatK %*% lambda3))

  center = sum(rMat*rowWeights)
  unitSum = sum(rMat^2*rowWeights)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(rMatK, 2, function(x){
    sum(rMat*x*rowWeights)
  })
  return(c(score,center, unitSum, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianRmat = function(beta, X, reg, muMarg, n,k, nLambda, Z, rowWeights, rMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg) * muMarg

  Jac = matrix(0, nrow= n + nLambda, ncol=n + nLambda)
  #The symmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:n, n+1] = rowWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:n, n+2] = 2 *rMat*rowWeights

  if(k>1){
    Jac[1:n,(n+3):(n+nLambda)] = apply(rMatK, 2, function(x){rowWeights*x})
  }
  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:n,1:n]) = c(tcrossprod(-mu*(1-Z), reg^2) + 2*rowWeights*beta[n+2])
  Jac

}
#--------------------------------------#
dZipMeanCmat = function(beta, X, muMarg, p,k, Z, nLambda, colWeights, reg, cMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  cMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  mu = exp(reg %*% cMat) * muMarg

  lambda1 = beta[p+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[p+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(p+3):length(beta)]}

  score = crossprod(reg,(1-Z)*(X-mu)) + colWeights*(lambda1 + lambda2*2*cMat + lambda3 %*% cMatK)

  center = sum(colWeights*cMat)
  unitSum = sum(colWeights*cMat^2)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(cMatK,1,function(x){
    sum(cMat*x*colWeights)
  })
  return(c(score,center, unitSum, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianCmat = function(beta, X, psis, rMat, colWeights, k, p, muMarg, Z, nLambda, reg, cMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  cMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  mu = exp(reg %*% cMat) * muMarg

  Jac = matrix(0, nrow= p + nLambda, ncol = p + nLambda)
  #The suXmmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:p,(p+1)] = colWeights
  #Jac[1:(p*k),(p*k+1):((p+1)*k)] = sapply(1:k, function(K){c(rep(0,(K-1)*p),colWeights,rep(0,(k-K)*p))})
  Jac[1:p,p+2] = colWeights*2 *cMat

  #dLag²/ds_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:p,(p+3):(p+nLambda)] = apply(cMatK, 1,function(x){
      colWeights*x
    })
  }
  #Symmetrize
  Jac = Jac + t(Jac)

  diag(Jac[1:p,1:p]) = c(crossprod(-mu*(1-Z), reg^2)) + 2*beta[p+2]*colWeights
  Jac
}
#--------------------------------------#
#Estimate the offsets for the zeroes
dZipZeroCol = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  expitZero = expit(outer(reg, beta, "+"))
  colSums((Z-expitZero))
}

#--------------------------------------#
ZipJacobianZeroCol = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  expZero = exp(outer(reg, beta, "+"))
  tmp=expZero/(1+expZero)^2
  -diag(colSums(tmp))
}
#--------------------------------------#
#Estimate the offsets for the zeroes
dZipZeroRow = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  expitZero = expit(outer(beta,reg, "+"))
  rowSums((Z-expitZero))
}

#--------------------------------------#
ZipJacobianZeroRow = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  expZero = exp(outer(beta,reg, "+"))
  tmp=expZero/(1+expZero)^2
  -diag(rowSums(tmp))
}
#--------------------------------------#
dZipMeanChi = function(beta, Z, reg, zeroMarg, nLambda){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  GZero = expit(reg*beta + logit(zeroMarg))
  sum((Z-GZero)*reg)
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianChi = function(beta, Z, reg, zeroMarg, nLambda){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  expZero = exp(reg* beta + logit(zeroMarg))
  tmp=expZero/(1+expZero)^2
  -sum(reg^2*tmp)
}

#--------------------------------------#
dZipMeanTmat = function(beta, reg, k,n, Z, zeroMarg, nLambda, rowWeights, tMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  tMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  muZero = expit(tMat %*% reg+logit(zeroMarg))

  lambda1 = beta[n+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[n+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(n+3):length(beta)]}
  score =
    tcrossprod((Z-muZero), reg) + c(rowWeights*(lambda1 + lambda2*2*tMat + tMatK %*% lambda3))


  center = sum(tMat*rowWeights)
  unitSum = sum(tMat^2*rowWeights)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(tMatK, 2, function(x){
    sum(tMat*x*rowWeights)
  })
  return(c(score,center, unitSum, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianTmat = function(beta, reg, k, n, Z, zeroMarg, nLambda, rowWeights, tMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param k: a scalar, dimension of the RC solution

  tMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  muZero = exp(tMat %*% reg + logit(zeroMarg))

  Jac = matrix(0, nrow= n + nLambda, ncol= n + nLambda)
  #The suymmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:n, n+1] = rowWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:n, n+2] = 2*tMat*rowWeights
  tmp=muZero/(1+muZero)^2
  #dLag²/dr_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:n,(n+3):(n+nLambda)] = apply(tMatK, 2, function(x){rowWeights*x})
  }
  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:n,1:n]) = c(-tcrossprod(tmp, reg^2) + 2*beta[n+2]*rowWeights)

  Jac
}
#--------------------------------------#
dZipMeanVmat = function(beta, reg, k, p, Z, zeroMarg, nLambda, colWeights, vMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  vMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  muZero = exp(reg %*% vMat + logit(zeroMarg))

  lambda1 = beta[p+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[p+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(p+3):length(beta)]}

  score =
    crossprod(reg,(Z-muZero/(1+muZero))) + colWeights*(lambda1 + lambda2*2*vMat + (lambda3 %*% vMatK))

  center = sum(colWeights*vMat)
  unitSum = sum(colWeights*vMat^2)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(vMatK, 1,function(x){
    sum(x*vMat*colWeights)})

  return(c(score, center, unitSum, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZipJacobianVmat = function(beta, reg, k, p, Z, zeroMarg, nLambda, colWeights, vMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  vMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  muZero = exp(reg %*% vMat +logit(zeroMarg))

  Jac = matrix(0, nrow= p + nLambda, ncol=p + nLambda)
  #The suXmmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:p,(p+1)] = colWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:p,p+2] = colWeights*2 *vMat
  tmp=muZero/(1+muZero)^2
  if(k>1){
    Jac[1:p,(p+3):(p+nLambda)] = apply(vMatK, 1,function(x){
      colWeights*x
    })
  }
  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:p,1:p]) = c(t( -crossprod(reg^2,tmp)) + colWeights *2*beta[p+2])

  Jac
}
#Weighing by abs or relabunds really doesn't matter, only chis and psis get more inflated and deflated
logit=function(x){log(x/(1-x))}

RCM_ZIP = function(X, k, rowWeights, colWeights, weightsChar, tol = 1e-3, maxItOut = 500, psiTol = 1e-4, chiTol=psiTol, verbose = TRUE, global ="dbldog", nleqslv.control = list(), method="Broyden", twoCores = FALSE, convNorm = 2,  maxItMean = 20, maxItZeroes=30, ZIPRCM=NULL){

  # @param X: a nxp data matrix
  # @param k: a scalar, number of dimensions in the RC(M) model
  # @param tol(optional): a scalar, the relative convergende tolerance for the row scores and column scores parameters, defaults to 1e-3
  # @param Psitol(optional): a scalar, the relative convergence tolerance for the psi parameters, defaults to 1e-4
  # @param maxItOut(optional): an integer, the maximum number of iteration in the outer loop, defaults to 50
  # @param libSizes(optional) : a vector of length n with (known) library sizes. If not provided, rowSums of x are used
  # @param verbose(optional): a boolean, should information on iterations be printed? Defaults to TRUE
  # @param method(optional): Method for jacobian estimation , see nleqslv. Defaults to Broyden. The difference with the newton method is that the Jacobian is not recalculated at every iteration
  # @param global(optional): global strategy for solving non-linear systems , see nleqslv
  # @param nleqslv.control: a list with control options, see nleqslv
  # @param lambdaRow: a vector of length 2*k+k*(k-1)/2 with inital estimates or the lagrange multipliers for the row scores
  # @param lambdaCol: a vector of length 2*k+k*(k-1)/2 with inital estimates or the lagrange multipliers for the column scores
  # @param rMatInit(optional): a nxk matrix with initial row scores. If not provided values from the singular value decomposition will be used as starting values
  # @param cMatInit(optional): a pxk matrix with initial column scores. If not provided values from the singular value decomposition will be used as starting values
  # @param psisInit(optional): a vector of length k with inital values for the importance parameters psi. If not provided values from the singular value decomposition will be used as starting values
  # @param dispFreq: a scalar, how many iterations the algorithm should wait before reestimationg the dispersions
  # @param convNorm: a scalar, the norm to use to determine convergence

  # @return A list with elements:
  # @return psis: a vector of length k with estimates for the importance parameters psi
  # @return thetas: a vector of length p with estimates for the overdispersion
  # @return rMat: a nxk matrix with estimated row scores
  # @return cMat: a pxk matrix with estimated column scores
  # @return converged: a boolean indicating if the algorithm converged
  # @return rowRec: a n x k x maxItOut array with a record of all rMat estimates through the iterations
  # @return colRec: a k x p x maxItOut array with a record of all cMat estimates through the iterations
  # @return psiRec.: a k x maxItOut array with a record of all psi estimates through the iterations

  abunds = colSums(X)/sum(X)
  libSizes = rowSums(X)
  n=NROW(X)
  p=NCOL(X)

  logLibSizesMLE = log(libSizes)
  logAbundsMLE = log(abunds)
  # logitZeroRows = logit(rep(0.25,n)) #Start with a 25 % chance on a structural zero
  logitZeroCols = logit(rep(0.25,p))

  initIter = 1

  while((initIter ==1) || ((initIter <= maxItOut) && (!convergenceInit))){
    cat("Starting iteration ",initIter, " \n")

    libsOld = logLibSizesMLE
    absOld = logAbundsMLE
    # zeroRowsOld = logitZeroRows
    zeroColsOld = logitZeroCols
    Z = matrix(0,n,p)

    #E-step
    expMu = exp(outer(logLibSizesMLE, logAbundsMLE, "+"))
    regZero = matrix(logitZeroCols, n, p, byrow=TRUE)
    Z[X==0] = (1+exp(-regZero-expMu))[X==0]^(-1)

    #M-step

    initIterMean = 1
    while((initIterMean ==1) || ((initIterMean <= maxItMean) && (!convergenceInitMean))){
      cat("Mean iteration ",initIterMean, "\n")
      libsOldIn = logLibSizesMLE
      absOldIn = logAbundsMLE
      libsTmp = try(nleqslv(fn = dZipMeanLibsizes, x = logLibSizesMLE, X = X, reg=logAbundsMLE, global=global, control = nleqslv.control, jac=ZipJacobianLibsizes, method=method, Z=Z)$x, silent=TRUE)
      if(class(libsTmp)!="try-error"){ logLibSizesMLE = libsTmp}
      absTmp = try(nleqslv(fn = dZipMeanAbunds, x = logAbundsMLE, X = X, reg=logLibSizesMLE, global=global, control = nleqslv.control, jac=ZipJacobianAbunds, method=method, Z=Z)$x, silent=TRUE)
      if(class(absTmp)!="try-error"){ logAbundsMLE = absTmp}

      initIterMean = initIterMean + 1

      convergenceInitMean = ((initIterMean <= maxItMean) &&
                               ((sum(abs(1-logLibSizesMLE/libsOldIn)^convNorm)/n)^(1/convNorm) < tol) &&
                               ((sum(abs(1-logAbundsMLE/absOldIn)^convNorm)/p)^(1/convNorm) < tol))
    }

    #     initIterZero = 1
    #   while((initIterZero ==1) || ((initIterZero <= maxItZeroes) && (!convergenceInitZeroes))){
    #         cat("Zero iteration ",initIterZero, "\n")
    # zeroRowsOldIn = logitZeroRows
    # zeroColsOldIn = logitZeroCols

    #   zeroRowsTmp = try(nleqslv(fn = dZipZeroRow, x = logitZeroRows,   reg=logitZeroCols, global=global, control = nleqslv.control, jac=ZipJacobianZeroRow, method=method, Z=Z)$x, silent=TRUE)
    #  if(class(zeroRowsTmp)!="try-error"){ logitZeroRows = zeroRowsTmp} # We don't let zero-inflation depend on the samples, there is no baseline absence rate for a sample

    zeroColsTmp = try(nleqslv(fn = dZipZeroCol, x = logitZeroCols, reg=rep.int(0L,n), global=global, control = nleqslv.control, jac=ZipJacobianZeroCol, method=method, Z=Z)$x, silent=TRUE)
    if(class(zeroColsTmp)!="try-error"){ logitZeroCols = zeroColsTmp}

    # initIterZero = initIterZero + 1

    #    convergenceInitZeroes = ((initIterZero <= maxItZeroes) &&
    #                    # ((sum(abs(1-logitZeroRows/zeroRowsOldIn)^convNorm)/p)^(1/convNorm) < tol) &&
    #                    ((sum(abs(1-logitZeroCols/zeroColsOldIn)^convNorm)/p)^(1/convNorm) < tol))
    #   }

    initIter = initIter + 1

    convergenceInit = ((initIter <= maxItOut) &&
                         ((sum(abs(1-logLibSizesMLE/libsOld)^convNorm)/n)^(1/convNorm) < tol) &&
                         ((sum(abs(1-logAbundsMLE/absOld)^convNorm)/p)^(1/convNorm) < tol) &&
                         # ((sum(abs(1-logitZeroRows/zeroRowsOld)^convNorm)/p)^(1/convNorm) < tol) &&
                         ((sum(abs(1-logitZeroCols/zeroColsOld)^convNorm)/p)^(1/convNorm) < tol))
  }
  muMarg = exp(outer(logLibSizesMLE, logAbundsMLE, "+"))
  zeroMarg = expit(matrix(logitZeroCols, n, p, byrow=TRUE))#The marginals to be used as expectation under independence. These are augmented with the previously estimated dimensions every time

  rowRec = rowRecZeroes = array(0,dim=c(n,k, maxItOut))
  colRec = colRecZeroes = thetaRec = array(0,dim=c(k,p, maxItOut))
  psiRec = matrix(0, nrow=k,ncol=maxItOut)
  convergence = rep(FALSE, k)
  iterOut = rep(1,k)

  #If previous fit provided with higher or equal dimension, stop here
  if((!is.null(ZIPRCM)) ){
    if(ZIPRCM$fit != "RCM_ZIP"){
      stop("Fit provided is not of same type as the one requested! \n")
    } else if((k <= ZIPRCM$k)) {
      # stop("Fit provided is already of the required dimension or higher! \n")
    } else{
      for(i in c("rMat","cMat","psis","lambdaCol","lambdaRow", "lambdaRowZero","lambdaColZero","tMat","vMat","chis","Z","zeroMarg")){
        assign(i, ZIPRCM[[i]])
      }
    }
    #Otherwise try to use intelligent starting values
  } else{
    #Depending on the weighting schemes, use other starting values
    svdX = svd(diag(1/sqrt(libSizes)) %*% (X-muMarg)*(1-zeroMarg) %*% diag(1/sqrt(colSums(X))))#switch(weightsChar,
    #                 "marginalmarginal" = svd(diag(1/libSizes) %*% (X-muMarg) %*% diag(1/colSums(X))),
    #                 "marginaluniform" = svd(diag(1/libSizes) %*% (X-muMarg)),
    #                 "uniformmarginal" = svd((X-muMarg) %*% diag(1/colSums(X))),
    #                 "uniformuniform" = svd(X-muMarg))
    rMat = svdX$u[,1:k,drop=FALSE]
    cMat = t(svdX$v[,1:k,drop=FALSE])
    psis = log(svdX$d[1:k])

    #   #Redistribute some weight to fit the constraints
    #   psis = c(psis *t(apply(cMat, 1, function(colS){
    #       sqrt(sum(colWeights * colS^2))
    #   })) * apply(rMat, 2, function(rowS){
    #       sqrt(sum(rowWeights * rowS^2))
    #   }))
    #
    # #Normalize
    # cMat = t(apply(cMat, 1, function(colS){
    #       colS/sqrt(sum(colWeights * colS^2))
    #   }))
    # rMat = apply(rMat, 2, function(rowS){
    #       rowS/sqrt(sum(rowWeights * rowS^2))
    #   })

    #Initial estimates for zeroes is also based on an svd

    Xzeroes = X==0

    svdZero = svd(diag(sqrt(1/sqrt(libSizes))) %*%(Xzeroes-zeroMarg)%*% diag(1/sqrt(colSums(X))))

    tMat = svdZero$u[,1:k, drop=FALSE]
    vMat = t(svdZero$v[,1:k, drop=FALSE])
    chis = log(svdZero$d[1:k])

    #Redistribute some weight to fit the constraints
    # chis = c(chis *t(apply(vMat, 1, function(colS){
    #       sqrt(sum(colWeights * colS^2))
    #   })) * apply(tMat, 2, function(rowS){
    #       sqrt(sum(rowWeights * rowS^2))
    #   }))
    #
    # #Normalize
    # vMat = t(apply(vMat, 1, function(colS){
    #       colS/sqrt(sum(colWeights * colS^2))
    #   }))
    # tMat = apply(tMat, 2, function(rowS){
    #       rowS/sqrt(sum(rowWeights * rowS^2))
    #   })

    lambdaRow = lambdaCol = lambdaColZero=lambdaRowZero = rep.int(0,2*k+k*(k-1)/2)
  }

  rowRec = rowRecZero = array(0,dim=c(NROW(X),k, maxItOut))
  colRec = colRecZero = array(0,dim=c(k,NCOL(X), maxItOut))
  psiRec = chiRec = matrix(0,ncol=maxItOut, nrow=k)

  if(!is.null(ZIPRCM)){ #If fit provided, replace lower dimension starting values
    Kprev = ZIPRCM$k
    rMat[,1:Kprev] = ZIPRCM$rMat
    rowRec[,1:Kprev,] = ZIPRCM$rowRec
    cMat[1:Kprev,] = ZIPRCM$cMat
    colRec[1:Kprev,,] = ZIPRCM$colRec
    psis[1:Kprev] = ZIPRCM$psis
    psiRec[1:Kprev,] = ZIPRCM$psiRec
    lambdaCol[1:(Kprev*(2+(Kprev-1)/2))] = ZIPRCM$lambdaCol
    lambdaRow[1:(Kprev*(2+(Kprev-1)/2))] = ZIPRCM$lambdaRow
    tMat[,1:Kprev] = ZIPRCM$tMat
    rowRecZero[,1:Kprev,] = ZIPRCM$rowRecZero
    vMat[1:Kprev,] = ZIPRCM$vMat
    colRecZero[1:Kprev,,] = ZIPRCM$colRecZero
    chis[1:Kprev] = ZIPRCM$chis
    chiRec[1:Kprev,] = ZIPRCM$chiRec
    lambdaColZero[1:(Kprev*(2+(Kprev-1)/2))] = ZIPRCM$lambdaColZero
    lambdaRowZero[1:(Kprev*(2+(Kprev-1)/2))] = ZIPRCM$lambdaRowZero
    convergence[1:Kprev] = ZIPRCM$converged
    iterOut[1:Kprev] = ZIPRCM$iter
    zeroMarg = ZIPRCM$zeroMarg
  }

  minK = ifelse(is.null(ZIPRCM),1,Kprev+1)
  for (KK in minK:k){

    cat("Dimension" ,KK, "is being esimated \n")

    #Modify offsets if needed
    if(KK>1){
      muMarg = muMarg * exp(rMat[,(KK-1), drop=FALSE] %*% (cMat[(KK-1),, drop=FALSE]*psis[(KK-1)]))
      zeroMarg = expit(logit(zeroMarg) + tMat[,(KK-1), drop=FALSE] %*% (vMat[(KK-1),, drop=FALSE]*chis[(KK-1)]))
    }
    #A lambda parameter
    nLambda = KK + 1

    #The location of the lambda parameters
    idK = seq_k(KK)
    Z = matrix(0, n,p)

    ## 2) Propagation

    while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) && (!convergence[KK])))
    {

      if(verbose && iterOut[KK]%%1 == 0){
        cat("\n","Outer Iteration", iterOut[KK], "\n","\n")
        if(iterOut[KK]!=1){
          cat("Old psi-estimate: ", psiOld, "\n")
          cat("New psi-estimate: ", psis[KK], "\n")
          cat("Old chi-estimates: ", chiOld, "\n")
          cat("New chi-estimates: ", chis[KK], "\n")
        }
      }
      ## 2)a. Store old parameters
      psiOld = psis[KK]
      rMatOld = rMat[,KK]
      cMatOld = cMat[KK,]

      chiOld = chis[KK]
      tMatOld = tMat[,KK]
      vMatOld = vMat[KK,]

      #Expectation
      expMu = muMarg* exp(outer(rMat[,KK], cMat[KK,]*psis[KK]))
      regZero = logit(zeroMarg) + outer(tMat[,KK], vMat[KK,]*chis[KK])
      Z[X==0] = (1+exp(-regZero-expMu))[X==0]^(-1)

      #Maximization
      Mlist = MstepZIP(Z = Z, X = X, rMat = rMat[,KK, drop=FALSE], cMat = cMat[KK,, drop=FALSE], tMat = tMat[,KK, drop=FALSE], vMat = vMat[KK,, drop=FALSE], k = KK, n=n, p=p, zeroMarg = zeroMarg, psis = psis[KK],chis = chis[KK], twoCores=twoCores, tol = tol, psiTol = psiTol, chiTol = chiTol, convNorm = convNorm, global=global, nLambda = nLambda, nleqslv.control = nleqslv.control, lambdaCol = lambdaCol[idK], lambdaRow=lambdaRow[idK], lambdaColZero = lambdaColZero[idK], lambdaRowZero = lambdaRowZero[idK], maxItMean = maxItMean, maxItZeroes = maxItZeroes, muMarg=muMarg, colWeights = colWeights, rowWeights = rowWeights, rMatK = rMat[,1:(KK-1), drop=FALSE], cMatK = cMat[1:(KK-1),, drop=FALSE], tMatK = tMat[,1:(KK-1), drop=FALSE], vMatK = vMat[1:(KK-1),, drop=FALSE])
      #
      #   cat("Mlist:", str(Mlist), "\n")

      #Assign outcomes to tracking vectors and to this environment
      lambdaCol[idK] = Mlist$lambdaCol
      lambdaRow[idK] = Mlist$lambdaRow
      lambdaColZero[idK] = Mlist$lambdaColZero
      lambdaRowZero[idK] = Mlist$lambdaRowZero

      rowRec[,KK, iterOut[KK]] = rMat[,KK] = Mlist$rMat
      colRec[KK,, iterOut[KK]] = cMat[KK,] = Mlist$cMat
      rowRecZero[,KK, iterOut[KK]] = tMat[,KK] = Mlist$tMat
      colRecZero[KK,, iterOut[KK]] = vMat[KK,] = Mlist$vMat
      psiRec[KK, iterOut[KK]] = psis[KK] = Mlist$psis
      chiRec[KK, iterOut[KK]] = chis[KK] = Mlist$chis

      ## 2)f. Change iterator
      iterOut[KK] = iterOut[KK] + 1

      ##Check convergence  (any numbered norm for row and column scores)
      convergence[KK] = (iterOut[KK] <= maxItOut) &&
        (all(abs(1-psis/psiOld) < psiTol)) &&
        ((sum((1-rMatOld/rMat)^convNorm)/n)^(1/convNorm) < tol) &&
        ((sum((1-cMatOld/cMat)^convNorm)/p)^(1/convNorm) < tol)  &&
        (all(abs(1-chis/chiOld) < chiTol)) &&
        (sum(abs(1-tMat/tMatOld)^convNorm)/n)^(1/convNorm) < tol &&
        (sum(abs(1-vMat/vMatOld)^convNorm)/p)^(1/convNorm) < tol
    } # END while-loop until convergence
  } # END for-loop over dimensions

  ## 3) Termination
  rownames(rMat) = rownames(X)
  colnames(cMat) = colnames(X)
  rownames(cMat) = colnames(rMat) = paste0("Dim",1:k)

  if(!convergence[KK] ){
    warning("Algorithm did not converge! Check for errors or consider changing tolerances or number of iterations")
  }
  return(list(converged = convergence,rMat=rMat, cMat=cMat, psis = psis, X=X,
              rowRec = rowRec, colRec = colRec, psiRec = psiRec, lambdaRow = lambdaRow, lambdaCol = lambdaCol, lambdaRowZero = lambdaRowZero, lambdaColZero = lambdaColZero, chis = chis, tMat = tMat, vMat = vMat, zeroMarg = zeroMarg, chiRec = chiRec, rowRecZero = rowRecZero, colRecZero = colRecZero, iter=iterOut-1, Z=Z, fit="RCM_ZIP", libSizesMLE = exp(logLibSizesMLE), abundsMLE = exp(logAbundsMLE), taxaZeroes = expit(logitZeroCols)))
}
```

### Zero-inflated negative binomial

```{r ZINB, purl=TRUE, echo=FALSE}
###Estimate the offsets
#--------------------------------------#
dZinbMeanLibsizes = function(beta, X, Z, reg, thetas){
  # @param beta: a vector of logged library size estimates
  # @param y: the nxp data matrix
  # @param reg: the current logged abundance estimates
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(outer(beta, reg, "+"))
  rowSums((1-Z)*(X-mu)/(1+t(t(mu)/thetas)))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianLibsizes = function(beta, X, reg, Z, thetas){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(outer(beta, reg, "+"))
  diag(rowSums(((1+t(t(X)/thetas))*mu/(1+t(t(mu)/thetas))^2)*(Z-1)))
}
#--------------------------------------#
dZinbMeanAbunds = function(beta, X, Z, reg, thetas){
  # @param beta: a vector of logged library size estimates
  # @param y: the nxp data matrix
  # @param reg: the current logged abundance estimates
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(outer(reg,beta, "+"))
  colSums((1-Z)*(X-mu)/(1+t(t(mu)/thetas)))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianAbunds = function(beta, X, reg, Z, thetas){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(outer(reg,beta, "+"))
  diag(colSums(((1+t(t(X)/thetas))*mu/(1+t(t(mu)/thetas))^2)*(Z-1)))
}

#--------------------------------------#
#Estimate the offsets for the zeroes
dZinbZeroCol = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  expitZero = expit(outer(reg, beta, "+"))
  colSums((Z-expitZero))
}

#--------------------------------------#
ZinbJacobianZeroCol = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  expZero = exp(outer(reg, beta, "+"))
  tmp=expZero/(1+expZero)^2
  -diag(colSums(tmp))
}

#--------------------------------------#
#Estimate the offsets for the zeroes
dZinbZeroRow = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  expitZero = expit(outer(beta,reg, "+"))
  rowSums((Z-expitZero))
}

#--------------------------------------#
ZinbJacobianZeroRow = function(beta, Z, reg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  expZero = exp(outer(beta,reg, "+"))
  tmp=expZero/(1+expZero)^2
  -diag(rowSums(tmp))
}

#-----------------------------------#
## A function to perform the E-step

EstepNB  = function(X, rMat, cMat, tMat, vMat, muMarg, zeroMarg, psis, chis, thetas){

  # @return: The values of Z
  expMu = exp(rMat %*% (cMat*psis)) * muMarg
  pZero =  expit(tMat %*% (vMat * chis) + logit(zeroMarg))

  Z = X
  Z[X>0] = 0
  thetaMat=matrix(thetas, ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
  d0=dnbinom(0,mu=expMu, size=thetaMat)
  Z[X==0] = (pZero/((1-pZero)*d0 + pZero))[X==0]
  Z
}
#-----------------------------------#

## A function to perform the M step: maximize the likelihoods. This will again be an iterative process, estimating the parameters step by step. estimation of poisson and zero-inflated part can occur independently, which opens up opportunities for parallelization.

MstepNB = function(Z, X, rMat, cMat, tMat, vMat, muMarg, k,  zeroMarg, psis, chis, lambdaCol, lambdaRow, lambdaColZero, lambdaRowZero, twoCores=TRUE, tol=1e-3, psiTol = 1e-4, chiTol = psiTol, convNorm = 2 , maxItMean=20 , maxItZeroes= 20,n, p, global=global, nleqslv.control= nleqslv.control, nLambda, thetas, dispFreq,rowWeights, colWeights, rMatK, cMatK, tMatK, vMatK){

  #Optimization of the mena and zero-inflated components are independent (see Lambert 1992), so fork here
  resList = mclapply(mc.cores= 1+twoCores, c(meanEstZINB, ZIestNB), function(fun){
    fun(X=X, rMat=rMat, cMat=cMat, tMat=tMat, chis=chis, vMat=vMat, zeroMarg = zeroMarg, lambdaCol=lambdaCol, lambdaRow=lambdaRow, lambdaRowZero=lambdaRowZero, lambdaColZero=lambdaColZero, psiTol=psiTol, chiTol=chiTol, tol=tol, convNorm = convNorm, nleqslv.control = nleqslv.control, global=global, nLambda=nLambda, k=k, Z=Z, muMarg=muMarg,n=n, p=p, psis=psis, maxItMean = maxItMean, maxItZeroes = maxItZeroes, thetas=thetas, dispFreq=dispFreq, rowWeights = rowWeights, colWeights=colWeights, rMatK = rMatK, cMatK = cMatK, tMatK = tMatK, vMatK = vMatK)
  })

  return(unlist(resList, recursive=FALSE))
}
#--------------------------------------#

# A function to estimate the mean component of the ZIP model
meanEstZINB = function(X, rMat, cMat, Z, muMarg,  k, global, nleqslv.control, tol, psiTol, thetas, lambdaCol, lambdaRow, convNorm, dispFreq, nLambda, n, p, psis, maxItMean = 10, maxItZeroes = 10, rowWeights, colWeights, rMatK, cMatK,...){
  #Mean component

  iter = 1
  while((iter==1 || !converged) && iter<=maxItMean){

    cat("Inner iteration(mean)", iter, "\n")

    psiOld = psis
    rMatOld = rMat
    cMatOld = cMat

    cat("Estimating overdispersions \n")
    if(iter==1 | iter %% dispFreq ==0){ #Again too slow and unnecessary to reestimate overdispersions every time
      thetasTry = try(estDisp(X=X, cMat=cMat, rMat=rMat, muMarg=muMarg, psis=psis, dispWeights=t(1-Z)), silent=TRUE)
      if(class(thetasTry)!="try-error") thetas=thetasTry
    }

    cat("Estimating psis ( k =",k,") \n")
    regPsis = rMat %*% cMat
    psisSol = try(abs(nleqslv(fn = dZinbMeanPsi, reg=regPsis, x = psis, X=X, Z=Z, muMarg=muMarg, global=global, control = nleqslv.control, jac=ZinbJacobianPsi, thetas=thetas)$x), silent=TRUE)
    if(class(psisSol)!="try-error") psis=psisSol

    cat("Estimating row scores mean \n")
    regRows = cMat*psis
    rMatSol = try(nleqslv(fn = dZinbMeanRmat, x = c(rMat, lambdaRow), X=X, reg =regRows, muMarg=muMarg, k=k, n=n, global=global, control = nleqslv.control, jac=ZinbJacobianRmat, Z=Z, nLambda=nLambda, thetas=thetas, rowWeights=rowWeights, rMatK = rMatK)$x, silent=TRUE)
    if(class(rMatSol)!="try-error"){
      rMat = matrix(rMatSol[1:n], byrow=FALSE, ncol=1, nrow=n)
      lambdaRow = rMatSol[(n+1):(n+nLambda)]
    }

    ## Column scores
    cat("Estimating column scores mean \n")
    regCols = rMat*psis
    cMatSol = try(nleqslv(fn = dZinbMeanCmat, x = c(cMat, lambdaCol), reg = regCols, X = X, muMarg = muMarg, k = k, p = p, global = global, control = nleqslv.control, jac = ZinbJacobianCmat, Z = Z, nLambda = nLambda, thetas = thetas, colWeights=colWeights, cMatK = cMatK)$x, silent=TRUE)
    if(class(cMatSol)!="try-error"){
      cMat = matrix(cMatSol[1:p], byrow=TRUE, nrow=1, ncol=p)
      lambdaCol = cMatSol[(p+1):(p+nLambda)]
    }

    converged = all(abs(psiOld-psis) < psiTol) &&  (sum(abs(1-rMat/rMatOld)^convNorm))^(1/convNorm) < tol &&  (sum(abs(1-cMat/cMatOld)^convNorm))^(1/convNorm) < tol
    iter = iter +1
  }

  return(list(cMat=cMat, rMat=rMat, iterMean = iter, psis=psis, convergedMean=converged, lambdaCol = lambdaCol, lambdaRow=lambdaRow))
}
#--------------------------------------#

# A function to estimate the zero inflated component of the ZIP model
ZIestNB = function(X, Z, k, global, nleqslv.control, tol,  chiTol, tMat, vMat, chis, zeroMarg, lambdaColZero, lambdaRowZero, convNorm, n, p,  nLambda, rowWeights, colWeights, thetas, tMatK, vMatK, maxItMean = 10, maxItZeroes = 10,  ...){

  iter = 1
  while((iter==1 || !converged) && iter<=maxItZeroes){
    chiOld = chis
    tMatOld = tMat
    vMatOld = vMat

    cat("Inner iteration(zeroes)", iter, "\n")

    # Zero component
    ## Chis
    cat("Estimating chis (zeroes) \n")
    regChis = tMat %*% vMat
    chisSol = try(abs(nleqslv(fn = dZinbZeroChi, reg=regChis, x = chis, Z=Z, global=global, control = nleqslv.control, zeroMarg=zeroMarg, jac=ZinbJacobianChi)$x), silent=TRUE)
    if(!inherits(chisSol,"try-error")){
      chis=chisSol
    }

    ## Row scoers
    cat("Estimating row scores zeroes \n")
    regRows = vMat*chis
    tMatSol = try(nleqslv(fn = dZinbZeroTmat, x = c(tMat, lambdaRowZero), k=k, n=n, global=global, control = nleqslv.control, zeroMarg=zeroMarg, reg=regRows, jac=ZinbJacobianTmat, Z=Z, nLambda=nLambda, rowWeights=rowWeights, tMatK= tMatK)$x, silent=TRUE)
    if(!inherits(tMatSol,"try-error")){
      tMat = matrix(tMatSol[1:n], byrow=FALSE, ncol=1, nrow=n)
      lambdaRowZero = tMatSol[(n+1):(nLambda+n)]
    }

    ## Column scores
    cat("Estimating column scores zeroes \n")
    regCols = tMat*chis
    vMatSol = try(nleqslv(fn = dZinbZeroVmat, x = c(vMat, lambdaColZero), reg=regCols, k=k, p=p, global=global, control = nleqslv.control, zeroMarg=zeroMarg, jac=ZinbJacobianVmat, Z=Z, nLambda=nLambda, colWeights=colWeights, vMatK = vMatK)$x, silent=TRUE)
    if(!inherits(vMatSol,"try-error")){
      vMat = matrix(vMatSol[1:p], byrow=TRUE, nrow=1, ncol=p)
      lambdaColZero = vMatSol[(p+1):(p+nLambda)]
    }

    converged = all ((chiOld-chis) < chiTol) &&  (sum(abs(1-tMat/tMatOld)^convNorm))^(1/convNorm) < tol &&  (sum(abs(1-vMat/vMatOld)^convNorm))^(1/convNorm) < tol
    iter = iter +1
  }
  return(list(vMat=vMat, tMat=tMat, iterZI = iter, chis=chis, convergedZI=converged, lambdaColZero=lambdaColZero, lambdaRowZero=lambdaRowZero))
}
#--------------------------------------#

dZinbMeanPsi = function(beta, X, muMarg, Z, reg, thetas){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  mu = exp(reg* beta) * muMarg
  sum(reg*(1-Z)*((X-mu)/(1+t(t(mu)/thetas))))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianPsi = function(beta, X, reg, muMarg, Z, thetas){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  mu = exp(reg*beta) * muMarg
  -sum(reg^2*(1-Z)*((1+t(t(X)/thetas))*mu/(1+t(t(mu)/thetas))^2))
}

#--------------------------------------#
dZinbMeanRmat = function(beta, X, reg, muMarg, k, n, Z, nLambda, thetas, rowWeights, rMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg) * muMarg

  lambda1 = beta[n+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[n+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(n+3):length(beta)]}

  score = tcrossprod((1-Z)*((X-mu)/(1+t(t(mu)/thetas))),reg) + c(rowWeights*(lambda1 + lambda2* 2*rMat + rMatK %*% lambda3))

  center = sum(rMat*rowWeights)
  unitSum = sum(rMat^2*rowWeights)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(rMatK, 2, function(x){
    sum(rMat*x*rowWeights)
  })
  return(c(score,centers, unitSums, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianRmat = function(beta, X, reg, muMarg, k, n, Z, nLambda, thetas, rowWeights, rMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg) * muMarg

  Jac = matrix(0, nrow = n + nLambda, ncol=n + nLambda)
  #The suymmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:n, n+1] = rowWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:n, n+2] = 2 *rMat*rowWeights

  if(k>1){
    Jac[1:n,(n+3):(n+nLambda)] = apply(rMatK, 2, function(x){rowWeights*x})
  }

  #Symmetrize
  Jac = Jac + t(Jac)
  tmp= ((1+t(t(X)/thetas))*mu/(1+t(t(mu)/thetas))^2)*(1-Z)
  diag(Jac[1:n,1:n]) = -tcrossprod(reg^2 ,tmp) + 2*rowWeights*beta[n+2]
  Jac
}

#--------------------------------------#
dZinbMeanCmat = function(beta, X, reg, muMarg, k, p, Z, nLambda, thetas, colWeights, cMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  cMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  mu = exp(reg %*% cMat) * muMarg

  lambda1 = beta[p+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[p+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(p+3):length(beta)]}

  score =
    crossprod(reg,((1-Z)*((X-mu)/(1+t(t(mu)/thetas))))) +
    colWeights*(lambda1 + lambda2*2*cMat + lambda3 %*% cMatK)

  center = sum(colWeights*cMat)
  unitSum = sum(colWeights*cMat^2)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = apply(cMatK,1,function(x){
    sum(cMat*x*colWeights)
  })
  return(c(score,centers, unitSums, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianCmat = function(beta, X, reg, colWeights, k, p, Z, nLambda, thetas, muMarg, cMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param X: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  cMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  mu = exp(reg %*% cMat) * muMarg

  tmp= ((1+t(t(X)/thetas))*mu/(1+t(t(mu)/thetas))^2)*(1-Z)
  Jac = matrix(0, nrow = p + nLambda, ncol = p + nLambda)
  #The symmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:p,(p+1)] = colWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:p,p+2] = colWeights*2 *cMat

  #dLag²/ds_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:p,(p+3):(p+nLambda)] = apply(cMatK, 1,function(x){
      colWeights*x
    })
  }
  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:p,1:p]) = -crossprod(tmp, reg^2) + 2*beta[p+2]*colWeights
  Jac
}

#--------------------------------------#
dZinbZeroChi = function(beta, Z, reg, zeroMarg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates
  GZero = expit(arrayprod(reg, beta) + logit(zeroMarg))
  sum((Z-GZero)*reg)
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianChi = function(beta, reg, Z, zeroMarg){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes
  expZero = exp(reg* beta + logit(zeroMarg))
  tmp = (expZero/(1+expZero)^2)
  # tmp[is.infinite(expZero)]=0
  -sum(reg^2*tmp)
}

#--------------------------------------#
dZinbZeroTmat = function(beta, reg,  k, n, Z, zeroMarg, nLambda, rowWeights, tMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  tMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)

  GZero = expit(tMat %*% reg + logit(zeroMarg))

  lambda1 = beta[n+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[n+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(n+3):length(beta)]}

  score = tcrossprod((Z-GZero), reg) +
    (lambda1+ lambda2*2*tMat + tMatK %*% lambda3)*rowWeights

  center = sum(tMat*rowWeights)
  unitSum = sum(tMat^2*rowWeights)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = colSums(tMatK*tMat*rowWeights)
  return(c(score,centers, unitSums, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianTmat = function(beta, reg, k, n, Z, zeroMarg, nLambda, rowWeights, tMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param abunds: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  tMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  Jac = matrix(0,nrow = n+nLambda, ncol = n+nLambda)

  expZero = exp(tMat %*% reg + logit(zeroMarg))

  tmp = (expZero/(1+expZero)^2)
  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:n, n+1] = rowWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:n, n+2] = 2*tMat*rowWeights
  tmp=expZero/(1+expZero)^2
  #dLag²/dr_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:n,(n+3):(n+nLambda)] = apply(tMatK, 2, function(x){rowWeights*x})
  }

  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:n,1:n]) = -tcrossprod(tmp, reg^2) + 2*beta[n+2]*rowWeights

  Jac
}

#--------------------------------------#
dZinbZeroVmat = function(beta, reg, colWeights, k, p, Z, zeroMarg, nLambda, vMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  # @return A vector of length r with the new psi estimates

  vMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)

  GZero = expit(reg %*% vMat + logit(zeroMarg))

  lambda1 = beta[p+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[p+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(p+3):length(beta)]}

  score =
    crossprod(reg,(Z-GZero)) +
    colWeights*(lambda1 + lambda2*2*vMat + (lambda3 %*% vMatK))


  center = sum(colWeights*vMat)
  unitSum = sum(colWeights*vMat^2)-1
  if(k==1){return(c(score,center, unitSum))}
  orthogons = colSums(t(vMatK)*vMat*colWeights)
  return(c(score,centers, unitSums, orthogons))
}

#--------------------------------------#
#A jacobian for the psi parameters
ZinbJacobianVmat = function(beta, reg, colWeights, k, p, Z, zeroMarg, nLambda, vMatK){
  # @param beta: a vector of r regression parameters to optimize: the r psi parameters
  # @param y: the nxp data matrix
  # @param reg: a nxpxr regressor array with r the number of regressors
  # @param theta: a vector of length p with the dispersion parameters
  # @param k: a scalar, dimension of the RC solution
  # @param colWeights: a vector of length p with the abundance parameters
  # @param libSizes (optional): a vector of length n with (known) library sizes

  vMat = matrix(beta[1:p], byrow=TRUE, ncol=p, nrow=1)
  expZero = exp(reg %*% vMat +logit(zeroMarg))
  tmp = (expZero/(1+expZero)^2)
  # tmp[is.infinite(expZero)]=0

  Jac = matrix(0, nrow = p + nLambda, ncol = p + nLambda)
  #The symmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:p,(p+1)] = colWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:p,p+2] = colWeights*2 *vMat
  if(k>1){
    Jac[1:p,(p+3):(p+nLambda)] = apply(vMatK, 1,function(x){
      colWeights*x
    })
  }

  #Symmetrize
  Jac = Jac + t(Jac)
  diag(Jac[1:p,1:p]) = -crossprod(reg^2,tmp) + colWeights *beta[p+2]*2

  Jac
}

RCM_ZINB = function(X, k, rowWeights , colWeights, weightsChar, tol = 1e-3, maxItOut = 500, psiTol = 1e-4, chiTol = psiTol, verbose = TRUE, ZINBRCM = NULL, global = "dbldog", nleqslv.control = list(), method="Broyden", twoCores = FALSE, convNorm = 2, maxItMean = 20, maxItZeroes = 20, dispFreq = 5){

  # @param X: a nxp data matrix
  # @param k: a scalar, number of dimensions in the RC(M) model
  # @param tol(optional): a scalar, the relative convergende tolerance for the row scores and column scores parameters, defaults to 1e-3
  # @param Psitol(optional): a scalar, the relative convergence tolerance for the psi parameters, defaults to 1e-4
  # @param maxItOut(optional): an integer, the maximum number of iteration in the outer loop, defaults to 50
  # @param libSizes(optional) : a vector of length n with (known) library sizes. If not provided, rowSums of x are used
  # @param verbose(optional): a boolean, should information on iterations be printed? Defaults to TRUE
  # @param method(optional): Method for jacobian estimation , see nleqslv. Defaults to Broyden. The difference with the newton method is that the Jacobian is not recalculated at every iteration
  # @param global(optional): global strategy for solving non-linear systems , see nleqslv
  # @param nleqslv.control: a list with control options, see nleqslv
  # @param lambdaRow: a vector of length 2*k+k*(k-1)/2 with inital estimates or the lagrange multipliers for the row scores
  # @param lambdaCol: a vector of length 2*k+k*(k-1)/2 with inital estimates or the lagrange multipliers for the column scores
  # @param rMatInit(optional): a nxk matrix with initial row scores. If not provided values from the singular value decomposition will be used as starting values
  # @param cMatInit(optional): a pxk matrix with initial column scores. If not provided values from the singular value decomposition will be used as starting values
  # @param psisInit(optional): a vector of length k with inital values for the importance parameters psi. If not provided values from the singular value decomposition will be used as starting values
  # @param dispFreq: a scalar, how many iterations the algorithm should wait before reestimationg the dispersions
  # @param convNorm: a scalar, the norm to use to determine convergence

  # @return A list with elements:
  # @return psis: a vector of length k with estimates for the importance parameters psi
  # @return thetas: a vector of length p with estimates for the overdispersion
  # @return rMat: a nxk matrix with estimated row scores
  # @return cMat: a pxk matrix with estimated column scores
  # @return converged: a boolean indicating if the algorithm converged
  # @return rowRec: a n x k x maxItOut array with a record of all rMat estimates through the iterations
  # @return colRec: a k x p x maxItOut array with a record of all cMat estimates through the iterations
  # @return psiRec.: a k x maxItOut array with a record of all psi estimates through the iterations

  libSizes = rowSums(X)
  abunds = (colSums(X)/sum(X))
  n=NROW(X)
  p=NCOL(X)
  logLibSizesMLE = log(libSizes)
  logAbundsMLE = log(abunds)
  # logitZeroRows = logit(rep(0.25,n)) #Start with a 25 % chance on a structural zero
  logitZeroCols = logit(rep(0.25,p))
  iterOut = rep(1,k)

  initIter = 1

  while((initIter ==1) || ((initIter <= maxItOut) && (!convergenceInit))){
    cat("Starting iteration ",initIter, " \n")

    libsOld = logLibSizesMLE
    absOld = logAbundsMLE
    # zeroRowsOld = logitZeroRows
    zeroColsOld = logitZeroCols
    Z = matrix(0,n,p)

    #M-step
    initIterMean = 1
    while((initIterMean ==1) || ((initIterMean <= maxItMean) && (!convergenceInitMean))){
      cat("Mean iteration ",initIterMean, "\n")
      libsOldIn = logLibSizesMLE
      absOldIn = logAbundsMLE
      thetas = estDisp(X = X, cMat = matrix(0,1,p), rMat = matrix(0,n,1),  muMarg=exp(outer(logLibSizesMLE, logAbundsMLE, "+")), psis = 0, dispWeights=t(1-Z))
      libsTmp = try(nleqslv(fn = dZinbMeanLibsizes, x = logLibSizesMLE, X = X, reg=logAbundsMLE, global=global, control = nleqslv.control, jac=ZinbJacobianLibsizes, method=method, Z=Z, thetas = thetas)$x, silent=TRUE)
      if(class(libsTmp)!="try-error"){ logLibSizesMLE = libsTmp}
      absTmp = try(nleqslv(fn = dZinbMeanAbunds, x = logAbundsMLE, X = X, reg=logLibSizesMLE, global=global, control = nleqslv.control, jac=ZinbJacobianAbunds, method=method, Z=Z, thetas = thetas)$x, silent=TRUE)
      if(class(absTmp)!="try-error"){logAbundsMLE = absTmp}

      initIterMean = initIterMean + 1

      convergenceInitMean = ((initIterMean <= maxItMean) &&
                               ((sum(abs(1-logLibSizesMLE/libsOldIn)^convNorm)/n)^(1/convNorm) < tol) &&
                               ((sum(abs(1-logAbundsMLE/absOldIn)^convNorm)/p)^(1/convNorm) < tol))
    }
    #     initIterZero = 1
    #   while((initIterZero ==1) || ((initIterZero <= maxItZeroes) && (!convergenceInitZeroes))){
    #         cat("Zero iteration ",initIterZero, "\n")
    #     zeroRowsOldIn = logitZeroRows
    #     zeroColsOldIn = logitZeroCols

    #   zeroRowsTmp = try(nleqslv(fn = dZinbZeroRow, x = logitZeroRows,   reg=logitZeroCols, global=global, control = nleqslv.control, jac=ZinbJacobianZeroRow, method=method, Z=Z)$x, silent=TRUE)
    #  if(class(zeroRowsTmp)!="try-error"){ logitZeroRows = zeroRowsTmp}
    zeroColsTmp = try(nleqslv(fn = dZinbZeroCol, x = logitZeroCols, reg = rep.int(0L, n), global=global, control = nleqslv.control, jac=ZinbJacobianZeroCol, method=method, Z=Z)$x, silent=TRUE)
    if(class(zeroColsTmp)!="try-error"){ logitZeroCols = zeroColsTmp}

    #     initIterZero = initIterZero + 1
    #
    #    convergenceInitZeroes = ((initIterZero <= maxItZeroes) &&
    #                    # ((sum(abs(1-logitZeroRows/zeroRowsOldIn)^convNorm)/p)^(1/convNorm) < tol) &&
    #                    ((sum(abs(1-logitZeroCols/zeroColsOldIn)^convNorm)/p)^(1/convNorm) < tol))
    #   }

    #E-step
    expMu = exp(outer(logLibSizesMLE, logAbundsMLE, "+"))
    pZero = expit(outer(logitZeroRows, logitZeroCols, "+"))
    thetaMat=matrix(thetas, ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    d0=dnbinom(0,mu=expMu, size=thetaMat)
    Z[X==0] = (pZero/((1-pZero)*d0 + pZero))[X==0]

    initIter = initIter + 1
    convergenceInit = ((initIter <= maxItOut) &&
                         ((sum(abs(1-logLibSizesMLE/libsOld)^convNorm)/n)^(1/convNorm) < tol) &&
                         ((sum(abs(1-logAbundsMLE/absOld)^convNorm)/p)^(1/convNorm) < tol) &&
                         # ((sum(abs(1-logitZeroRows/zeroRowsOld)^convNorm)/p)^(1/convNorm) < tol) &&
                         ((sum(abs(1-logitZeroCols/zeroColsOld)^convNorm)/p)^(1/convNorm) < tol))
  }
  muMarg = exp(outer(logLibSizesMLE, logAbundsMLE, "+"))
  zeroMarg = expit(outer(logitZeroRows, logitZeroCols, "+"))#The marginals to be used as expectation. These are augmented with the previously estimated dimensions every time

  rowRec = rowRecZeroes = array(0,dim=c(n,k, maxItOut))
  colRec = colRecZeroes = thetaRec = array(0,dim=c(k,p, maxItOut))
  psiRec = matrix(0, nrow=k,ncol=maxItOut)
  convergence = rep(FALSE, k)
  iterOut = rep(1,k)

  #If previous fit provided with higher or equal dimension, stop here
  if((!is.null(ZINBRCM)) ){
    if(ZINBRCM$fit != "RCM_ZINB"){
      stop("Fit provided is not of same type as the one requested! \n")
    } else if((k <= ZINBRCM$k)) {
      # stop("Fit provided is already of the required dimension or higher! \n")
    } else{
      for(i in c("rMat","cMat","psis","lambdaCol","lambdaRow", "lambdaRowZero","lambdaColZero","tMat","vMat","chis","Z","zeroMarg","thetas")){
        assign(i, ZINBRCM[[i]])
      }
    }
    #Otherwise try to use intelligent starting values
  } else{
    #Use a more heavily weighted
    svdX = svd(diag(1/libSizes) %*% (X-muMarg)*(1-zeroMarg) %*% diag(1/colSums(X)))

    rMat = svdX$u[,1:k,drop=FALSE]
    cMat = t(svdX$v[,1:k,drop=FALSE])
    psis = svdX$d[1:k]
    #
    #   #Redistribute some weight to fit the constraints
    #   psis = c(psis *t(apply(cMat, 1, function(colS){
    #       sqrt(sum(colWeights * colS^2))
    #   })) * apply(rMat, 2, function(rowS){
    #       sqrt(sum(rowWeights * rowS^2))
    #   }))
    #
    # #Normalize
    # cMat = t(apply(cMat, 1, function(colS){
    #       colS/sqrt(sum(colWeights * colS^2))
    #   }))
    # rMat = apply(rMat, 2, function(rowS){
    #       rowS/sqrt(sum(rowWeights * rowS^2))
    #   })
    #
    #   #Initial estimates for zeroes is also based on an svd
    #
    Xzeroes = X==0

    svdZero = svd(diag(sqrt(1/libSizes)) %*%(Xzeroes-zeroMarg)%*% diag(1/sqrt(colSums(X))))

    tMat = svdZero$u[,1:k, drop=FALSE]
    vMat = t(svdZero$v[,1:k, drop=FALSE])
    chis = svdZero$d[1:k]
    #
    # #Redistribute some weight to fit the constraints
    # chis = c(chis *t(apply(vMat, 1, function(colS){
    #       sqrt(sum(colWeights * colS^2))
    #   })) * apply(tMat, 2, function(rowS){
    #       sqrt(sum(rowWeights * rowS^2))
    #   }))
    #
    # #Normalize
    # vMat = t(apply(vMat, 1, function(colS){
    #       colS/sqrt(sum(colWeights * colS^2))
    #   }))
    # tMat = apply(tMat, 2, function(rowS){
    #       rowS/sqrt(sum(rowWeights * rowS^2))
    #   })

    #     tMat = rMat = matrix(1/n, n, k)
    #     vMat = cMat = t(matrix(1/p, p, k))
    #     psis = chis = rep(1,k)

    lambdaRow = lambdaCol = lambdaColZero=lambdaRowZero = rep.int(0,2*k+k*(k-1)/2)
  }

  rowRec = rowRecZero = array(0,dim=c(NROW(X),k, maxItOut))
  colRec = colRecZero = array(0,dim=c(k,NCOL(X), maxItOut))
  psiRec = chiRec = matrix(0,ncol=maxItOut, nrow=k)

  if(!is.null(ZINBRCM)){ #If fit provided, replace lower dimension starting values
    Kprev = ZINBRCM$k
    rMat[,1:Kprev] = ZINBRCM$rMat
    rowRec[,1:Kprev,] = ZINBRCM$rowRec
    cMat[1:Kprev,] = ZINBRCM$cMat
    colRec[1:Kprev,,] = ZINBRCM$colRec
    psis[1:Kprev] = ZINBRCM$psis
    psiRec[1:Kprev,] = ZINBRCM$psiRec
    lambdaCol[1:(Kprev*(2+(Kprev-1)/2))] = ZINBRCM$lambdaCol
    lambdaRow[1:(Kprev*(2+(Kprev-1)/2))] = ZINBRCM$lambdaRow
    tMat[,1:Kprev] = ZINBRCM$tMat
    rowRecZero[,1:Kprev,] = ZINBRCM$rowRecZero
    vMat[1:Kprev,] = ZINBRCM$vMat
    colRecZero[1:Kprev,,] = ZINBRCM$colRecZero
    chis[1:Kprev] = ZINBRCM$chis
    chiRec[1:Kprev,] = ZINBRCM$chiRec
    lambdaColZero[1:(Kprev*(2+(Kprev-1)/2))] = ZINBRCM$lambdaColZero
    lambdaRowZero[1:(Kprev*(2+(Kprev-1)/2))] = ZINBRCM$lambdaRowZero
    convergence[1:Kprev] = ZINBRCM$converged
    iterOut[1:Kprev] = ZINBRCM$iter
    zeroMarg = ZINBRCM$zeroMarg
  }

  minK = ifelse(is.null(ZINBRCM),1,Kprev+1)
  for (KK in minK:k){

    cat("Dimension" ,KK, "is being esimated \n")

    #Modify offsets if needed
    if(KK>1){
      muMarg = muMarg * exp(rMat[,(KK-1), drop=FALSE] %*% (cMat[(KK-1),, drop=FALSE]*psis[(KK-1)]))
      zeroMarg = expit(logit(zeroMarg) + tMat[,(KK-1), drop=FALSE] %*% (vMat[(KK-1),, drop=FALSE]*chis[(KK-1)]))
    }
    #A lambda parameter
    nLambda = KK + 1

    #The location of the lambda parameters
    idK = seq_k(KK)
    Z = matrix(0, n, p)

    ## 2) Propagation

    while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) && (!convergence[KK])))
    {

      if(verbose && iterOut[KK]%%1 == 0){
        cat("\n","Outer Iteration", iterOut[KK], "\n","\n")
        if(iterOut[KK]!=1){
          cat("Old psi-estimate: ", psiOld, "\n")
          cat("New psi-estimate: ", psis[KK], "\n")
        }
      }
      ## 2)a. Store old parameters
      psiOld = psis[KK]
      rMatOld = rMat[,KK]
      cMatOld = cMat[KK,]

      chiOld = chis[KK]
      tMatOld = tMat[,KK]
      vMatOld = vMat[KK,]

      #Expectation
      Z = EstepNB (X = X, rMat = rMat[,KK, drop=FALSE], cMat = cMat[KK,, drop=FALSE], tMat= tMat[,KK, drop=FALSE], vMat = vMat[KK,, drop=FALSE], muMarg = muMarg, zeroMarg = zeroMarg, psis = psis[KK], chis = chis[KK], thetas =thetas)

      #Maximization
      Mlist = MstepNB(Z = Z, X = X, rMat = rMat[,KK, drop=FALSE], cMat = cMat[KK,, drop=FALSE], tMat = tMat[,KK, drop=FALSE], vMat = vMat[KK,, drop=FALSE], k = KK, n=n, p=p,muMarg = muMarg, zeroMarg = zeroMarg, psis = psis [KK],chis = chis[KK], twoCores=twoCores, tol = tol, psiTol = psiTol, chiTol = chiTol, convNorm = convNorm, global=global, nLambda = nLambda, nleqslv.control = nleqslv.control, lambdaCol = lambdaCol[idK], lambdaRow=lambdaRow[idK], lambdaColZero = lambdaColZero[idK], lambdaRowZero = lambdaRowZero[idK], maxItMean = maxItMean, maxItZeroes = maxItZeroes, dispFreq=dispFreq, thetas = thetas, rowWeights = rowWeights, colWeights = colWeights, rMatK = rMat[,1:(KK-1), drop=FALSE], cMatK = cMat[1:(KK-1),, drop=FALSE], tMatK = tMat[,1:(KK-1), drop=FALSE], vMatK = vMat[1:(KK-1),, drop=FALSE] )

      #Assign outcomes to tracking vectors and to this environment
      lambdaCol[idK] = Mlist$lambdaCol
      lambdaRow[idK] = Mlist$lambdaRow
      lambdaColZero[idK] = Mlist$lambdaColZero
      lambdaRowZero[idK] = Mlist$lambdaRowZero

      rowRec[,KK, iterOut[KK]] = rMat[,KK] = Mlist$rMat
      colRec[KK,, iterOut[KK]] = cMat[KK,] = Mlist$cMat
      rowRecZero[,KK, iterOut[KK]] = tMat[,KK] = Mlist$tMat
      colRecZero[KK,, iterOut[KK]] = vMat[KK,] = Mlist$vMat
      psiRec[KK, iterOut[KK]] = psis[KK] = Mlist$psis
      chiRec[KK, iterOut[KK]] = chis[KK] = Mlist$chis

      ## 2)f. Change iterator
      iterOut[KK] = iterOut[KK] + 1

      ##Check convergence  (any numbered norm for row and column scores)
      convergence[KK] = (iterOut[KK] <= maxItOut) &&
        (all(abs(1-psis/psiOld) < psiTol)) &&
        ((sum((1-rMatOld/rMat)^convNorm)/n)^(1/convNorm) < tol) &&
        ((sum((1-cMatOld/cMat)^convNorm)/p)^(1/convNorm) < tol)  &&
        (all(abs(1-chis/chiOld) < chiTol)) &&
        (sum(abs(1-tMat/tMatOld)^convNorm)/n)^(1/convNorm) < tol &&
        (sum(abs(1-vMat/vMatOld)^convNorm)/p)^(1/convNorm) < tol
    } # END while-loop
  } # End for-loop

  ## 3) Termination
  rownames(rMat) = rownames(tMat) = rownames(X)
  colnames(cMat) = colnames(vMat) = colnames(X)
  rownames(cMat) = rownames(vMat) = colnames(tMat) = colnames(rMat) = paste0("Dim",1:k)

  if(!convergence ){
    warning("Algorithm did not converge! Check for errors or consider changing tolerances or number of iterations")
  }
  return(list(converged = convergence, rMat=rMat, cMat=cMat, psis = psis, X=X,
              rowRec = rowRec, colRec = colRec, psiRec = psiRec, lambdaRow = lambdaRow, lambdaCol = lambdaCol, lambdaRowZero = lambdaRowZero, lambdaColZero = lambdaColZero, chis = chis, tMat = tMat, vMat = vMat, chiRec = chiRec, rowRecZero = rowRecZero, colRecZero = colRecZero, iter=iterOut-1, Z=Z, thetas=thetas, libSizesMLE = exp(logLibSizesMLE), abundsMLE = exp(logAbundsMLE), taxaZeroes = expit(logitZeroCols), rowWeights=rowWeights, colWeights=colWeights, iter=iterOut-1, fit="RCM_ZINB"))
}
```

### Zero-inflated poisson

  #### ZIP without signal

  ##### Generate data

  ```{r ZIP without signal, eval=FALSE, purl = FALSE}
load("/home/stijn/PhD/American Gut/AGphylo.RData")
#First estimate the ZIP parameters

if(!file.exists("AGzipParams.RData")){
  otuTab = otu_table(AGphylo)@.Data
  logLibs = log(sample_sums(AGphylo))
  ZIPfits = mclapply(mc.cores=4,1:ncol(otuTab),  function(i){
    zeroinfl(otuTab[,i]~offset(logLibs)| 1)
  })
  zeroProbs = sapply(ZIPfits, function(x){expit(x$coef$zero)})
  ZIPmeans = sapply(ZIPfits, function(x){exp(x$coef$count)})
  save(zeroProbs, ZIPmeans, file="AGzipParams.RData")
} else {load(file="AGzipParams.RData")}
#ZeroProbs are very high

#Define parameters
zeroProbs = zeroProbs[zeroProbs<0.9]
NsamplesZIPnoSig= 200
NtaxaZIPnoSig = 800
idSampleZIP = sample(size=NtaxaZIPnoSig, 1:length(ZIPmeans))
lambdasZIPnoSig = ZIPmeans[idSampleZIP]
lambdasZIPnoSig= lambdasZIPnoSig/sum(lambdasZIPnoSig)
zeroesZIPnoSig = sample(zeroProbs, NtaxaZIPnoSig)
libSizesZIPnoSig =c(rep(1e5, floor(NsamplesZIPnoSig/2)), rep(1e6, floor(NsamplesZIPnoSig/2)))

#Mean and zero matrices
meanMatZIPnoSig = outer(libSizesZIPnoSig, lambdasZIPnoSig)
zeroMatZIPnoSig = matrix(zeroesZIPnoSig, nrow=NsamplesZIPnoSig, ncol=NtaxaZIPnoSig, byrow = TRUE)

#Data generation
dataMatZIPnoSig = matrix(rzipois(n=prod(dim(meanMatZIPnoSig)),lambda = meanMatZIPnoSig, pstr0 = zeroMatZIPnoSig), ncol=NtaxaZIPnoSig, nrow=NsamplesZIPnoSig)
save(dataMatZIPnoSig, file = "zipNSData.RData")
```

##### Fit RC(M)

```{r RC(M) fit, eval=FALSE, purl = FALSE}
nleqslv.controlZIP = list(trace=TRUE, maxit=250, cndtol = .Machine$double.eps, allowSingular=TRUE)
if(!file.exists("syntZIPnoSig.RData")){
  syntZIPnoSigUnifUnifJob = mcparallel(RCM(dataMatZIPnoSig, method = "ZIP", k=2, nleqslv.control = nleqslv.controlZIP, colWeights = "uniform", rowWeights = "uniform", maxItOut=5e0))
  syntZIPnoSigUnifUnif =  mccollect(syntZIPnoSigUnifUnifJob, FALSE)[[1]]

  syntZIPnoSigUnifMargJob = mcparallel(RCM(dataMatZIPnoSig, method = "ZIP",k=2, nleqslv.control = nleqslv.controlZIP, colWeights = "marginal", rowWeights = "uniform", maxItOut=5e2))
  syntZIPnoSigUnifMarg =  mccollect(syntZIPnoSigUnifMargJob, FALSE)[[1]]

  syntZIPnoSigMargUnifJob = mcparallel(RCM(dataMatZIPnoSig, method = "ZIP",k=2, nleqslv.control = nleqslv.controlZIP, colWeights = "uniform", rowWeights = "marginal", maxItOut=5e2))
  syntZIPnoSigMargUnif =  mccollect(syntZIPnoSigMargUnifJob, FALSE)[[1]]

  syntZIPnoSigMargMargJob = mcparallel(RCM(dataMatZIPnoSig, method = "ZIP",k=2, nleqslv.control = nleqslv.controlZIP, colWeights = "marginal", rowWeights = "marginal", maxItOut=5e1))
  syntZIPnoSigMargMarg =  mccollect(syntZIPnoSigMargMargJob, FALSE)[[1]]

  save(syntZIPnoSigUnifMarg, syntZIPnoSigUnifUnif, syntZIPnoSigMargMarg, syntZIPnoSigMargUnif, dataMatZIPnoSig, file="syntZIPnoSig.RData")
} else {load("syntZIPnoSig.RData")}
#No convergence, mean cannot be fitted
```

##### Plot the results

```{r RC(M) zip, eval=FALSE, purl = FALSE}
load(file="toyDataNoSigZIPserver.RData")
solListWSZIP = noSigListZIP[sapply(noSigListZIP, class)=="list"]
#Runtimes and convergence
sapply(solListWSZIP,function(x){x$converged})
sapply(solListWSZIP,function(x){x$runtime})
sapply(solListWSZIP,function(x){x$iter})

par(mfrow=c(2,2))
sapply(solListWSZIP,function(x){with(x, plot(psiRec[1,1:iter[1]]))})
```

#### ZIP with signal

```{r ZIP with signal, eval=FALSE, purl = FALSE}
load("/home/stijn/PhD/American Gut/AGphylo.RData")
#First estimate/load the ZIP parameters

if(!file.exists("AGzipParams.RData")){
  otuTab = otu_table(AGphylo)@.Data
  logLibs = log(sample_sums(AGphylo))
  ZIPfits = mclapply(mc.cores=4,1:ncol(otuTab),  function(i){
    zeroinfl(otuTab[,i]~offset(logLibs)| 1)
  })
  zeroProbs = sapply(ZIPfits, function(x){expit(x$coef$zero)})
  ZIPmeans = sapply(ZIPfits, function(x){exp(x$coef$count)})
  save(zeroProbs, ZIPmeans, file="AGzipParams.RData")
} else {load(file="AGzipParams.RData")}
#ZeroPorbs are very low

#Defina parameters
NsamplesZIPsig= 150
NtaxaZIPsig = 500
lambdasZIPref = lambdasZIPsig1 =lambdasZIPsig2 = lambdasZIPsig12 =  sample(ZIPmeans, NtaxaZIPsig)
zeroesZIPref = zeroesZIPsig1 = zeroesZIPsig2 = zeroesZIPsig12 =sample(zeroProbs, NtaxaZIPsig)
libSizesZIPsig =c(rep(1e4, floor(NsamplesZIPsig/2)), rep(1e5, floor(NsamplesZIPsig/2)))

#Define the signal
NtaxaSignal1ZIP = 20
NtaxaSignal2ZIP = 20

NtaxaSignal1zeroZIP = 20
NtaxaSignal2zeroZIP = 20

Signal1taxZIP = 4
Signal2taxZIP = 3

Signal1zeroTaxZIP = 5
Signal2zeroTaxZIP = 3

NsamplesSignal1ZIP = 20
NsamplesSignal2ZIP = 20
NsamplesSignal12ZIP = 20

NsamplesSignal1zeroZIP = 20
NsamplesSignal2zeroZIP = 20
NsamplesSignal12zeroZIP = 20

idSig1TaxZIP = 1:NtaxaSignal1ZIP
idSig2TaxZIP = sample(1:NtaxaZIPsig,NtaxaSignal2ZIP)

lambdasZIPsig1[idSig1TaxZIP] = lambdasZIPsig1[idSig1TaxZIP] * Signal1taxZIP
lambdasZIPsig2[idSig2TaxZIP] = lambdasZIPsig2[idSig2TaxZIP] * Signal2taxZIP
lambdasZIPsig12[c(idSig1TaxZIP, idSig2TaxZIP)] = c(lambdasZIPsig1[1:NtaxaSignal1ZIP],lambdasZIPsig2[idSig2TaxZIP])

idSig1TaxZIPzeroes = 1:NtaxaSignal1zeroZIP + NtaxaSignal1ZIP/2
idSig2TaxZIPzeroes = idSig1TaxZIPzeroes + NtaxaSignal1zeroZIP/2 + NtaxaSignal2zeroZIP/2

zeroesZIPsig1[idSig1TaxZIPzeroes] = expit(logit(zeroesZIPsig1[idSig1TaxZIPzeroes]) * Signal1zeroTaxZIP)
zeroesZIPsig2[idSig2TaxZIPzeroes] = expit(logit(zeroesZIPsig2[idSig2TaxZIPzeroes]) * Signal1zeroTaxZIP)
zeroesZIPsig12[c(idSig1TaxZIPzeroes, idSig2TaxZIPzeroes)] = c(zeroesZIPsig1[idSig1TaxZIPzeroes],zeroesZIPsig2[idSig2TaxZIPzeroes])

NsamplesZIPsigRef =  NsamplesZIPsig - NsamplesSignal1ZIP - NsamplesSignal2ZIP-  NsamplesSignal12ZIP

#Define mean and zero matrices
meanMatrefZIP = outer(sample(libSizesZIPsig, NsamplesZIPsigRef), lambdasZIPref)
meanMatSig1ZIP = outer(sample(libSizesZIPsig, NsamplesSignal1ZIP), lambdasZIPsig1)
meanMatSig2ZIP = outer(sample(libSizesZIPsig, NsamplesSignal2ZIP), lambdasZIPsig2)
meanMatSig12ZIP = outer(sample(libSizesZIPsig, NsamplesSignal12ZIP), lambdasZIPsig12)

zeroMatrefZIP = matrix(zeroesZIPref, nrow=nrow(meanMatrefZIP),ncol=ncol(meanMatrefZIP), byrow = TRUE)
zeroMatSig1ZIP = matrix(zeroesZIPsig1, nrow=nrow(meanMatSig1ZIP),ncol=ncol(meanMatSig1ZIP), byrow = TRUE)
zeroMatSig2ZIP = matrix(zeroesZIPsig2, nrow=nrow(meanMatSig2ZIP),ncol=ncol(meanMatSig2ZIP), byrow = TRUE)
zeroMatSig12ZIP = matrix(zeroesZIPsig12, nrow=nrow(meanMatSig12ZIP),ncol=ncol(meanMatSig12ZIP), byrow = TRUE)

#Generate the data
dataMatZIPref = matrix(rzipois(n=prod(dim(meanMatrefZIP)),lambda = meanMatrefZIP, pstr0 = zeroMatrefZIP), ncol=ncol(meanMatrefZIP), nrow=nrow(meanMatrefZIP))
dataMatZIPsig1 = matrix(rzipois(n=prod(dim(meanMatSig1ZIP)),lambda = meanMatSig1ZIP, pstr0 = zeroMatSig1ZIP), ncol=ncol(meanMatSig1ZIP), nrow=nrow(meanMatSig1ZIP))
dataMatZIPsig2 = matrix(rzipois(n=prod(dim(meanMatSig2ZIP)),lambda = meanMatSig2ZIP, pstr0 = zeroMatSig2ZIP), ncol=ncol(meanMatSig2ZIP), nrow=nrow(meanMatSig2ZIP))
dataMatZIPsig12 = matrix(rzipois(n=prod(dim(meanMatSig12ZIP)),lambda = meanMatSig12ZIP, pstr0 = zeroMatSig12ZIP), ncol=ncol(meanMatSig12ZIP), nrow=nrow(meanMatSig12ZIP))

#Bind the data
dataMatZIPSig = rbind(dataMatZIPref, dataMatZIPsig1, dataMatZIPsig2, dataMatZIPsig12)
#mean(dataMatZIPsig==0) #Correct zero fraction

#Save signals
sampleSigZIP = factor(c(rep("Reference",NsamplesZIPsigRef), rep("Signal1", NsamplesSignal1ZIP), rep("Signal2", NsamplesSignal2ZIP), rep("Signal 1 and 2", NsamplesSignal12ZIP)))
taxaSigZIPtmp = rep("Reference",ncol(dataMatSigNB))
taxaSigZIPtmp[idSig1TaxZIP] = "Signal 1"
taxaSigZIPtmp[idSig2TaxZIP] = "Signal 2"
taxaSigZIP =factor(taxaSigZIPtmp)
names(taxaSigZIP) = colnames(dataMatZIPSig) = names(rhosNBSigref)
names(sampleSigZIP) = rownames(dataMatZIPSig) = 1:NsamplesZIPsig
```

#### Apply the RC(M) algorithm

```{r ZIP signal fit, eval=FALSE, purl=FALSE}
if(!file.exists("syntZIPSig.RData")){
  syntZIPSigUnifJob = mcparallel(RCM(dataMatZIPSig, method = "ZIP",k=2, nleqslv.control = list(trace=TRUE, maxit=250), colWeights = "uniform"))
  syntZIPSigUnif =  mccollect(syntZIPSigUnifJob, FALSE)[[1]]

  syntZIPSigMargJob = mcparallel(RCM(dataMatZIPSig, method = "ZIP",k=2, nleqslv.control = list(trace=TRUE, maxit=250), colWeights = "marginal"))
  syntZIPSigMarg =  mccollect(syntZIPSigMargJob, FALSE)[[1]]

  save(syntZIPSigUnif, syntZIPSigMarg, sampleSigZIP, taxaSigNB, file="syntZIPSig.RData")
} else {load("syntZIPSig.RData")}
```

##### Plot the results

### Zero-inflated negative binomial

#### Without signal

##### Generate data

```{r ZINB without signal, eval=FALSE, purl=FALSE}
load("/home/stijn/PhD/American Gut/AGphylo.RData")
#First estimate the ZIP parameters
expit=function(x){
  tmp = exp(x)/(1+exp(x))
  tmp[is.na(tmp)]=1 #Adjust for overflow
  tmp}

if(!file.exists("AGzinbParams.RData")){
  otuTab = otu_table(AGphylo)@.Data
  # otuTab = otuTab[, colMeans(otuTab==0)<0.95]
  logLibs = log(sample_sums(AGphylo))
  ZINBfits = mclapply(mc.cores=4,1:ncol(otuTab),  function(i){
    try(zeroinfl(otuTab[,i]~offset(logLibs)| 1, dist="negbin"), silent=TRUE)
  })
  zeroProbsNB = sapply(ZINBfits, function(x){if(is.list(x)) expit(x$coef$zero) else NA})
  ZINBmeans = sapply(ZINBfits, function(x){if(is.list(x)) exp(x$coef$count) else NA})
  thetasZINB = sapply(ZINBfits, function(x){if(is.list(x)) x$theta else NA})
  naID=is.na(zeroProbsNB)
  zeroProbsNB= zeroProbsNB[!naID]
  ZINBmeans= ZINBmeans[!is.na(ZINBmeans)]
  thetasZINB= thetasZINB[!is.na(thetasZINB)]
  save(thetasZINB, zeroProbsNB, ZINBmeans, file="AGzinbParams.RData")
} else {load(file="AGzinbParams.RData")}
#ZeroProbs are lower than for the ZIP

#Define parameters
zeroProbsNB = zeroProbsNB[zeroProbsNB<0.9]
NsamplesZINBnoSig= 200
NtaxaZINBnoSig = 800
idSampleZINB = sample(size=NtaxaZINBnoSig, 1:length(ZINBmeans))
MeansZINB = ZINBmeans[idSampleZINB]
ThetasZINB = thetasZINB[idSampleZINB]
zeroesZINBnoSig = sample(zeroProbsNB, NtaxaZINBnoSig)
libSizesZINBnoSig =c(rep(1e5, floor(NsamplesZINBnoSig/2)), rep(1e6, floor(NsamplesZINBnoSig/2)))

#Mean, OD and zero matrices
meanMatZINBnoSig = outer(libSizesZINBnoSig, MeansZINB)
ODmatZINB = outer(rep(1,NsamplesZINBnoSig) ,ThetasZINB )
zeroMatZINBnoSig = matrix(zeroesZINBnoSig, nrow=NsamplesZINBnoSig, ncol=NtaxaZINBnoSig, byrow = TRUE)

#Data generation
dataMatZINBnoSig = matrix(rzinegbin(n=prod(dim(meanMatZINBnoSig)),munb = meanMatZINBnoSig, size = ODmatZINB, pstr0 = zeroMatZINBnoSig), ncol=NtaxaZINBnoSig, nrow=NsamplesZINBnoSig)
save(dataMatZINBnoSig, file = "zinbNoSigData.RData")
```

##### Apply the RC(M) algorithm with the ZINB distribution

```{r ZINB signal fit, eval=FALSE, purl=FALSE}
nleqslv.control.zinb = list(trace = TRUE, maxit = 250, cndtol = .Machine$double.eps, allowSingular=TRUE)
if(!file.exists("syntZINBnoSig.RData")){
  syntZINBSigunifunifJob = mcparallel(RCM(dataMatZINBnoSig, method = "ZINB", k=3, nleqslv.control = nleqslv.control.zinb, colWeights = "uniform", rowWeights = "uniform", twoCores=FALSE, maxItOut = 5e0))
  syntZINBSigunifunif =  mccollect(syntZINBSigunifunifJob, FALSE)[[1]]

  syntZINBSigunifmargJob = mcparallel(RCM(dataMatZINBnoSig, method = "ZINB", k=3, nleqslv.control = nleqslv.control.zinb, colWeights = "marginal", rowWeights = "uniform", twoCores=FALSE))
  syntZINBSigunifmarg =  mccollect(syntZINBSigunifmargJob, FALSE)[[1]]

  save(syntZINBSigunifunif, syntZINBSigunifmarg, file="syntZINBnoSig.RData")
} else {load("syntZINBnoSig.RData")}
```

##### Plot the ZINB results

```{r RC(M) ZINB, eval=FALSE, purl=FALSE}
load(file="toyDataNoSigZINBserver.RData")
solListWSZINB = noSigListZINB[sapply(noSigListZINB, class)=="list"]
#Runtimes and convergence
sapply(solListWSZINB,function(x){x$converged})
sapply(solListWSZINB,function(x){x$runtime})
sapply(solListWSZINB,function(x){x$iter})
par(mfrow=c(2,2))
sapply(solListWSZINB,function(x){with(x, plot(psiRec[1,1:iter[1]]))})
```

#### Restore defaults

```{r restore defaults, eval=FALSE, purl = FALSE}
palette(palStore)
```

--->

  <!---
  ## Toy data

  ### Negative binomial

  We generate some data as before with the NB distribution but differing library sizes, apply our algorithm and plot the results.

#### NB without signal

##### Create the data

```{r NB no signal, purl = FALSE, eval = FALSE}
#Negative binomial, no signal. Set parameters
NsamplesNBNS = 300
NtaxaNBNS = 1100
thetasNBNS = sample(thetas,NtaxaNBNS)
thetasNBNS = thetasNBNS[1/thetasNBNS<60]
rhosNBNS = rhos[names(thetasNBNS)]
rhosNBNS = rhosNBNS/sum(rhosNBNS)
NtaxaNBNS = length(rhosNBNS)
libSizesNBNS = c(rep(1e4, floor(NsamplesNBNS/2)), rep(1e5, floor(NsamplesNBNS/2)))

#Create means and overdispersion matrices
meanNBNS = outer(libSizesNBNS, rhosNBNS)
thetaMatNBNS =  matrix(thetasNBNS, nrow=NsamplesNBNS, ncol=NtaxaNBNS, byrow=TRUE)

#Define a function to make NB data
makeNBdata = function(meanMat, thetaMat){apply(array(data = c(meanMat, thetaMat), dim=c(nrow(meanMat), ncol(meanMat), 2)), c(1,2), function(x){rnbinom(1,mu=x[1], size=x[2])})}

#Generate the data
dataMatNBNS = makeNBdata(meanNBNS, thetaMatNBNS)
```

##### Fit the RC(M) model

```{r fit RCM, purl = FALSE}
#Set control parameters
nleqslv.control = list(trace=FALSE, maxit = 500, cndtol=.Machine$double.eps)
#Fit the RC(M) model
if(!file.exists("toyDataNS.RData")){

  syntNBNSmargmarg_3Job = mcparallel(RCM(dataMatNBNS, distribution="NB", k=3, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "marginal",rowWeights ="marginal", prevCutOff=0.01))
  syntNBNSmargmarg_3 = mccollect(syntNBNSmargmarg_3Job, FALSE)[[1]]

  syntNBNSunifmarg_3Job = mcparallel(RCM(dataMatNBNS, distribution="NB", k=3, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "marginal",rowWeights ="uniform",prevCutOff=0.01))
  syntNBNSunifmarg_3 = mccollect(syntNBNSunifmarg_3Job, FALSE)[[1]]

  syntNBNSunifunif_3Job = mcparallel(RCM(dataMatNBNS, distribution="NB", k=3, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "uniform",rowWeights ="uniform", prevCutOff=0.01))
  syntNBNSunifunif_3 = mccollect(syntNBNSunifunif_3Job, FALSE)[[1]]

  syntNBNSmargunif_3Job = mcparallel(RCM(dataMatNBNS, distribution="NB", k=3, nleqslv.control= nleqslv.control, maxItOut=1e3, rowWeights = "marginal",colWeights ="uniform", prevCutOff=0.01))
  syntNBNSmargunif_3 = mccollect(syntNBNSmargunif_3Job, FALSE)[[1]]

  syntNBNSmargmarg_2JobMLE = mcparallel(RCM(dataMatNBNS, distribution="NB", k=2, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "marginal",rowWeights ="marginal", prevCutOff=0.01, marginEst = "MLE"))
  syntNBNSmargmarg_2MLE = mccollect(syntNBNSmargmarg_2JobMLE, FALSE)[[1]]

  syntNBNSunifmarg_2JobMLE = mcparallel(RCM(dataMatNBNS, distribution="NB", k=2, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "marginal",rowWeights ="uniform",prevCutOff=0.01, marginEst = "MLE"))
  syntNBNSunifmarg_2MLE = mccollect(syntNBNSunifmarg_2JobMLE, FALSE)[[1]]

  syntNBNSunifunif_2JobMLE = mcparallel(RCM(dataMatNBNS, distribution="NB", k=2, nleqslv.control= nleqslv.control, maxItOut=1e3, colWeights = "uniform",rowWeights ="uniform", prevCutOff=0.01, marginEst = "MLE"))
  syntNBNSunifunif_2MLE = mccollect(syntNBNSunifunif_2JobMLE, FALSE)[[1]]

  syntNBNSmargunif_2JobMLE = mcparallel(RCM(dataMatNBNS, distribution="NB", k=2, nleqslv.control= nleqslv.control, maxItOut=1e3, rowWeights = "marginal",colWeights ="uniform", prevCutOff=0.01, marginEst = "MLE"))
  syntNBNSmargunif_2MLE = mccollect(syntNBNSmargunif_2JobMLE, FALSE)[[1]]

  #Save the results
  save( syntNBNSunifmarg_3,syntNBNSmargunif_3,syntNBNSmargmarg_3,syntNBNSunifunif_3, dataMatNBNS,test1, test2, test3, test1unif, test2unif, test3unif, syntNBNSunifunif_2MLE, syntNBNSunifmarg_2MLE, syntNBNSmargunif_2MLE, syntNBNSmargmarg_2MLE ,file="toyDataNS.RData" ) #syntNBNSunif,syntNBNSunif, yntNBNSmarg,syntNBNSmarg,
} else {load("toyDataNS.RData")}
```

```{r Plot results no signal, purl=FALSE, include=FALSE}
solListNS = list("unifunif" = syntNBNSunifunif_1B1_3, "unifmarg" = syntNBNSunifmarg_1B1_3, "margunif" = syntNBNSmargunif_1B1_3, "margmarg" = syntNBNSmargmarg_1B1_3, "unifunifMLE" = syntNBNSunifunif_1B1_2MLE, "unifmargMLE" = syntNBNSunifmarg_1B1_2MLE, "margunifMLE" = syntNBNSmargunif_1B1_2MLE, "margmargMLE" = syntNBNSmargmarg_1B1_2MLE)# "unifmargunifRmarg" = test1,  "margmargunifRmarg"=test2, "margunifunifRmarg" = test3, "unifmargunifRunif" = test1unif, "margmargunifRunif"=test2unif, "margunifunifRunif" = test3unif,
#Runtimes and convergence
sapply(solListNS,function(x){x$converged})
```

```{r SamplePlotsNoSignal, results='hide', purl=FALSE}
par(mfrow=c(2,4))
#unifmarg
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {
  rMatPsi = rMat %*% diag(psis)
  dfCol = data.frame(Dim1=rMatPsi[,1], Dim2=rMatPsi[,2], col=rowSums(X))
  ggplot(data=dfCol, aes(x=Dim1, y=Dim2, col=col)) +geom_point(size=2) +ggtitle(Y) + scale_colour_continuous(name="Library sizes", low="red",high = "green")
})
})
```

```{r loadings vs. library sizes dim1, results='hide', purl=FALSE}
#Dimension 1
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y,rMat[,1] *psis[1] ,rowSums(X), ylab="Library sizes",xlab="Dim1", sub = paste0("Cor = ",round(cor(rMat[,1], rowSums(X)),2)) )})})
```

```{r loadings vs. library sizes dim2, results='hide', purl = FALSE}
#Dimension 2
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y,rMat[,2] *psis[2],rowSums(X), ylab="Library sizes",xlab="Dim2", sub = paste0("Cor = ",round(cor(rMat[,2], rowSums(X)),2)) )})})
```

```{r loadings vs. library sizes dim3, results='hide', purl = FALSE}
#Dimension 3
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){try(with(solListNS[[Y]], {plot(main=Y,rMat[,3] *psis[3],rowSums(X), ylab="Library sizes",xlab="Dim3", sub = paste0("Cor = ",round(cor(rMat[,3], rowSums(X)),2)))}))})
```

In all weighing schemes, the row scores of the first dimension increase with library size (not in absolute values but they form a gradient). This is likely a result of a discrepancy between the MLE estimates of the sequencing depths under the NB model and the library size (row sums of our count table $\sum_{j=1}^px_{ij}$, which are also just an _estimate_ of the sequencing depth, which does not account for different variances between species but treats them all equally).  The sum of NB-distributed variables is also approximately NB-distributed, although this is not mathematically guaranteed. The first dimension of the MLE row scores estimates then tries to correct for this aberration so that its effect disappears in the higher dimensions. We plot the estimates of the sequencing depths of the two different estimators

```{r plotSequencingDepthEstimators, purl = FALSE}
par(mfrow=c(1,1), pty="s")
with(solListNS[["unifunifMLE"]], {plot(rowSums(X),libSizesMLE, xlab ="Library sizes", ylab = "MLE sequencing depth",log="xy")})
abline(0,1);abline(h=c(1e4,1e5), col="red");abline(v=c(1e4,1e5), col="red")
MSElibs = mean((rowSums(solListNS[["unifunifMLE"]]$X)-rep(c(1e4,1e5), each=150))^2)
MSEmle = mean((solListNS[["unifunifMLE"]]$libSizesMLE-rep(c(1e4,1e5), each=150))^2)
```

There exist considerable differences between both estimation strategies, the MLE one seems less variable. Still the library sizes lead to the smallest MSE.

###### Taxon plots

Now let's take a look at the taxon plots

```{r NBNS taxon plots, results='hide', purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y,t(cMat), ylab="Dim2",xlab="Dim1")}) })
```

The taxon plots of the fits with uniform weights for the taxa are used ($z_j = 1/p$) suffer from outliers when the abundances are estimated as column sums. For the other weighting schemes it is hard to establish anay effect for now. Contrarily to my expectations it are the weights for the centering rather than the normalization requirement that cause the outliers. With MLE estimation of the abundances, the outliers are gone in all weighting schemes!

The taxon plot of second and third dimension

```{r taxonPlot2and3D, results='hide', purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){try(with(solListNS[[Y]], {plot(main=Y,t(cMat[2:3,]), xlab = "Dim2",ylab="Dim3")}))})
```

No more outliers in the higher dimensions

Are the taxon scores related to the abundances in any way?

```{r taxonScesvsAbunds1D, results='hide', purl = FALSE}
#Abundances
par(mfrow=c(2,4))
sapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y,abs(cMat[1,]), colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="abs(Dim1)", sub = paste0("Cor = ",round(cor(abs(cMat[1,]), log(colSums(X))),2)))})})
```

Marginal centering weights causes the column scores to slightly correlate with the abundances

```{r taxonScesvsAbunds23D, results='hide', purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y,abs(cMat[2,]), colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="abs(Dim2)", sub = paste0("Cor = ",round(cor(abs(cMat[2,]), log(colSums(X))),2)))})})
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){try(with(solListNS[[Y]], {plot(main=Y,cMat[3,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="abs(Dim3)", sub = paste0("Cor = ",round(cor(abs(cMat[3,]), log(colSums(X))),2)))}))})
```

No more problems in higher dimensions.

Do the taxon scores relate to the overdispersions?

```{r taxonScesvsthetas, results='hide', purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y, abs(cMat[1,]), thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="abs(Dim1)", sub = paste0("Cor = ",round(cor(abs(cMat[1,]), log(thetas)),2)))})})
par(mfrow=c(2,4))
lapply(names(solListNS), function(Y){with(solListNS[[Y]], {plot(main=Y, cMat[2,], thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="abs(Dim2)", sub = paste0("Cor = ",round(cor(abs(cMat[2,]), log(thetas)),2)))})})
```

Large overdispersions get larger scores in all dimensions. This need not be a problem, and may stand in relationship to abundances. Taxa ill fit by the NB may try to compensate this by a large OD and more extreme scores?

How do the MLE estimates for the mean abundances differ from the column sums?

```{r MLEvsColsums, results='hide', purl = FALSE}
par(mfrow=c(1,1))
with(solListNS[["unifunifMLE"]], {plot(main="MLE abundances vs column sum abundances", colSums(X)/sum(X), abundsMLE, log="xy", cex=0.5, ylab="MLE abundances", xlab="Column sum abundance")})
abline(0,1)
```

There are subtle differences but less than with the library sizes. Why is the MLE less different from the marginal sum in this case? We can think of the following mathematical explanation: when calculating library sizes or abundances by row and column sums, each observations receives the same weight. However, when we estimate the margins through ML, we solve the following score equations:

$$\frac{\partial L(\mathbf{X}_i|\mathbf{u}_{i}, \mathbf{u}_j, \boldsymbol{\theta})}{\partial u_{i}} = \sum_{j=1}^p \frac{y_{ij}-\mu_{ij}}{1+\frac{\mu_{ij}}{\theta_j}}$$

$$\frac{\partial L(\mathbf{X}_i|\mathbf{u}_{i}, \mathbf{u}_j, \boldsymbol{\theta})}{\partial u_{j}} = \sum_{i=1}^n \frac{y_{ij}-\mu_{ij}}{1+\frac{\mu_{ij}}{\theta_j}}.$$

Note that since we are estimating offsets the value of the regressor is 1. In this case the difference of $y_{ij}$ with the expected value $\mu_{ij}$ is weighted by a factor $\frac{1}{1+\frac{\mu_{ij}}{\theta_j}} = \frac{\theta_j}{\theta_j+\mu_{ij}}$. When estimating $u_j$, $\theta$ is a constant but when estimating $u_i$ it is different for every observation $y_{ij}$. Hence the weights put on every observation differ much more when estimating the sample offsets than when estimating the taxon offsets. That is the reason why the MLE differs more from the marginal sum for the library size than for the abundances.

##### Influence function

Take a look at the influence function values for the influence on the first psi parameters

```{r NS influence measures: psis, results="hide", purl = FALSE, eval = FALSE}
infl1 = lapply(solListNS, function(x, i){
infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
id=abs(infl) > quantile( abs(infl),0.995)
Xid = x$X[id]
abunds = colSums(x$X)/sum(x$X)
libsizes = rowSums(x$X)
list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},1)
# par(mfrow=c(3,5))
# lapply(names(infl1), function(x){
#   with(infl1[[x]], plot(abunds, colSums(id), log="x", main=x, ylab= "Number of influential observations"))
# })
par(mfrow=c(2,4))
lapply(names(infl1), function(x){
with(infl1[[x]], plot(log(abs(c(t(infl))))~rep(cut(log(abunds),8),nrow(infl)), main=x, ylab = "log(|Influence|)",xlab = "log Abundance")) #Very heavy!
})
```

Larger abundances, larger influences in all weighting schemes (remember it is a log(abs())-plot).

Influence vs libSizes

```{r Influence vs libSizes, results="hide", purl = FALSE, eval = FALSE}
par(mfrow=c(2,4))
lapply(names(infl1), function(x){
with(infl1[[x]], plot(log(abs(c(infl)))~rep(cut(log(libsizes),8),ncol(infl)), main=x, ylab = "log(|Influence|)",xlab = "log library size")) #Very heavy!
})
```

No clear trend of influence in function of the library sizes

How are the most influential observations as compared to the bulk of the observations?

```{r Influence vs non-influential observations, purl = FALSE, eval = FALSE}
lapply(names(infl1), function(x){
with(infl1[[x]], {rbind(expectationsAll = quantile(outer(libsizes, abunds)),
expectationsInfluentialObservation =quantile(outer(libsizes, abunds)[id])) })
})

# sapply(names(infl1), function(x){
#   rbind(zeroesAll = mean(solListNS[[x]]$X==0),
#        zeroesInfluentialObservation = mean(solListNS[[x]]$X[infl1[[x]]$Xid]==0)) })
```

The most influential observations are zero counts (or low counts) in highly abundant species and high libsizes (i.e. high expected values)!

Plot the influence vs. observed abundances

```{r X vs Infl(X), results="hide", purl = FALSE, eval = FALSE}
#Expectations
par(mfrow=c(2,4))
lapply(names(infl1), function(x){
id = sample(1:prod(dim(solListNS[[x]]$X)), size=10000)
plot(log10(solListNS[[x]]$X+1)[id]setwd(WD), infl1[[x]]$infl[id], ylab= "Influence", xlab = "log10(Counts +1)", main=x)
abline(v=0, col="red")
})
```

Versus expected abundances

```{r Infl vs observed abundances, results="hide", purl = FALSE, eval = FALSE}
par(mfrow=c(2,4))
lapply(names(infl1), function(x){
id = sample(1:prod(dim(solListNS[[x]]$X)), size=10000)
plot(outer(rowSums(solListNS[[x]]$X), colSums(solListNS[[x]]$X)/sum(solListNS[[x]]$X))[id], infl1[[x]]$infl[id], ylab = "Influence", xlab = "Expected counts", main=x)
abline(v=0, col="red")
})
```

Versus Raw residuals (Observed - Expected count)

```{r vs Pearson residuals, results = "hide", purl = FALSE, eval = FALSE}
par(mfrow=c(2,4))
lapply(names(infl1), function(x){
id = sample(1:prod(dim(solListNS[[x]]$X)), size=10000) #Take a susbet for computability
plot((solListNS[[x]]$X-outer(rowSums(solListNS[[x]]$X), colSums(solListNS[[x]]$X)/sum(solListNS[[x]]$X)))[id], infl1[[x]]$infl[id], ylab = "Influence", xlab = "Raw residuals", main=x)#sqrt(diag(1/rowSums(solListNS[[x]]$X)))%*%%*%sqrt(diag(1/colSums(solListNS[[x]]$X)))
abline(v=0, col="red")
})
```


```{r psi2Infl, eval=FALSE, purl = FALSE, eval = FALSE}
infl2 = lapply(solListNS, function(x, i){
infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
id=abs(infl) > quantile( abs(infl),0.995)
Xid = x$X[id]
abunds = colSums(x$X)/sum(x$X)
libsizes = rowSums(x$X)
list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},2)
par(mfrow=c(2,2))
lapply(names(infl2), function(x){
with(infl2[[x]], plot(abunds, colSums(id), log="x", main=x))
})
#In the second dimension, larger abundance means larger influence on the psis too, and no more outliers
lapply(names(infl2), function(x){
with(infl2[[x]], plot(rep(abunds, nrow(infl)), c(t(infl)), log="x", main=x)) #Very heavy!
}) #Outlier in unifmarg weigthing scheme
lapply(names(infl2), function(x){
with(infl2[[x]], plot(rep(libsizes, ncol(infl)), c(infl), log="x", main=x)) #Very heavy!
})
#Expectations
lapply(names(infl2), function(x){
with(infl2[[x]], {rbind(quantile(outer(libsizes, abunds)),
quantile(outer(libsizes, abunds)[id])) })
})
lapply(names(infl2), function(x){
cat(x, "\n")
with(infl2[[x]], table(Xid))
})
```

Under the uniform weighting scheme for the taxa, one single observation will get an enormous influence in the first dimension, this is in all likelihood the outlier. Note that we cen derive this from the influence function that does not even depend on the weights!

Influence on the colScores

```{r NBNS influence colscores, eval=FALSE, purl = FALSE, eval = FALSE}
inflCol1 = lapply(solListNS, function(x){
try(with(x, NBcolInfl(X, psis, cMat, rMat,thetas , colWeights , k=1 , lambdaCol)))
})
#Only look at most extreme colscores
#Cannot calculate influence functions for the interesting cases (with the outliers)

```

```{r NBNS influence rowscores, eval=FALSE, purl = FALSE, eval = FALSE}
inflRowList1 = lapply(solListNS, function(x){
try(with(x, NBrowInfl(X, psis, cMat, rMat,thetas , rowWeights , k=1 , lambdaRow)))
})
#Only look at most extreme colscores
#Cannot calculate influence functions for the interesting cases (with the outliers)
inflRowList1=inflRowList1[sapply(inflRowList1, class)=="list"]
#The influence on the first row score
inflRow1 = lapply(inflRowList1, function(x){getInflRow(x$score, x$InvJac, 1)})
par(mfrow=c(1,2))
libSizes = rowSums(solListNS[[1]]$X)
lapply(inflRow1, function(x){
plot(y=x[,1], libSizes, log="x")
})
```

##### Conclusions from the null plots

- No signal in either plot.
- For the uniform taxon weighting the first dimension of the plot is dominated by outliers
- The first dimension row scores are related to the library sizes in all weighting schemes
- The least abundant species get higher column scores in all weighting schemes

#### NB with signal

The signal is introduced by sampling from different distributions in which each time the abundance of a small set of species was modified.

With modified abundances (4 groups)

```{r NB with signal abundance,eval=FALSE , purl = FALSE}
rhosNBSigref = rhosNBSig1 = rhosNBSig2 = rhosNBSig3 = rhosNBSig4 = rhosNBNS

load("/home/stijn/PhD/American Gut/AGphylo.RData")
libSizesAG = sample_sums(AGphylo)

NsamplesSignal1 = NsamplesSignal2 =  20
NsamplesSignal3 = NsamplesSignal4 =  15

NtaxaSignalNB = 40

Signal1NB = 18
Signal2NB = 12
Signal3NB = 10
Signal4NB = 7

idSig1NB = 1:NtaxaNBNS %in% sample(1:NtaxaNBNS, NtaxaSignalNB) #Random sampling should ensure orthogonality
idSig2NB = 1:NtaxaNBNS %in% sample(1:NtaxaNBNS, NtaxaSignalNB)
idSig3NB = 1:NtaxaNBNS %in% sample(1:NtaxaNBNS, NtaxaSignalNB) #Random sampling should ensure orthogonality
idSig4NB = 1:NtaxaNBNS %in% sample(1:NtaxaNBNS, NtaxaSignalNB)
#Apply the signals

rhosNBSig1[idSig1NB] = rhosNBSig1[idSig1NB]*Signal1NB
rhosNBSig2[idSig2NB] = rhosNBSig2[idSig2NB]*Signal2NB
rhosNBSig3[idSig3NB] = rhosNBSig3[idSig3NB]*Signal3NB
rhosNBSig4[idSig4NB] = rhosNBSig4[idSig4NB]*Signal4NB

#Renormalize
renorm=function(x){x/sum(x)}
rhosNBSig1=renorm(rhosNBSig1);rhosNBSig2=renorm(rhosNBSig2);
rhosNBSig3=renorm(rhosNBSig3);rhosNBSig4=renorm(rhosNBSig4);

#Generate data
Nref = (NsamplesNBNS-NsamplesSignal1-NsamplesSignal2-NsamplesSignal3-NsamplesSignal4)
meanMatRefNB = outer(libSizesAG[sample(size=Nref, 1:NsamplesNBNS)], rhosNBSigref)
meanMatSig1NB = outer(libSizesAG[sample(size=NsamplesSignal1, 1:NsamplesNBNS)], rhosNBSig1)
meanMatSig2NB = outer(libSizesAG[sample(size=NsamplesSignal2, 1:NsamplesNBNS)], rhosNBSig2)
meanMatSig3NB = outer(libSizesAG[sample(size=NsamplesSignal3, 1:NsamplesNBNS)], rhosNBSig3)
meanMatSig4NB = outer(libSizesAG[sample(size=NsamplesSignal4, 1:NsamplesNBNS)], rhosNBSig4)
thetaMatSigNB = thetaMatNBNS

dataMatRefNB = makeNBdata(meanMatRefNB, thetaMatSigNB)
dataMatSig1NB = makeNBdata(meanMatSig1NB, thetaMatSigNB)
dataMatSig2NB = makeNBdata(meanMatSig2NB, thetaMatSigNB)
dataMatSig3NB = makeNBdata(meanMatSig3NB, thetaMatSigNB)
dataMatSig4NB = makeNBdata(meanMatSig4NB, thetaMatSigNB)

dataMatSigNBab = rbind(dataMatRefNB, dataMatSig1NB, dataMatSig2NB, dataMatSig3NB,dataMatSig4NB)

#Save signals
sampleSigNBab = factor(c(rep("Reference", Nref), rep("Signal1", NsamplesSignal1), rep("Signal2", NsamplesSignal2), rep("Signal 3", NsamplesSignal3), rep("Signal 4", NsamplesSignal4)))
taxaSigNBtmp = rep("Reference",ncol(dataMatSigNBab))
taxaSigNBtmp[idSig1NB] = "Signal 1"
taxaSigNBtmp[idSig2NB] = "Signal 2"
taxaSigNBtmp[idSig3NB] = "Signal 3"
taxaSigNBtmp[idSig4NB] = "Signal 4"
taxaSigNBtmp[idSig1NB & idSig2NB] = "Signal 1 and 2"
taxaSigNBtmp[idSig1NB & idSig3NB] = "Signal 1 and 3"
taxaSigNBtmp[idSig1NB & idSig4NB] = "Signal 1 and 4"
taxaSigNBtmp[idSig2NB & idSig3NB] = "Signal 2 and 3"
taxaSigNBtmp[idSig2NB & idSig4NB] = "Signal 2 and 4"
taxaSigNBtmp[idSig3NB & idSig4NB] = "Signal 3 and 4"
taxaSigNBtmp[idSig1NB & idSig2NB & idSig3NB] = "Signal 1 and 2 and 3"
taxaSigNBtmp[idSig4NB & idSig2NB & idSig3NB] = "Signal 2 and 3 and 4"
taxaSigNBtmp[idSig1NB & idSig2NB & idSig4NB] = "Signal 1 and 2 and 4"
taxaSigNBtmp[idSig1NB & idSig4NB & idSig3NB] = "Signal 1 and 3 and 4"
taxaSigNBab = factor(taxaSigNBtmp, levels = c("Reference","Signal 1","Signal 2","Signal 3","Signal 4","Signal 1 and 2","Signal 1 and 3","Signal 1 and 4","Signal 2 and 3","Signal 2 and 4","Signal 3 and 4", "Signal 1 and 2 and 3", "Signal 1 and 2 and 4", "Signal 1 and 3 and 4", "Signal 2 and 3 and 4"))
names(taxaSigNBab) = colnames(dataMatSigNBab) = names(rhosNBSigref)
names(sampleSigNBab) = rownames(dataMatSigNBab) = 1:NsamplesNBNS

save(dataMatSigNBab, taxaSigNBab, sampleSigNBab, file = "abData.RData")
```

#### The abundance based approach: fit

```{r Fit3Dab, eval=FALSE, purl = FALSE}
nleqslv.control.ab = list(trace=TRUE, maxit = 500, cndtol=.Machine$double.eps)
if(!file.exists("toyDataSigab.RData")){

syntNBSigmargmarg_abJob = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "marginal"))
syntNBSigmargmarg_ab = mccollect(syntNBSigmargmarg_abJob, FALSE)[[1]] #Check

syntNBSigunifmarg_abJob = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "uniform"))
syntNBSigunifmarg_ab = mccollect(syntNBSigunifmarg_abJob, FALSE)[[1]]

syntNBSigmargunif_abJob = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, rowWeights="marginal", colWeights = "uniform"))
syntNBSigmargunif_ab = mccollect(syntNBSigmargunif_abJob, FALSE)[[1]]

syntNBSigunifunif_abJob = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="uniform", rowWeights = "uniform"))
syntNBSigunifunif_ab = mccollect(syntNBSigunifunif_abJob, FALSE)[[1]] #Check

syntNBSigmargmarg_abJobMLE = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "marginal", marginEst = "MLE"))
syntNBSigmargmarg_abMLE = mccollect(syntNBSigmargmarg_abJobMLE, FALSE)[[1]]

syntNBSigunifmarg_abJobMLE = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "uniform", marginEst = "MLE"))
syntNBSigunifmarg_abMLE = mccollect(syntNBSigunifmarg_abJobMLE, FALSE)[[1]]

syntNBSigmargunif_abJobMLE = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, rowWeights="marginal", colWeights = "uniform", marginEst = "MLE"))
syntNBSigmargunif_abMLE = mccollect(syntNBSigmargunif_abJobMLE, FALSE)[[1]]

syntNBSigunifunif_abJobMLE = mcparallel(RCM(dataMatSigNBab, distribution="NB", k=3, nleqslv.control= nleqslv.control.ab, maxItOut=2e3, prevCutOff=0.01, colWeights="uniform", rowWeights = "uniform", marginEst = "MLE"))
syntNBSigunifunif_abMLE = mccollect(syntNBSigunifunif_abJobMLE, FALSE)[[1]]

save(dataMatSigNBab, syntNBSigunifmarg_ab, syntNBSigunifunif_ab,syntNBSigmargmarg_ab,syntNBSigmargunif_ab,taxaSigNBab, sampleSigNBab, testAb1, testAb2, testAb3, testAb1unif, testAb2unif, testAb3unif,testAbInv,testAbLibs,dataMatSigNBabag,syntNBSigmargmarg_abMLE,syntNBSigunifunif_abMLE,syntNBSigmargunif_abMLE,syntNBSigunifmarg_abMLE, file="toyDataSigab.RData")#syntNBSigunif, syntNBSigmarg, syntNBSigmarg_3,
} else {load("toyDataSigab.RData")}
```

```{r ab plots, include=FALSE, purl = FALSE}
load("toyDataSigABserver.RData")
# solListWSab = list("unifunif" = syntNBSigunifunif_ab, "unifmarg" = syntNBSigunifmarg_ab, "margunif" = syntNBSigmargunif_ab, "margmarg" = syntNBSigmargmarg_ab, "unifmargunifRmarg" = testAb1,  "margmargunifRmarg"=testAb2, "margunifunifRmarg" = testAb3, "unifmargunifRunif" = testAb1unif, "margmargunifRunif"=testAb2unif, "margunifunifRunif" = testAb3unif, "inverse"= testAbInv, "testEqualLibs" =   testAbLibs, "unifunifMLE" = syntNBSigunifunif_abMLE, "unifmargMLE" = syntNBSigunifmarg_abMLE, "margunifMLE" = syntNBSigmargunif_abMLE, "margmargMLE" = syntNBSigmargmarg_abMLE )
names(sigListAB) = c("unifunif","unifunifMLE","margunif","margunifMLE","unifmarg","unifmargMLE","margmarg","margmargMLE")
solListWSab = sigListAB

solListWSab = solListWSab[sapply(solListWSab, class)=="list"]
#Runtimes and convergence
sapply(solListWSab,function(x){x$converged})
sapply(solListWSab,function(x){x$runtime})
sapply(solListWSab,function(x){x$iter})
```

The first abbreviation refers to the weighting scheme for the rows(samples), the second to the weighting scheme for the columns. E.g "margunif" means $w_i = x_{i.}$ and $z_j = 1/p$. The MLE epitheton indicates that the independence model was estimated by MLE.

#### Sample plots

\newpage

```{r abSignalPlots, results='hide', fig.show="hold", fig.width = 12, fig.height = 9, purl = FALSE, eval=FALSE}
cols = c("grey", "red","blue","purple","green","brown","cyan","orange","black",  "magenta","yellow","pink", "darkblue","grey75","cadetblue")
palette(cols)
par(mfrow=c(1,1), mar = c(4,4,4,4), mai=c(0,0,0,0))
lapply(names(solListWSab), function(Y){plotRCM(solListWSab[[Y]], samColour = sampleSigNBab, main=Y, biplot=FALSE, Dim = c(1,2), libInset=c(-0.1,0))})#, libLeg = Y ==names(solListWSab)[1]
```

##### 1st and 3rd dimensions

\newpage

```{r abSignalPlots23D, results='hide', fig.show="hold", fig.width = 12, fig.height = 9, purl = FALSE}
palette(cols)
par(mfrow=c(1,1))
lapply(names(solListWSab), function(Y){try(plotRCM(solListWSab[[Y]], samColour = sampleSigNBab, main=Y, biplot=FALSE, Dim = c(1,3), libInset=c(-0.8,0)))})
```

Only when the margins are estimated with MLE are the first two dimensions used to separate the signal well. Especially the uniform weighting for the samples and marginal for the columns (as was our first intuition) appears to perform well. Also uniform weighting for both yields good results.

#### Sample scores vs. libsizes

```{r abSignalPlots 13D, results='hide', fig.show="hold", purl = FALSE}
par(mfrow=c(2,4))
#Maybe a function of the library sizes
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {
rMatPsi = rMat %*% diag(psis)
dfCol = data.frame(Dim1=rMatPsi[,1], Dim2=rMatPsi[,2], col=log(rowSums(X)))
ggplot(data=dfCol, aes(x=Dim1, y=Dim2, col=col)) +geom_point(size=3) +ggtitle(Y) + scale_colour_continuous(name="Library sizes", low="red",high = "green")
})
})
```

Look at the scores' relationship to the library sizes more directly.

```{r ab libsizes, results='hide', fig.show="hold", eval=FALSE, purl = FALSE}
#Look at the loadings in function of the library sizes and abunds
par(pty = "m", mfrow = c(1,1), mar = c(5,4,4,4))
logP="x"
#Libsizes
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,rMat[,1] *psis[1],x=rowSums(X), log=logP, xlab ="Dim1",ylab = "Library sizes", sub = paste0("Cor = ",round(cor(rMat[,1] ,rowSums(X)),2) ) )})})
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,rMat[,2] *psis[2],x=rowSums(X), log=logP, xlab ="Dim2", ylab = "Library sizes", sub = paste0("Cor = ",round(cor(rMat[,2] ,rowSums(X)),2) ) )})})
# lapply(names(solListWSab), function(Y){try(with(solListWSab[[Y]], {plot(main=Y,rMat[,3] *psis[3],x=rowSums(X), log=logP, xlab ="Dim3", ylab = "Library sizes", sub = paste0("Cor = ",round(cor(rMat[,1] ,rowSums(X)),2) ) )}))})
```

On these plots it is clear that when the margins are estimated using the library sizes, the scores correlate with the library sizes.

In the MLE framework for the offsets, dependence on the library sizes has disappeared.

Compare MLE and rowSums estimators for the sequencing depth

```{r plotSequencingDepthEstimatorsAb, purl = FALSE}
par(mfrow=c(1,1), pty="s")
with(solListWSab[["unifunifMLE"]], {plot(rowSums(X),libSizesMLE, xlab ="Library sizes", ylab = "MLE sequencing depth",log="xy")})
abline(0,1);#abline(h=c(1e4,1e5), col="red");abline(v=c(1e4,1e5), col="red")
# MSElibs = mean((rowSums(solListWSab[["unifunifMLE"]]$X)-rep(c(1e4,1e5), each=150))^2)
# MSEmle = mean((solListWSab[["unifunifMLE"]]$libSizesMLE-rep(c(1e4,1e5), each=150))^2)
```

The MLE estimate is generally larger than the rowSums estimate (presumably since the former puts less weight on zero counts (large overdispersion)), this may explain the gradient of the row scores in function of the library sizes.

##### Taxa plots

Now let's take a look at the taxon plots

```{r NBab taxon plots, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(1,1), pty="s")
palette(cols)
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,t(cMat), ylab="Dim2",xlab="Dim1", col = taxaSigNBab)}); legend("topright",legend = levels(taxaSigNBab), col=cols,pch= 1, cex=0.7, inset = c(-0.80,0), xpd=TRUE) })
```

For the uniform weighting of the samples we see outliers for the taxa

Are the taxon scores related to the abundances in any way?

```{r taxonScesvsAbunds1Dab, results='hide', eval=FALSE, purl = FALSE}
#Abundances
par(mfrow=c(1,1), mar=c(5,4,4,4))
sapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,cMat[1,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim1", sub = paste0("Cor = ", round(cor(cMat[1,], colSums(X)),2)))})})
```

The MLE offset estimation does reduce correlation with the mean abundances.

```{r taxonScesvsAbunds23Dab, results='hide', eval=FALSE, purl = FALSE}
# Second and third dimensions
par(mfrow=c(1,1), mar=c(5,4,4,4))
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,cMat[2,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim2", sub = paste0("Cor = ", round(cor(cMat[2,], colSums(X)),2)))})})
par(mfrow=c(1,1), mar=c(5,4,4,4))
lapply(names(solListWSab), function(Y){try(with(solListWSab[[Y]], {plot(main=Y,cMat[3,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim3", sub = paste0("Cor = ", round(cor(cMat[3,], colSums(X)),2)))}))})
```

No more problems with correlation in higher dimensions.

Do the taxon scores relate to the overdispersions?

```{r taxonScesvsthetas1, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(1,1))
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,abs(cMat[1,]), thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="Dim1", sub = paste0("Cor = ", round(cor(abs(cMat[1,]), thetas),2)))})})
par(mfrow=c(1,1))
lapply(names(solListWSab), function(Y){with(solListWSab[[Y]], {plot(main=Y,cMat[2,], thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="Dim2", sub = paste0("Cor = ", round(cor(abs(cMat[2,]), thetas),2)))})})
```

No real relationship with overdispersions

```{r lagrange multipliers, eval=FALSE, purl = FALSE}
##### Lagrange multipliers
par(mfrow=c(2,4))
lapply(names(solListWSab), function(Y){cat(Y, "\n");solListWSab[[Y]]$lambdaCol})
```

##### Biplots

First two dimensions

```{r ab biplots, results='hide', fig.show="hold", fig.width = 12, fig.height = 9, purl = FALSE}
#par(mfrow=c(2,4), mai = c(2,2,2,4))
par(mfrow=c(1,1))
palette(cols)
lapply(names(solListWSab), function(Y){plotRCM(solListWSab[[Y]], samColour = sampleSigNBab, main=Y, biplot=TRUE, Dim = c(1,2), taxColour=taxaSigNBab, taxLegPos = "bottomright", abundLeg=TRUE, taxCol  =cols, taxInset = c(-0.0,-0.0), libCex = 0.65, libInset = c(0.0,-0.0), arrowFrac = 0.05)})
```

Fisrt and third dimension

```{r ab23 biplots, results='hide', fig.show="hold", fig.width = 12, fig.height = 9, purl = FALSE}
#par(mfrow=c(2,4), mai = c(2,2,2,4))
par(mfrow=c(1,1))
palette(cols)
lapply(names(solListWSab), function(Y){plotRCM(solListWSab[[Y]], samColour = sampleSigNBab, main=Y, biplot=TRUE, Dim = c(1,3), taxColour=taxaSigNBab, taxLegPos = "bottomright", abundLeg=TRUE, taxCol  =cols, taxInset = c(-0,-0), libCex = 0.65, libInset = c(-0.0,-0.0), arrowFrac = 0.05)})
```

##### Influence function

Take a look at the influence function values on the psis

```{r Ab influence measures: psis, eval=FALSE, purl = FALSE}
infl1 = lapply(solListWSab, function(x, i){
infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
id=abs(infl) > quantile( abs(infl),0.995)
Xid = x$X[id]
abunds = colSums(x$X)/sum(x$X)
libsizes = rowSums(x$X)
list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},1)
par(mfrow=c(1,2))
lapply(names(infl1), function(x){
with(infl1[[x]], plot(abunds, colSums(id), log="x", main=x, ylab="Number of very influential observations"))
}) #High abundances, larger influence
lapply(names(infl1), function(x){
idTmp = sample(seq_along(infl1[[x]]$abunds), 200) #For speed
with(infl1[[x]], plot(rep(abunds[idTmp], nrow(infl)), c(t(infl[,idTmp])), log="x", main=x, xlab="Abundance",ylab = "Influence")) #Very heavy!
})
lapply(names(infl1), function(x){
idTmp = sample(seq_along(infl1[[x]]$libsizes), 200) #For speed
with(infl1[[x]], plot(rep(libsizes[idTmp], ncol(infl)), c(infl[idTmp,]), log="", main=x, xlab ="Library sizes",ylab = "Influence")) #Very heavy!
})
#Expectations
lapply(names(infl1), function(x){
with(infl1[[x]], {rbind(quantile(outer(libsizes, abunds)),
quantile(outer(libsizes, abunds)[id])) })
}) #Influential observations have high expectations
lapply(names(infl1), function(x){
cat(x, "\n")
with(infl1[[x]], table(Xid))
}) # ... but are very often zero!

infl2 = lapply(solListWSab, function(x, i){
infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
id=abs(infl) > quantile( abs(infl),0.995)
Xid = x$X[id]
abunds = colSums(x$X)/sum(x$X)
libsizes = rowSums(x$X)
list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},2)
par(mfrow=c(1,2))
lapply(names(infl2), function(x){
with(infl2[[x]], plot(abunds, colSums(id), log="x", main=x))
})
#In the second dimension, larger abundance means larger influence on the psis too, and no outliers
lapply(names(infl2), function(x){
with(infl2[[x]], plot(rep(abunds, nrow(infl)), c(t(infl)), log="x", main=x, xlab="Abundance",ylab = "Influence")) #Very heavy!
}) #Outlier in unifmarg weigthing scheme
lapply(names(infl2), function(x){
with(infl2[[x]], plot(rep(libsizes, ncol(infl)), c(infl), log="x", main=x, xlab="Library size",ylab = "Influence")) #Very heavy! #No trend in terms of library sizes visible from here, but their distribution should be log-normal perhaps
})
#Expectations
lapply(names(infl2), function(x){
with(infl2[[x]], {rbind(quantile(outer(libsizes, abunds)),
quantile(outer(libsizes, abunds)[id])) })
})
lapply(names(infl2), function(x){
cat(x, "\n")
with(infl2[[x]], table(Xid))
})
```

The most influential observations are zero counts in highly abundant species and high libsizes (i.e. high expectations)! Note that we cen derive this from the influence function that does not even depend on the weights!

Influence on the colScores

```{r Ab influence colscores, eval=FALSE, purl = FALSE}
inflColList1 = lapply(solListWSab, function(x){
try(with(x, NBcolInfl(X, psis, cMat, rMat,thetas , colWeights = if(is.matrix(colWeights)) colWeights else{ cbind(colWeights, colWeights, colWeights)} , k=1 , lambdaCol)))
})
#Look at the signal from the first group
Id1 = which(taxaSigNBab=="Signal 1")[1]

inflCol1 = lapply(inflColList1, function(x){getInflCol(x$score, x$InvJac, Id1)})
lapply(inflCol1, function(x){
boxplot(x[,Id1]~sampleSigNBab, las=2, ylab="Influence")
})
```

As expected, the samples with the signal have the highest impact on the score

```{r Ab influence rowscores Ab, eval=FALSE, purl = FALSE}
inflRowList1 = lapply(solListWSab, function(x){
try(with(x, NBrowInfl(X, psis, cMat, rMat,thetas , rowWeights , k=1 , lambdaRow)))
})
#Only look at most extreme colscores
#Cannot calculate influence functions for the interesting cases (with the outliers)
inflRowList1=inflRowList1[sapply(inflRowList1, class)=="list"]
#The influence on the first row score

par(mfrow=c(1,2))
idRow = 10
libSizes = rowSums(solListWSab[[1]]$X)
inflRow1 = lapply(inflRowList1, function(x){getInflRow(x$score, x$InvJac, idRow)})
lapply(names(inflRow1), function(x){
plot(y=inflRow1[[x]][,idRow], libSizes, log="x", main=x, col=(libSizes ==libSizes[idRow])+1, ylab="Influence")
})
```

Evidently, the observation from the sample itself has the largest influence. In the marginal weighting scheme for the libsizes, the larger libsizes sometimes get a larger influence, although the effect may be small.


```{r tolerances, eval=FALSE, purl = FALSE}
#Achieved tolerances
tolAchPsi = with(syntNBSigmargmarg_ab, lapply(1:k, function(y){sapply( 2:iter[y], function(x){abs(1-psiRec[y,x]/psiRec[y,x-1])})}))
par(mfrow=c(1,3))
lapply(tolAchPsi, function(x){plot(x, log="y");abline(h=0.01, col="red");abline(h=0.001, col="blue")})

with(syntNBSigmargmarg_ab,plot(psiRec[3,1:iter[3]])) #Weird jumps, don't stop too early!
  ```

#### A log-uniform approach

We sample abundances and library sizes log-uniformly to get a better view of the scores in function of them

```{r NB with signal lu, purl = FALSE, eval = FALSE}
makeNBdata=function(meanMat, thetaMat){apply(array(data= c(meanMat, thetaMat), dim=c(nrow(meanMat), ncol(meanMat), 2)), c(1,2), function(x){rnbinom(1,mu=x[1], size=x[2])})}
NtaxaLogUnif = 1000
NsamplesLogUnif = 300
rhosRefLU = rhosSig1LU = rhosSig2LU = rhosSig3LU = rhosSig4LU = 10^(-runif(NtaxaLogUnif, 2.5,6))

NsamplesSignal1LU = NsamplesSignal2LU =  25
NsamplesSignal3LU = NsamplesSignal4LU =  20

NtaxaSignalNBLU = 40

Signal1NBLU = 10
Signal2NBLU = 8
Signal3NBLU = 7.5
Signal4NBLU = 7

idSig1NBLU = 1:NtaxaLogUnif %in% sample(1:NtaxaLogUnif, NtaxaSignalNBLU) #Random sampling should ensure orthogonality
idSig2NBLU = 1:NtaxaLogUnif %in% sample(1:NtaxaLogUnif, NtaxaSignalNBLU)
idSig3NBLU = 1:NtaxaLogUnif %in% sample(1:NtaxaLogUnif, NtaxaSignalNBLU) #Random sampling should ensure orthogonality
idSig4NBLU = 1:NtaxaLogUnif %in% sample(1:NtaxaLogUnif, NtaxaSignalNBLU)
#Apply the signals

rhosSig1LU[idSig1NBLU] = rhosSig1LU[idSig1NBLU]*Signal1NBLU
rhosSig2LU[idSig2NBLU] = rhosSig2LU[idSig2NBLU]*Signal2NBLU
rhosSig3LU[idSig3NBLU] = rhosSig3LU[idSig3NBLU]*Signal3NBLU
rhosSig4LU[idSig4NBLU] = rhosSig4LU[idSig4NBLU]*Signal4NBLU

#Renormalize
renorm=function(x){x/sum(x)}
rhosSig1LU=renorm(rhosSig1LU);rhosSig2LU=renorm(rhosSig2LU);
rhosSig3LU=renorm(rhosSig3LU);rhosSig4LU=renorm(rhosSig4LU);

libSizesNBLU = 10^runif(NsamplesLogUnif, 3.5,5.5)

#Generate data
Nref = (NsamplesLogUnif-NsamplesSignal1LU-NsamplesSignal2LU-NsamplesSignal3LU-NsamplesSignal4LU)
meanMatRefNBLU = outer(libSizesNBLU[sample(size=Nref, 1:NsamplesLogUnif)], rhosRefLU)
meanMatSig1NBLU = outer(libSizesNBLU[sample(size=NsamplesSignal1LU, 1:NsamplesLogUnif)], rhosSig1LU)
meanMatSig2NBLU = outer(libSizesNBLU[sample(size=NsamplesSignal2LU, 1:NsamplesLogUnif)], rhosSig2LU)
meanMatSig3NBLU = outer(libSizesNBLU[sample(size=NsamplesSignal3LU, 1:NsamplesLogUnif)], rhosSig3LU)
meanMatSig4NBLU = outer(libSizesNBLU[sample(size=NsamplesSignal4LU, 1:NsamplesLogUnif)], rhosSig4LU)
load("/home/stijn/PhD/American Gut/AGpars.RData")
thetaMatSigNBLU = matrix(sample(thetas, NtaxaLogUnif), byrow=TRUE, ncol=NtaxaLogUnif, nrow=NsamplesLogUnif)

dataMatRefNBLU = makeNBdata(meanMatRefNBLU, thetaMatSigNBLU)
dataMatSig1NBLU = makeNBdata(meanMatSig1NBLU, thetaMatSigNBLU)
dataMatSig2NBLU = makeNBdata(meanMatSig2NBLU, thetaMatSigNBLU)
dataMatSig3NBLU = makeNBdata(meanMatSig3NBLU, thetaMatSigNBLU)
dataMatSig4NBLU = makeNBdata(meanMatSig4NBLU, thetaMatSigNBLU)

dataMatSigNBLU = rbind(dataMatRefNBLU, dataMatSig1NBLU, dataMatSig2NBLU, dataMatSig3NBLU, dataMatSig4NBLU)

#Save signals
sampleSigNBLU = factor(c(rep("Reference", Nref), rep("Signal1", NsamplesSignal1LU), rep("Signal2", NsamplesSignal2LU), rep("Signal 3", NsamplesSignal3LU), rep("Signal 4", NsamplesSignal4LU)))
taxaSigNBLUtmp = rep("Reference",ncol(dataMatSigNBLU))
taxaSigNBLUtmp[idSig1NBLU] = "Signal 1"
taxaSigNBLUtmp[idSig2NBLU] = "Signal 2"
taxaSigNBLUtmp[idSig3NBLU] = "Signal 3"
taxaSigNBLUtmp[idSig4NBLU] = "Signal 4"
taxaSigNBLUtmp[idSig1NBLU & idSig2NBLU] = "Signal 1 and 2"
taxaSigNBLUtmp[idSig1NBLU & idSig3NBLU] = "Signal 1 and 3"
taxaSigNBLUtmp[idSig1NBLU & idSig4NBLU] = "Signal 1 and 4"
taxaSigNBLUtmp[idSig2NBLU & idSig3NBLU] = "Signal 2 and 3"
taxaSigNBLUtmp[idSig2NBLU & idSig4NBLU] = "Signal 2 and 4"
taxaSigNBLUtmp[idSig3NBLU & idSig4NBLU] = "Signal 3 and 4"
taxaSigNBLUtmp[idSig1NBLU & idSig2NBLU & idSig3NBLU] = "Signal 1 and 2 and 3"
taxaSigNBLUtmp[idSig4NBLU & idSig2NBLU & idSig3NBLU] = "Signal 2 and 3 and 4"
taxaSigNBLUtmp[idSig1NBLU & idSig2NBLU & idSig4NBLU] = "Signal 1 and 2 and 4"
taxaSigNBLUtmp[idSig1NBLU & idSig4NBLU & idSig3NBLU] = "Signal 1 and 3 and 4"
taxaSigNBLUtmp[idSig1NBLU & idSig4NBLU & idSig2NBLU & idSig3NBLU] = "Signal 1,2,3 and 4"
taxaSigNBLU = factor(taxaSigNBLUtmp)
names(taxaSigNBLU) = colnames(dataMatSigNBLU) = names(rhosRefLU)
names(sampleSigNBLU) = rownames(dataMatSigNBLU) = 1:NsamplesLogUnif
```

Fit the RC(M) model to the logUnif data

```{r Fit3Dlu, eval=FALSE, purl = FALSE}
nleqslv.control.lu = list(trace=TRUE, maxit = 500, cndtol=.Machine$double.eps)
if(!file.exists("toyDataSigLU.RData")){

  syntNBSigmargmarg_LUJob = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=2, nleqslv.control= nleqslv.control.lu, maxItOut=2e4, prevCutOff=0.01, colWeights="marginal", rowWeights = "marginal", marginEst="MLE"))
  syntNBSigmargmarg_LU = mccollect(syntNBSigmargmarg_LUJob, FALSE)[[1]]

  syntNBSigunifmarg_LUJob = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=2, nleqslv.control= nleqslv.control.lu, maxItOut=2e4, prevCutOff=0.01, colWeights="marginal", rowWeights = "uniform", marginEst="MLE"))
  syntNBSigunifmarg_LU = mccollect(syntNBSigunifmarg_LUJob, FALSE)[[1]]

  syntNBSigmargunif_LUJob = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=2, nleqslv.control= nleqslv.control.lu, maxItOut=2e4, prevCutOff=0.01, rowWeights="marginal", colWeights = "uniform", marginEst="MLE"))
  syntNBSigmargunif_LU = mccollect(syntNBSigmargunif_LUJob, FALSE)[[1]]

  syntNBSigunifunif_LUJob = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=2, nleqslv.control= nleqslv.control.lu, maxItOut=2e4, prevCutOff=0.01, colWeights="uniform", rowWeights = "uniform", marginEst="MLE"))
  syntNBSigunifunif_LU = mccollect(syntNBSigunifunif_LUJob, FALSE)[[1]]

  #     Dim = dim(dataMatSigNBLU)
  #   unifWeightsLU = rep(1/Dim[2], Dim[2])
  #   margWeightsLU = colSums(dataMatSigNBLU)/sum(dataMatSigNBLU)
  #   colWeightsTestLU1 = cbind(unifWeightsLU, margWeightsLU, unifWeightsLU)
  #   colWeightsTestLU2 = cbind(margWeightsLU, margWeightsLU, unifWeightsLU)
  #   colWeightsTestLU3 = cbind(margWeightsLU, unifWeightsLU, unifWeightsLU)

  syntNBSigmargmarg_LUJobMLE = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=3, nleqslv.control= nleqslv.control.lu, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "marginal", marginEst = "MLE"))
  syntNBSigmargmarg_LUMLE = mccollect(syntNBSigmargmarg_LUJobMLE, FALSE)[[1]]

  syntNBSigunifmarg_LUJobMLE = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=3, nleqslv.control= nleqslv.control.lu, maxItOut=2e3, prevCutOff=0.01, colWeights="marginal", rowWeights = "uniform", marginEst = "MLE"))
  syntNBSigunifmarg_LUMLE = mccollect(syntNBSigunifmarg_LUJobMLE, FALSE)[[1]]

  syntNBSigmargunif_LUJobMLE = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=3, nleqslv.control= nleqslv.control.lu, maxItOut=2e3, prevCutOff=0.01, rowWeights="marginal", colWeights = "uniform", marginEst = "MLE"))
  syntNBSigmargunif_LUMLE = mccollect(syntNBSigmargunif_LUJobMLE, FALSE)[[1]]

  syntNBSigunifunif_LUJobMLE = mcparallel(RCM(dataMatSigNBLU, distribution="NB", k=3, nleqslv.control= nleqslv.control.lu, maxItOut=2e3, prevCutOff=0.01, colWeights="uniform", rowWeights = "uniform", marginEst = "MLE"))
  syntNBSigunifunif_LUMLE = mccollect(syntNBSigunifunif_LUJobMLE, FALSE)[[1]]

  save(dataMatSigNBLU, syntNBSigunifmarg_LU, syntNBSigunifunif_LU,syntNBSigmargmarg_LU,syntNBSigmargunif_LU, sampleSigNBLU, taxaSigNBLU,  file="toyDataSigLU.RData")#syntNBSigunif, syntNBSigmarg, syntNBSigmarg_3,testLU1, testLU2, testLU3, testLU1unif, testLU2unif, testLU3unif,testLUInv,testLULibs,
} else {load("toyDataSigLU.RData")}
```

LU plots

```{r LU plots, include=FALSE, purl = FALSE}
load("toyDataSigLUserver.RData")
# solListWSLU = list("unifunif" = syntNBSigunifunif_LU, "unifmarg" = syntNBSigunifmarg_LU, "margunif" = syntNBSigmargunif_LU, "margmarg" = syntNBSigmargmarg_LU, "unifmargunifRmarg" = testLU1,  "margmargunifRmarg"=testLU2, "margunifunifRmarg" = testLU3, "unifmargunifRunif" = testLU1unif, "margmargunifRunif"=testLU2unif, "margunifunifRunif" = testLU3unif, "inverse"= testLUInv, "testEqualLibs" =   testLULibs, "unifunifMLE" = syntNBSigunifunif_LUMLE, "unifmargMLE" = syntNBSigunifmarg_LUMLE, "margunifMLE" = syntNBSigmargunif_LUMLE, "margmargMLE" = syntNBSigmargmarg_LUMLE )
names(sigListLU) = c("unifunif","unifunifMLE","margunif","margunifMLE","unifmarg","unifmargMLE","margmarg","margmargMLE")
solListWSLU = sigListLU

solListWSLU = solListWSLU[sapply(solListWSLU, class)=="list"]
#Runtimes and convergence
sapply(solListWSLU,function(x){x$converged})
sapply(solListWSLU,function(x){x$runtime})
sapply(solListWSLU,function(x){x$iter})
```

The first Abreviation refers to the weighting scheme for the rows(samples), the second to the weighting scheme for the columns. E.g "margunif" means $w_i = x_{i.}$ and $z_j = 1/p$. The MLE epitheton indicates that the independence model was estimated by MLE.

##### Sample plots

```{r LUSignalPlots, eval=FALSE, purl = FALSE}
cols = c("grey", "red","blue","purple","green","brown","cyan","black", "orange", "magenta","yellow","pink", "olive","grey75","grey85")
palette(cols)
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){plotRCM(solListWSLU[[Y]], samColour = sampleSigNBLU, main=Y, biplot=FALSE, Dim = c(1,2), libInset=c(-0.8,0))})
```

```{r LUSignalPlots23D, eval=FALSE, purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){try(plotRCM(solListWSLU[[Y]], samColour = sampleSigNBLU, main=Y, biplot=FALSE, Dim = c(1,3), libInset=c(-0.8,0)))})
```

Only when the margins are estimated with MLE are the first two dimensions used to separate the signal well. Especially the uniform weighting for the samples and marginal for the columns (as was our first intuition) appears to perform well. Also uniform weighting for both yields good results.

```{r LUSignalPlots 13D, eval=FALSE, purl = FALSE}
par(mfrow=c(2,4))
#Maybe a function of the library sizes
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {
  rMatPsi = rMat %*% diag(psis)
  dfCol = data.frame(Dim1=rMatPsi[,1], Dim2=rMatPsi[,2], col=log(rowSums(X)))
  ggplot(data=dfCol, aes(x=Dim1, y=Dim2, col=col)) + geom_point(size=3) + ggtitle(Y)+ scale_colour_continuous(name = "Library sizes", low="red",high = "green")
})
})
```

On these plots it is clear that when the margins are estimated using the library sizes, the scores correlate with the library sizes.

Look at the scores' relationship to the library sizes more directly.

```{r LU libsizes, eval=FALSE, results = "hide", purl = FALSE}
#Look at the loadings in function of the library sizes and abunds
par(pty = "m", mfrow = c(2,4), mar=c(5,4,4,4))
logP="y"
#Libsizes
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,rMat[,1] *psis[1],rowSums(X), log=logP, xlab ="Dim1", ylab = "Library sizes", sub = paste0("Cor = ", round(cor(rMat[,1], rowSums(X)),2)) )})})
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,rMat[,2] *psis[2],rowSums(X), log=logP, xlab ="Dim2", ylab = "Library sizes" , sub = paste0("Cor = ", round(cor(rMat[,2], rowSums(X)),2)))})})
lapply(names(solListWSLU), function(Y){try(with(solListWSLU[[Y]], {plot(main=Y,rMat[,3] *psis[3],rowSums(X), log=logP, xlab ="Dim3", ylab = "Library sizes" , sub = paste0("Cor = ", round(cor(rMat[,3], rowSums(X)),2)))}))})
```

In the MLE framework for the offsets, dependence on the library sizes has disappeared

##### Taxa plots

Now let's take a look at the taxon plots

```{r NBLU taxon plots, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,t(cMat), ylab="Dim2",xlab="Dim1", col = taxaSigNBLU)}); legend("topright",legend = levels(taxaSigNBLU), col=cols,pch= 1, cex=0.7, inset = c(-0.80,0), xpd=TRUE) })
```

For the uniform weighting of the samples we see outliers for the taxa

Are the taxon scores related to the abundances in any way?

```{r taxonScesvsabunds1DLU, results='hide', eval=FALSE, purl = FALSE}
#abundances
par(mfrow=c(2,4), mar=c(5,4,4,4))
sapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,cMat[1,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim1", sub = paste0("Cor = ", round(cor(cMat[1,], colSums(X)),2)))})})
```

The MLE offset estimation does reduce correlation with the mean abundances.

Second and third dimensions

```{r taxonScesvsabunds23DLU, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(2,4), mar=c(5,4,4,4))
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,cMat[2,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim2", sub = paste0("Cor = ", round(cor(cMat[2,], colSums(X)),2)))})})
par(mfrow=c(2,4), mar=c(5,4,4,4))
lapply(names(solListWSLU), function(Y){try(with(solListWSLU[[Y]], {plot(main=Y,cMat[3,], colSums(X), log="y", cex=0.5, ylab="Mean  abundances",xlab="Dim3", sub = paste0("Cor = ", round(cor(cMat[3,], colSums(X)),2)))}))})
```

No more problems with correlation in higher dimensions.

Do the taxon scores relate to the overdispersions?

```{r taxonScesvsthetas2, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,abs(cMat[1,]), thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="Dim1", sub = paste0("Cor = ", round(cor(abs(cMat[1,]), thetas),2)))})})
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){with(solListWSLU[[Y]], {plot(main=Y,cMat[2,], thetas, log="y", cex=0.5, ylab="Overdispersion", xlab="Dim2", sub = paste0("Cor = ", round(cor(abs(cMat[2,]), thetas),2)))})})
```

No real relationship with overdispersions

```{r lagrange multiplierslu, eval=FALSE, purl = FALSE}
##### Lagrange multipliers
par(mfrow=c(2,4))
lapply(names(solListWSLU), function(Y){cat(Y, "\n");solListWSLU[[Y]]$lambdaCol})
```

##### Biplots

```{r LU biplots, results='hide', eval=FALSE, purl = FALSE}
par(mfrow=c(2,4), mai = c(2,2,2,4))
par(mfrow=c(1,1))
lapply(names(solListWSLU), function(Y){plotRCM(solListWSLU[[Y]], samColour = sampleSigNBLU, main=Y, biplot=TRUE, Dim = c(1,3), taxColour=tmp, taxLegPos = "bottomright", abundLeg=TRUE, taxCol  =cols, taxInset = c(-1.45,-0.75), libCex = 0.65, libInset = c(-0.85,0), arrowFrac = 0.05)})
```

##### Influence function

Take a look at the influence function values on the psis

```{r LU influence measures: psis, eval=FALSE, purl = FALSE}
infl1 = lapply(solListWSLU, function(x, i){
  infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
  id=LUs(infl) > quantile( LUs(infl),0.995)
  Xid = x$X[id]
  abunds = colSums(x$X)/sum(x$X)
  libsizes = rowSums(x$X)
  list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},1)
par(mfrow=c(1,2))
lapply(names(infl1), function(x){
  with(infl1[[x]], plot(abunds, colSums(id), log="x", main=x, ylab="Number of very influential observations"))
}) #High abundances, larger influence
lapply(names(infl1), function(x){
  idTmp = sample(seq_along(infl1[[x]]$abunds), 200) #For speed
  with(infl1[[x]], plot(rep(abunds[idTmp], nrow(infl)), c(t(infl[,idTmp])), log="x", main=x, xlab="abundance",ylab = "Influence")) #Very heavy!
})
lapply(names(infl1), function(x){
  idTmp = sample(seq_along(infl1[[x]]$libsizes), 200) #For speed
  with(infl1[[x]], plot(rep(libsizes[idTmp], ncol(infl)), c(infl[idTmp,]), log="", main=x, xlab ="Library sizes",ylab = "Influence")) #Very heavy!
})
#Expectations
lapply(names(infl1), function(x){
  with(infl1[[x]], {rbind(quantile(outer(libsizes, abunds)),
                          quantile(outer(libsizes, abunds)[id])) })
}) #Influential observations have high expectations
lapply(names(infl1), function(x){
  cat(x, "\n")
  with(infl1[[x]], table(Xid))
}) # ... but are very often zero!

infl2 = lapply(solListWSLU, function(x, i){
  infl=with(x, NBpsiInfl(psi = psis[i], X = X, cMat = cMat[i,,drop=FALSE], rMat = rMat[,i, drop=FALSE], muMarg = outer(rowSums(X), colSums(X)/sum(X)), theta = thetas))
  id=LUs(infl) > quantile( LUs(infl),0.995)
  Xid = x$X[id]
  abunds = colSums(x$X)/sum(x$X)
  libsizes = rowSums(x$X)
  list(infl = infl, id = id, Xid = Xid, abunds=abunds, libsizes=libsizes)
},2)
par(mfrow=c(1,2))
lapply(names(infl2), function(x){
  with(infl2[[x]], plot(abunds, colSums(id), log="x", main=x))
})
#In the second dimension, larger abundance means larger influence on the psis too, and no outliers
lapply(names(infl2), function(x){
  with(infl2[[x]], plot(rep(abunds, nrow(infl)), c(t(infl)), log="x", main=x, xlab="abundance",ylab = "Influence")) #Very heavy!
}) #Outlier in unifmarg weigthing scheme
lapply(names(infl2), function(x){
  with(infl2[[x]], plot(rep(libsizes, ncol(infl)), c(infl), log="x", main=x, xlab="Library size",ylab = "Influence")) #Very heavy! #No trend in terms of library sizes visible from here, but their distribution should be log-normal perhaps
})
#Expectations
lapply(names(infl2), function(x){
  with(infl2[[x]], {rbind(quantile(outer(libsizes, abunds)),
                          quantile(outer(libsizes, abunds)[id])) })
})
lapply(names(infl2), function(x){
  cat(x, "\n")
  with(infl2[[x]], tLUle(Xid))
})
```

The most influential observations are zero counts in highly abundant species and high libsizes (i.e. high expectations)! Note that we cen derive this from the influence function that does not even depend on the weights!

  Influence on the colScores

```{r LU influence colscores, eval=FALSE, purl = FALSE}
inflColList1 = lapply(solListWSLU, function(x){
  try(with(x, NBcolInfl(X, psis, cMat, rMat,thetas , colWeights = if(is.matrix(colWeights)) colWeights else{ cbind(colWeights, colWeights, colWeights)} , k=1 , lambdaCol)))
})
#Look at the signal from the first group
Id1 = which(taxaSigNBLU=="Signal 1")[1]

inflCol1 = lapply(inflColList1, function(x){getInflCol(x$score, x$InvJac, Id1)})
lapply(inflCol1, function(x){
  boxplot(x[,Id1]~sampleSigNBLU, las=2, ylab="Influence")
})
```

As expected, the samples with the signal have the highest impact on the score

```{r LU influence rowscores LU, eval=FALSE, purl = FALSE}
inflRowList1 = lapply(solListWSLU, function(x){
  try(with(x, NBrowInfl(X, psis, cMat, rMat,thetas , rowWeights , k=1 , lambdaRow)))
})
#Only look at most extreme colscores
#Cannot calculate influence functions for the interesting cases (with the outliers)
inflRowList1=inflRowList1[sapply(inflRowList1, class)=="list"]
#The influence on the first row score

par(mfrow=c(1,2))
idRow = 10
libSizes = rowSums(solListWSLU[[1]]$X)
inflRow1 = lapply(inflRowList1, function(x){getInflRow(x$score, x$InvJac, idRow)})
lapply(names(inflRow1), function(x){
  plot(y=inflRow1[[x]][,idRow], libSizes, log="x", main=x, col=(libSizes ==libSizes[idRow])+1, ylab="Influence")
})
```

Evidently, the observation from the sample itself has the largest influence. In the marginal weighting scheme for the libsizes, the larger libsizes sometimes get a larger influence, although the effect may be small.

Achieved tolerances

```{r toleranceslu, eval=FALSE, purl = FALSE}
tolAchPsi = with(syntNBSigmargmarg_LU, lapply(1:k, function(y){sapply( 2:iter[y], function(x){LUs(1-psiRec[y,x]/psiRec[y,x-1])})}))
par(mfrow=c(1,3))
lapply(tolAchPsi, function(x){plot(x, log="y");LUline(h=0.01, col="red"); abline(h=0.001, col="blue")})

with(syntNBSigmargmarg_LU,plot(psiRec[3,1:iter[3]])) #Weird jumps, don't stop too early!
```
--->
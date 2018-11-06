#' A modified version of the vgam.fit() function from the VGAM pacakage
#' @param x the model matrix
#' @param y the outcome vector
#' @param psi the importance parameter
#' @param w a vector of weights
#' @param mf the model frame
#' @param coefstart vector of initial parameter estimates
#' @param offset the offset, based on previous dimensions
#' @param family the distributional family
#' @param control a list of control parameters
#' @param Terms,smooth.labels,function.name,... Other, orginal arguments of the vgam.fit() function
#'
#' @importFrom stats .getXlevels is.empty.model model.offset model.response

vgam.fit.edit = function (x, y, psi, w = rep_len(1, nrow(x)), mf, coefstart = NULL, offset = 0, family, control = vgam.control(), Terms, smooth.labels, function.name = "vgam", ...)
{
  if (is.null(criterion <- control$criterion))
    criterion <- "coefficients"
  eff.n <- nrow(x)
  specialCM = constraints = etastart = mustart = NULL
  post <- list()
  check.rank <- control$Check.rank
  epsilon <- control$epsilon
  maxit <- control$maxit
  save.weights <- control$save.weights
  trace <- control$trace
  bf.maxit <- control$bf.maxit
  bf.epsilon <- control$bf.epsilon
  se.fit <- control$se.fit
  minimize.criterion <- control$min.criterion
  fv <- NULL
  n <- nrow(x)
  old.coeffs <- coefstart
  intercept.only <- ncol(x) == 1 && colnames(x) == "(Intercept)"
  y.names <- predictors.names <- NULL
  n.save <- n
  if (length(slot(family, "initialize")))
    eval(slot(family, "initialize"))
  if (length(etastart)) {
    eta <- etastart
    mu <- if (length(mustart))
      mustart
    else slot(family, "linkinv")(eta, extra = NULL)
  }
  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra = NULL)
    }
    else {
      warning("argument 'mustart' assigned a value ", "but there is no 'linkfun' slot to use it")
    }
  }
  validparams <- validfitted <- TRUE
  if (length(body(slot(family, "validparams"))))
    validparams <- slot(family, "validparams")(eta, y = y,
                                               extra = NULL)
  if (length(body(slot(family, "validfitted"))))
    validfitted <- slot(family, "validfitted")(mu, y = y,
                                               extra = NULL)
  if (!(validparams && validfitted))
    stop("could not obtain valid initial values. ", "Try using 'etastart', 'coefstart' or 'mustart', else ",
         "family-specific arguments such as 'imethod'.")
  if (length(family@constraints))
    eval(slot(family, "constraints"))
  Hlist <- process.constraints(constraints, x = x, M = M, specialCM = specialCM, Check.cm.rank = control$Check.cm.rank)
  ncolHlist <- unlist(lapply(Hlist, ncol))
  M <- NCOL(eta)
    smooth.frame <- mf
    assignx <- attr(x, "assign")
    which <- assignx[smooth.labels]
    bf <- "s.vam"
    bf.call <- parse(text = paste("s.vam(x, z, wz, tfit$smomat, which, tfit$smooth.frame,",
                                  "bf.maxit, bf.epsilon, trace, se = se.fit, X.vlm.save, ",
                                  "Hlist, ncolHlist, M = M, qbig = qbig, Umat = U, ",
                                  "all.knots = control$all.knots, nk = control$nk)",
                                  sep = ""))[[1]]
    qbig <- sum(ncolHlist[smooth.labels])
    smomat <- matrix(0, n, qbig)
    dy <- if (is.matrix(y))
      dimnames(y)[[1]]
    else names(y)
    d2 <- if (is.null(predictors.names))
      paste("(Additive predictor ", 1:M, ")", sep = "")
    else predictors.names
    dimnames(smomat) <- list(dy, vlabel(smooth.labels, ncolHlist[smooth.labels],
                                        M))
    tfit <- list(smomat = smomat, smooth.frame = smooth.frame)

  X.vlm.save <- lm2vlm.model.matrix(x, Hlist, xij = control$xij,
                                    Xm2 = NULL)
  if (length(coefstart)) {
    eta <- if (ncol(X.vlm.save) > 1) {
      matrix(X.vlm.save %*% coefstart, n, M, byrow = TRUE) +
        offset
    }
    else {
      matrix(X.vlm.save * coefstart, n, M, byrow = TRUE) +
        offset
    }
    if (M == 1)
      eta <- c(eta)
    mu <- slot(family, "linkinv")(eta, extra = NULL)
  }
  if (criterion != "coefficients") {
    tfun <- slot(family, criterion)
  }
  iter <- 1
  new.crit <- switch(criterion, coefficients = 1, tfun(mu = mu,
                                                       y = y, w = w, res = FALSE, eta = eta, extra = NULL))
  old.crit <- ifelse(minimize.criterion, 10 * new.crit + 10,
                     -10 * new.crit - 10)
  deriv.mu <- eval(slot(family, "deriv"))
  wz <- eval(slot(family, "weight"))
  if (control$checkwz)
    wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)
  U <- vchol(wz, M = M, n = n, silent = !trace)
  tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
  z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset
  nrow.X.vlm <- nrow(X.vlm.save)
  ncol.X.vlm <- ncol(X.vlm.save)
  fully.cvged <- FALSE
  for (iter.outer in 1:control$Maxit.outer) {
    if (fully.cvged)
      break
    iter <- 1
    one.more <- TRUE
    while (one.more) {
      tfit <- eval(bf.call)
      fv <- tfit$fitted.values
      new.coeffs <- tfit$coefficients
      if (length(slot(family, "middle")))
        eval(slot(family, "middle"))
      eta <- fv * psi + offset #Added psi
      mu <- slot(family, "linkinv")(eta, extra = NULL)
      if (length(family@middle2))
        eval(family@middle2)
      old.crit <- new.crit
      new.crit <- switch(criterion, coefficients = new.coeffs,
                         tfun(mu = mu, y = y, w = w, res = FALSE, eta = eta,
                              extra = NULL))
      if (trace) {
        cat("VGAM ", bf, " loop ", iter, ": ", criterion,
            "= ")
        UUUU <- switch(criterion, coefficients = format(new.crit,
                                                        dig = round(1 - log10(epsilon))), format(new.crit,
                                                                                                 dig = max(4, round(-0 - log10(epsilon) + log10(sqrt(eff.n))))))
        switch(criterion, coefficients = {
          if (length(new.crit) > 2) cat("\n")
          cat(UUUU, fill = TRUE, sep = ", ")
        }, cat(UUUU, fill = TRUE, sep = ", "))
      }
      one.more <- eval(control$convergence)
      if (!is.logical(one.more))
        one.more <- FALSE
      if (one.more) {
        iter <- iter + 1
        deriv.mu <- eval(slot(family, "deriv"))
        wz <- eval(slot(family, "weight"))
        if (control$checkwz)
          wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)
        U <- vchol(wz, M = M, n = n, silent = !trace)
        tvfor <- vforsub(U, as.matrix(deriv.mu), M = M,
                         n = n)
        z <- eta + vbacksub(U, tvfor, M = M, n = n) -
          offset
      }
      else {
        fully.cvged <- TRUE
      }
      old.coeffs <- new.coeffs
    }
  }
  if (maxit > 1 && iter >= maxit && !control$noWarning)
    warning("convergence not obtained in ", maxit, " IRLS iterations")
  if (control$Maxit.outer > 1 && iter.outer >= control$Maxit.outer &&
      !control$noWarning)
    warning("convergence not obtained in ", control$Maxit.outer,
            " outer iterations")
  dnrow.X.vlm <- labels(X.vlm.save)
  xnrow.X.vlm <- dnrow.X.vlm[[2]]
  ynrow.X.vlm <- dnrow.X.vlm[[1]]
  if (length(slot(family, "fini")))
    eval(slot(family, "fini"))
  if (M > 1)
    fv <- matrix(fv, n, M)
  final.coefs <- new.coeffs
  asgn <- attr(X.vlm.save, "assign")
  names(final.coefs) <- xnrow.X.vlm
  if (!is.null(tfit$rank)) {
    rank <- tfit$rank
  }
  else {
    rank <- NCOL(x)
  }
  cnames <- xnrow.X.vlm
  R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
  R[lower.tri(R)] <- 0
  attributes(R) <- list(dim = c(ncol.X.vlm, ncol.X.vlm), dimnames = list(cnames,
                                                                         cnames), rank = rank)
  dim(fv) <- c(n, M)
  dn <- labels(x)
  yn <- dn[[1]]
  xn <- dn[[2]]
  wresiduals <- z - fv
  if (M == 1) {
    fv <- as.vector(fv)
    wresiduals <- as.vector(wresiduals)
    names(wresiduals) <- names(fv) <- yn
  }
  else {
    dimnames(wresiduals) <- dimnames(fv) <- list(yn, predictors.names)
  }
  if (is.matrix(mu)) {
    if (length(dimnames(y)[[2]])) {
      y.names <- dimnames(y)[[2]]
    }
    if (length(dimnames(mu)[[2]])) {
      y.names <- dimnames(mu)[[2]]
    }
    dimnames(mu) <- list(yn, y.names)
  }
  else {
    names(mu) <- names(fv)
  }
  tfit$fitted.values <- NULL
  fit <- structure(c(tfit, list(assign = asgn, constraints = Hlist,
                                control = control, fitted.values = mu, formula = as.vector(attr(Terms,
                                                                                                "formula")), iter = iter, offset = offset, rank = rank,
                                R = R, terms = Terms)))
  if (!se.fit) {
    fit$varmat <- NULL
  }
  if (M == 1) {
    wz <- as.vector(wz)
  }
  fit$weights <- if (save.weights)
    wz
  else NULL
  misc <- list(colnames.x = xn, colnames.X.vlm = xnrow.X.vlm,
               criterion = criterion, function.name = function.name,
               intercept.only = intercept.only, predictors.names = predictors.names,
               M = M, n = n,  nonparametric = TRUE,
               nrow.X.vlm = nrow.X.vlm, orig.assign = attr(x, "assign"),
               p = ncol(x), ncol.X.vlm = ncol.X.vlm, ynames = colnames(y)) #new.assign = VGAM:::new.assign(x, NewHlist),
    misc$smooth.labels <- smooth.labels
  crit.list <- list()
  if (criterion != "coefficients")
    crit.list[[criterion]] <- fit[[criterion]] <- new.crit
  for (ii in names(.min.criterion.VGAM)) {
    if (ii != criterion && any(slotNames(family) == ii) &&
        length(body(slot(family, ii)))) {
      fit[[ii]] <- crit.list[[ii]] <- (slot(family, ii))(mu = mu,
                                                         y = y, w = w, res = FALSE, eta = eta, extra = NULL)
    }
  }
  if (w[1] != 1 || any(w != w[1]))
    fit$prior.weights <- w
  if (length(slot(family, "last")))
    eval(slot(family, "last"))
   fit$nl.chisq <- 0#VGAM:::vgam.nlchisq(fit$qr, fit$resid, wz = wz,smomat = fit$smomat, deriv = deriv.mu, U = U, smooth.labels,                                        attr(x, "assign"), M = M, n = n, constraints = Hlist)
  fit$qr <- NULL
  fit$misc <- NULL
  structure(c(fit, list(predictors = fv, contrasts = attr(x,
                                                          "contrasts"), control = control, crit.list = crit.list,
                       family = family, iter = iter, misc = misc,
                        post = post, x = x, y = y)), vclass = slot(family, "vfamily"))
}

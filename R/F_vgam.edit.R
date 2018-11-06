#' Some modified functions of the VGAM pacakage, to account for the particular model structure
vgam.fit.edit = function (x, y, psi, w = rep_len(1, nrow(x)), mf, Xm2 = NULL, Ym2 = NULL,
          etastart = NULL, mustart = NULL, coefstart = NULL, offset = 0,
          family, control = vgam.control(), qr.arg = FALSE, constraints = NULL,
          extra = NULL, Terms, nonparametric, smooth.labels, function.name = "vgam",
          sm.osps.list = NULL, ...)
{
  mgcvvgam <- length(sm.osps.list) > 0
  if (is.null(criterion <- control$criterion))
    criterion <- "coefficients"
  eff.n <- nrow(x)
  specialCM <- NULL
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
    else slot(family, "linkinv")(eta, extra = extra)
  }
  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra = extra)
    }
    else {
      warning("argument 'mustart' assigned a value ", "but there is no 'linkfun' slot to use it")
    }
  }
  validparams <- validfitted <- TRUE
  if (length(body(slot(family, "validparams"))))
    validparams <- slot(family, "validparams")(eta, y = y,
                                               extra = extra)
  if (length(body(slot(family, "validfitted"))))
    validfitted <- slot(family, "validfitted")(mu, y = y,
                                               extra = extra)
  if (!(validparams && validfitted))
    stop("could not obtain valid initial values. ", "Try using 'etastart', 'coefstart' or 'mustart', else ",
         "family-specific arguments such as 'imethod'.")
  M <- NCOL(eta)
  if (length(family@constraints))
    eval(slot(family, "constraints"))
  Hlist <- process.constraints(constraints, x = x, M = M, specialCM = specialCM,
                               Check.cm.rank = control$Check.cm.rank)
  ncolHlist <- unlist(lapply(Hlist, ncol))
  if (nonparametric) {
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
  }
  else {
    bf.call <- expression(vlm.wfit(xmat = X.vlm.save, z,
                                   Hlist = NULL, U = U, matrix.out = FALSE, is.vlmX = TRUE,
                                   qr = qr.arg, xij = NULL))
    bf <- "vlm.wfit"
  }
  X.vlm.save <- lm2vlm.model.matrix(x, Hlist, xij = control$xij,
                                    Xm2 = Xm2)
  if (mgcvvgam) {
    Xvlm.aug <- get.X.VLM.aug(constraints = constraints,
                              sm.osps.list = sm.osps.list)
    first.sm.osps <- TRUE
  }
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
    mu <- slot(family, "linkinv")(eta, extra = extra)
  }
  if (criterion != "coefficients") {
    tfun <- slot(family, criterion)
  }
  iter <- 1
  new.crit <- switch(criterion, coefficients = 1, tfun(mu = mu,
                                                       y = y, w = w, res = FALSE, eta = eta, extra = extra))
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
  if (!nonparametric && nrow.X.vlm < ncol.X.vlm)
    stop("There are ", ncol.X.vlm, " parameters but only ",
         nrow.X.vlm, " observations")
  if (mgcvvgam) {
    bf.call <- expression(vlm.wfit(xmat = X.vlm.save, z,
                                   Hlist = Hlist, U = U, matrix.out = FALSE, is.vlmX = TRUE,
                                   qr = qr.arg, xij = NULL, Xvlm.aug = Xvlm.aug, sm.osps.list = sm.osps.list,
                                   constraints = constraints, first.sm.osps = first.sm.osps,
                                   control = control, trace = trace))
    bf <- "vlm.wfit"
  }
  fully.cvged <- FALSE
  for (iter.outer in 1:control$Maxit.outer) {
    if (fully.cvged)
      break
    if (trace && mgcvvgam) {
      cat("VGAM outer iteration ", iter.outer, " =============================================\n")
      flush.console()
    }
    iter <- 1
    one.more <- TRUE
    sm.osps.list$fixspar <- sm.osps.list$orig.fixspar
    while (one.more) {
      tfit <- eval(bf.call)
      if (mgcvvgam) {
        first.sm.osps <- tfit$first.sm.osps
        Xvlm.aug <- tfit$Xvlm.aug
        sm.osps.list <- tfit$sm.osps.list
        if (control$Maxit.outer > 1)
          sm.osps.list$fixspar <- rep_len(TRUE, length(sm.osps.list$fixspar))
        magicfit <- tfit$magicfit
      }
      fv <- tfit$fitted.values
      if (mgcvvgam) {
        fv <- head(fv, n * M)
      }
      new.coeffs <- tfit$coefficients
      if (length(slot(family, "middle")))
        eval(slot(family, "middle"))
      eta <- fv * psi + offset #Added psi
      mu <- slot(family, "linkinv")(eta, extra = extra)
      if (length(family@middle2))
        eval(family@middle2)
      old.crit <- new.crit
      new.crit <- switch(criterion, coefficients = new.coeffs,
                         tfun(mu = mu, y = y, w = w, res = FALSE, eta = eta,
                              extra = extra))
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
      flush.console()
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
        fully.cvged <- if (mgcvvgam)
          (iter <= 2)
        else TRUE
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
  if (!nonparametric && check.rank && rank < ncol.X.vlm)
    stop("vgam() only handles full-rank models (currently)")
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
  if (qr.arg) {
    fit$qr <- tfit$qr
    dimnames(fit$qr$qr) <- dnrow.X.vlm
  }
  if (!mgcvvgam && !se.fit) {
    fit$varmat <- NULL
  }
  if (M == 1) {
    wz <- as.vector(wz)
  }
  fit$weights <- if (save.weights)
    wz
  else NULL
  NewHlist <- process.constraints(constraints, x, M, specialCM = specialCM,
                                  by.col = FALSE)
  misc <- list(colnames.x = xn, colnames.X.vlm = xnrow.X.vlm,
               criterion = criterion, function.name = function.name,
               intercept.only = intercept.only, predictors.names = predictors.names,
               M = M, n = n, new.assign = VGAM:::new.assign(x, NewHlist), nonparametric = nonparametric,
               nrow.X.vlm = nrow.X.vlm, orig.assign = attr(x, "assign"),
               p = ncol(x), ncol.X.vlm = ncol.X.vlm, ynames = colnames(y))
  if (!mgcvvgam && se.fit && length(fit$s.xargument)) {
    misc$varassign <- VGAM:::varassign(Hlist, names(fit$s.xargument))
  }
  if (nonparametric) {
    misc$smooth.labels <- smooth.labels
  }
  if (mgcvvgam) {
    misc$Xvlm.aug <- Xvlm.aug
    misc$sm.osps.list <- sm.osps.list
    misc$magicfit <- magicfit
    misc$iter.outer <- iter.outer
  }
  crit.list <- list()
  if (criterion != "coefficients")
    crit.list[[criterion]] <- fit[[criterion]] <- new.crit
  for (ii in names(.min.criterion.VGAM)) {
    if (ii != criterion && any(slotNames(family) == ii) &&
        length(body(slot(family, ii)))) {
      fit[[ii]] <- crit.list[[ii]] <- (slot(family, ii))(mu = mu,
                                                         y = y, w = w, res = FALSE, eta = eta, extra = extra)
    }
  }
  if (w[1] != 1 || any(w != w[1]))
    fit$prior.weights <- w
  if (length(slot(family, "last")))
    eval(slot(family, "last"))
  if (!is.null(fit$smomat)) {
    fit$nl.chisq <- VGAM:::vgam.nlchisq(fit$qr, fit$resid, wz = wz,
                                 smomat = fit$smomat, deriv = deriv.mu, U = U, smooth.labels,
                                 attr(x, "assign"), M = M, n = n, constraints = Hlist)
  }
  if (!qr.arg) {
    fit$qr <- NULL
  }
  fit$misc <- NULL
  structure(c(fit, list(predictors = fv, contrasts = attr(x,
                                                          "contrasts"), control = control, crit.list = crit.list,
                        extra = extra, family = family, iter = iter, misc = misc,
                        post = post, x = x, y = y)), vclass = slot(family, "vfamily"))
}

vgam.edit = function (formula, psi, family = stop("argument 'family' needs to be assigned"),
                      data = list(), weights = NULL, subset = NULL, na.action = na.fail,
                      etastart = NULL, mustart = NULL, coefstart = NULL, control = vgam.control(criterion = "loglikelihood",...),
                      offset = NULL, method = "vgam.fit", model = FALSE, x.arg = TRUE,
                      y.arg = TRUE, contrasts = NULL, constraints = NULL, extra = list(),
                      form2 = NULL, qr.arg = FALSE, smart = TRUE, ...)
{
  dataname <- as.character(substitute(data))
  function.name <- "vgam"
  ocall <- match.call()
  if (smart)
    setup.smart("write")
  if (missing(data))
    data <- environment(formula)
  mtsave <- terms(formula, specials = c("s", "sm.os", "sm.ps"),
                  data = data)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method, model.frame = return(mf), vgam.fit = 1, stop("invalid 'method': ",
                                                              method))
  mt <- attr(mf, "terms")
  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")
  x <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(y), 0)
  attr(x, "assign") <- VGAM:::attrassigndefault(x, mt)
  if (!is.null(form2)) {
    if (!is.null(subset))
      stop("argument 'subset' cannot be used when ", "argument 'form2' is used")
    retlist <- shadowvgam(formula = form2, family = family,
                          data = data, na.action = na.action, control = vgam.control(...),
                          method = method, model = model, x.arg = x.arg, y.arg = y.arg,
                          contrasts = contrasts, constraints = constraints,
                          extra = extra, qr.arg = qr.arg)
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2
    if (length(Ym2)) {
      if (NROW(Ym2) != NROW(y))
        stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2)) {
      if (NROW(Xm2) != NROW(x))
        stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  }
  else {
    Xm2 <- Ym2 <- NULL
  }
  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0
  mf2 <- mf
  if (!missing(subset)) {
    mf2$subset <- NULL
    mf2 <- eval(mf2, parent.frame())
    spars2 <- lapply(mf2, attr, "spar")
    dfs2 <- lapply(mf2, attr, "df")
    sx2 <- lapply(mf2, attr, "s.xargument")
    for (ii in seq_along(mf)) {
      if (length(sx2[[ii]])) {
        attr(mf[[ii]], "spar") <- spars2[[ii]]
        attr(mf[[ii]], "dfs2") <- dfs2[[ii]]
        attr(mf[[ii]], "s.xargument") <- sx2[[ii]]
      }
    }
    rm(mf2)
  }
  w <- model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  }
  else if (NCOL(w) == 1 && any(w < 0))
    stop("negative weights not allowed")
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }
  eval(vcontrol.expression)
  n <- dim(x)[1]
  if (length(slot(family, "first")))
    eval(slot(family, "first"))
  aa <- attributes(mtsave)
  smoothers <- aa$specials
  mgcv.sm.os <- length(smoothers$sm.os) > 0
  mgcv.sm.ps <- length(smoothers$sm.ps) > 0
  mgcv.sm.PS <- length(smoothers$sm.PS) > 0
  any.sm.os.terms <- mgcv.sm.os
  any.sm.ps.terms <- mgcv.sm.ps || mgcv.sm.PS
  mgcv.s <- length(smoothers$s) > 0
  if ((any.sm.os.terms || any.sm.ps.terms) && mgcv.s)
    stop("cannot include both s() and any of sm.os() or ",
         "sm.ps() (or sm.PS()) terms in the formula")
  if (any.sm.os.terms && any.sm.ps.terms)
    stop("cannot include both sm.os() and ", "sm.ps() (or sm.PS()) terms in the formula")
  nonparametric <- length(smoothers$s) > 0
  if (nonparametric) {
    ff <- apply(aa$factors[smoothers[["s"]], , drop = FALSE],
                2, any)
    smoothers[["s"]] <- if (any(ff))
      seq(along = ff)[aa$order == 1 & ff]
    else NULL
    smooth.labels <- aa$term.labels[unlist(smoothers)]
  }
  else {
    function.name <- "vglm"
  }
  are.sm.os.terms <- length(smoothers$sm.os) > 0
  are.sm.ps.terms <- (length(smoothers$sm.ps) + length(smoothers$sm.PS)) >
    0
  if (are.sm.os.terms || are.sm.ps.terms) {
    control$criterion <- "coefficients"
    if (length(smoothers$sm.os) > 0) {
      ff.sm.os <- apply(aa$factors[smoothers[["sm.os"]],
                                   , drop = FALSE], 2, any)
      smoothers[["sm.os"]] <- if (any(ff.sm.os))
        seq(along = ff.sm.os)[aa$order == 1 & ff.sm.os]
      else NULL
      smooth.labels <- aa$term.labels[unlist(smoothers)]
    }
    if (length(smoothers$sm.ps) > 0) {
      ff.sm.ps <- apply(aa$factors[smoothers[["sm.ps"]],
                                   , drop = FALSE], 2, any)
      smoothers[["sm.ps"]] <- if (any(ff.sm.ps))
        seq(along = ff.sm.ps)[aa$order == 1 & ff.sm.ps]
      else NULL
      smooth.labels <- aa$term.labels[unlist(smoothers)]
    }
    assignx <- attr(x, "assign")
    which.X.sm.osps <- assignx[smooth.labels]
    Data <- mf[, names(which.X.sm.osps), drop = FALSE]
    attr(Data, "class") <- NULL
    S.arg <- lapply(Data, attr, "S.arg")
    sparlist <- lapply(Data, attr, "spar")
    ridge.adj <- lapply(Data, attr, "ridge.adj")
    fixspar <- lapply(Data, attr, "fixspar")
    ps.int <- lapply(Data, attr, "ps.int")
    knots <- lapply(Data, attr, "knots")
    term.labels <- aa$term.labels
  }
  sm.osps.list <- if (any.sm.os.terms || any.sm.ps.terms)
    list(indexterms = if (any.sm.os.terms) ff.sm.os else ff.sm.ps,
         intercept = aa$intercept, which.X.sm.osps = which.X.sm.osps,
         S.arg = S.arg, sparlist = sparlist, ridge.adj = ridge.adj,
         term.labels = term.labels, fixspar = fixspar, orig.fixspar = fixspar,
         ps.int = ps.int, knots = knots, assignx = assignx)
  else NULL
  fit <- vgam.fit.edit(x = x, y = y, psi = psi, w = w, mf = mf, Xm2 = Xm2,
                  Ym2 = Ym2, etastart = etastart, mustart = mustart, coefstart = coefstart,
                  offset = offset, family = family, control = control,
                  constraints = constraints, extra = extra, qr.arg = qr.arg,
                  Terms = mtsave, nonparametric = nonparametric, smooth.labels = smooth.labels,
                  function.name = function.name, sm.osps.list = sm.osps.list,
                  ...)
  if (is.Numeric(fit$nl.df) && any(fit$nl.df < 0)) {
    fit$nl.df[fit$nl.df < 0] <- 0
  }
  if (!is.null(fit[["smooth.frame"]])) {
    fit <- fit[-1]
  }
  else {
  }
  fit$smomat <- NULL
  fit$call <- ocall
  if (model)
    fit$model <- mf
  if (!x.arg)
    fit$x <- NULL
  if (!y.arg)
    fit$y <- NULL
  if (nonparametric)
    fit$misc$smooth.labels <- smooth.labels
  fit$misc$dataname <- dataname
  if (smart)
    fit$smart.prediction <- get.smart.prediction()
  answer <- new(if (any.sm.os.terms || any.sm.ps.terms)
    "pvgam"
    else "vgam", assign = attr(x, "assign"), call = fit$call,
    coefficients = fit$coefficients, constraints = fit$constraints,
    criterion = fit$crit.list, df.residual = fit$df.residual,
    dispersion = 1, family = fit$family, misc = fit$misc,
    model = if (model)
      mf
    else data.frame(), R = fit$R, rank = fit$rank, residuals = as.matrix(fit$residuals),
    ResSS = fit$ResSS, smart.prediction = as.list(fit$smart.prediction),
    terms = list(terms = fit$terms))
  if (!smart)
    answer@smart.prediction <- list(smart.arg = FALSE)
  if (qr.arg) {
    class(fit$qr) <- "list"
    slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts")))
    slot(answer, "contrasts") <- attr(x, "contrasts")
  if (length(fit$fitted.values))
    slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
    list(aaa)
  else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$weights))
    slot(answer, "weights") <- as.matrix(fit$weights)
  if (x.arg)
    slot(answer, "x") <- x
  if (length(fit$misc$Xvlm.aug)) {
    slot(answer, "ospsslot") <- list(Xvlm.aug = fit$misc$Xvlm.aug,
                                     sm.osps.list = fit$misc$sm.osps.list, magicfit = fit$misc$magicfit,
                                     iter.outer = fit$misc$iter.outer)
    fit$misc$Xvlm.aug <- NULL
    fit$misc$sm.osps.list <- NULL
    fit$misc$magicfit <- NULL
    fit$misc$iter.outer <- NULL
  }
  if (x.arg && length(Xm2))
    slot(answer, "Xm2") <- Xm2
  if (y.arg && length(Ym2))
    slot(answer, "Ym2") <- as.matrix(Ym2)
  if (!is.null(form2))
    slot(answer, "callXm2") <- retlist$call
  answer@misc$formula <- formula
  answer@misc$form2 <- form2
  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(fit$y)
  answer@misc$formula <- formula
  slot(answer, "control") <- fit$control
  if (length(fit$extra)) {
    slot(answer, "extra") <- fit$extra
  }
  slot(answer, "iter") <- fit$iter
  slot(answer, "post") <- fit$post
  fit$predictors <- as.matrix(fit$predictors)
  dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]],
                                   fit$misc$predictors.names)
  slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights))
    slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)
  if (nonparametric) {
    slot(answer, "Bspline") <- fit$Bspline
    slot(answer, "nl.chisq") <- fit$nl.chisq
    if (is.Numeric(fit$nl.df))
      slot(answer, "nl.df") <- fit$nl.df
    slot(answer, "spar") <- fit$spar
    slot(answer, "s.xargument") <- fit$s.xargument
    if (length(fit$varmat)) {
      slot(answer, "var") <- fit$varmat
    }
  }
  if (length(fit$effects))
    slot(answer, "effects") <- fit$effects
  if (nonparametric && is.buggy.vlm(answer)) {
    warning("some s() terms have constraint matrices that have columns",
            " which are not orthogonal;", " try using sm.os() or sm.ps() instead of s().")
  }
  else {
  }
  answer
}

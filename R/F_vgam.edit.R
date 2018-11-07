#' A modified varsion of the vgam() function in the VGAM pacakage, to account for the particular model structure
#' @param formula The model formula
#' @param psi the importance parameter
#' @param family the family object
#' @param data the data frame
#' @param coefstart The starting values for the coefficients
#' @param control a list of control arguments
#' @param offset the offset, depending on lower dimensions
#' @param method A character string, the fitting method
#' @param ... Additional arguments, passed on to the vgam.fit.edit() function
#'
#' @return an object of class vgam
vgam.edit = function (formula, psi, family = stop("argument 'family' needs to be assigned"),
                      data = list(), coefstart = NULL, control = vgam.control(criterion = "loglikelihood",...),
                      offset = NULL, method = "vgam.fit",...)
{
  dataname <- as.character(substitute(data))
  function.name <- "vgam"
  ocall <- match.call()
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
  attrassigndefault = function (mmat, tt)
  {
    if (!inherits(tt, "terms"))
      stop("need terms object")
    aa <- attr(mmat, "assign")
    if (is.null(aa))
      stop("argument is not really a model matrix")
    ll <- attr(tt, "term.labels")
    if (attr(tt, "intercept") > 0)
      ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
    split(order(aa), aaa)
  }
  attr(x, "assign") <- attrassigndefault(x, mt)
  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0
  mf2 <- mf
  w <- rep_len(1, nrow(mf))
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
  fit <- vgam.fit.edit(x = x, y = y, psi = psi, w = w, mf = mf,
                  offset = offset, family = family, control = control,
               extra = list(), qr.arg = FALSE, coefstart = coefstart,
                  Terms = mtsave, nonparametric = nonparametric, smooth.labels = smooth.labels,
                  function.name = function.name,
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
  fit$model <- mf
  if (nonparametric)
    fit$misc$smooth.labels <- smooth.labels
  fit$misc$dataname <- dataname
  answer <- new("vgam", assign = attr(x, "assign"), call = fit$call,
    coefficients = fit$coefficients, constraints = fit$constraints,
    criterion = fit$crit.list, df.residual = fit$df.residual,
    dispersion = 1, family = fit$family, misc = fit$misc, R = fit$R, rank = fit$rank, residuals = as.matrix(fit$residuals),
    ResSS = fit$ResSS, smart.prediction = as.list(fit$smart.prediction),
    terms = list(terms = fit$terms))
  answer@smart.prediction <- list(smart.arg = FALSE)
  if (length(fit$fitted.values))
    slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
    list(aaa)
  else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$misc$Xvlm.aug)) {
    slot(answer, "ospsslot") <- list(Xvlm.aug = fit$misc$Xvlm.aug,
                                     sm.osps.list = fit$misc$sm.osps.list, magicfit = fit$misc$magicfit,
                                     iter.outer = fit$misc$iter.outer)
    fit$misc$Xvlm.aug <- NULL
    fit$misc$sm.osps.list <- NULL
    fit$misc$magicfit <- NULL
    fit$misc$iter.outer <- NULL
  }
  answer@misc$formula <- formula
  answer@misc$form2 <- NULL
  answer@misc$formula <- formula
  slot(answer, "control") <- fit$control
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
  else {
  }
  list(coef = coef(answer), spline = answer@Bspline[[1]])
}


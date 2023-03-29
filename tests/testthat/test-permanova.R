context("permanova")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(150)],
                    prune_samples(sample_names(Zeller)[seq_len(100)], Zeller))
rcmFit = RCM(tmpPhy)
test_that("permanova runs as expected", {
  expect_is(permanova(rcmFit, "Diagnosis"), "list")
  expect_is(permanova(rcmFit, # Bogus, user-supplied grouping factor
                        sample(c("a", "b"), replace = TRUE, nsamples(tmpPhy))),
              "list")
})

test_that("permanova throws errros for wrong input", {
  expect_error(permanova(tmpPhy, get_variable(tmpPy, "Diagnosis")))
  expect_error(permanova(rcmFit, groups = c("a", "b")))
  expect_error(permanova(rcmFit, groups = rep("a", nsamples(tmpPhy))))
  expect_error(permanova(rcmFit, groups = c("b", rep("a", nsamples(tmpPhy)-1))))
})

test_that("permanova trhows warning for low number of permutations", {
  expect_warning(permanova(rcmFit, "Diagnosis", nPerm = 20))
})

#Introduce some NAs
tmpPhyNA = transform_sample_counts(tmpPhy, fun = function(x){
    x[sample(length(x), size = 3)] = NA
    x
})
misUnconstr <- suppressWarnings(RCM(tmpPhyNA, k = 2, allowMissingness = TRUE))
suppressWarnings(misConstrLin <- RCM(tmpPhyNA, k = 2, allowMissingness = TRUE,
                    covariates = c("Diagnosis", "Country", "Gender", "BMI"),
                    confounders = "Age"))
suppressWarnings(misConstrNP <- RCM(tmpPhyNA, k = 2, allowMissingness = TRUE,
                   covariates = c("Diagnosis", "Country", "Gender", "BMI"),
                   confounders = "Age", responseFun = "nonparametric"))
test_that("RCM allows for missingness", {
  expect_is(misUnconstr, "RCM")
  expect_is(misConstrLin, "RCM")
  expect_is(misConstrNP, "RCM")
})
test_that("All plotting functions still work with missing data", {
  expect_silent(plot(misUnconstr))
  expect_silent(plot(misConstrLin))
  expect_silent(plot(misConstrNP))
  expect_silent(plot(misUnconstr, samColour = "Deviance"))
  expect_silent(plot(misUnconstr, taxCol = "Deviance", plotType = "species"))
  expect_silent(plot(misUnconstr, inflVar = "psi"))
  expect_silent(plot(misConstrLin, inflVar = "BMI"))
  expect_warning(plotRespFun(misConstrNP))
})

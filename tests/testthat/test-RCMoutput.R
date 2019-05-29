context("RCM output")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(150)],
                    prune_samples(sample_names(Zeller)[seq_len(100)], Zeller))

test_that("RCM returns element of class RCM", {
  expect_is(RCM(tmpPhy, k = 1), "RCM")
})

test_that("RCM returns phyloseq object", {
  expect_is(RCM(tmpPhy, k = 1)$physeq, "phyloseq")
})

test_that("RCM throws warning when not converged", {
  expect_warning(RCM(tmpPhy, k = 1, maxItOut = 2L))
})

#Introduce some NAs
tmpPhyNA = transform_sample_counts(tmpPhy, fun = function(x){
    x[sample(length(x), size = 3)] = NA
    x
})

test_that("RCM allows for missingness", {
  expect_is(RCM(tmpPhyNA, k = 2, allowMissingness = TRUE), "RCM")
  expect_is(RCM(tmpPhyNA, k = 2, allowMissingness = TRUE,
                covariates = c("Diagnosis", "Country", "Gender", "BMI"),
                confounders = "Age"), "RCM")
  expect_is(RCM(tmpPhyNA, k = 2, allowMissingness = TRUE,
                covariates = c("Diagnosis", "Country", "Gender", "BMI"),
                confounders = "Age", responseFun = "nonparametric"), "RCM")
})

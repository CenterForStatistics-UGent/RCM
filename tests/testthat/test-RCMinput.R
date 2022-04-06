context("RCM input")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(100)],
 prune_samples(sample_names(Zeller)[seq_len(50)], Zeller))

test_that("RCM throws warnings for integer variables", {
 expect_warning(RCM(tmpPhy, covariates = c("Age","Diagnosis"), k = 1))
})

test_that("RCM throws errors for wrong input type", {
  expect_error(RCM("Zeller", covariates = c("Diagnosis", "Country", "Gender"), k = 1))
})

test_that("RCM throws errors when only one covariate with one level supplied", {
  expect_error(suppressWarnings(RCM(Zeller, covariates = "Age", k = 1)))
})

test_that("RCM throws warning when less covariate combinations than samples supplied", {
  expect_warning(RCM(Zeller, covariates = c("Diagnosis", "Country"), k = 1))
})

test_that("RCM throws errors when NAs present in data matrix", {
  expect_error(RCM(matrix(c(1,2,3,NA),2,2), covariates = "Age", k = 1))
})

sample_data(Zeller)$bogusVariable = get_variable(Zeller, "Age") +
  as.integer(get_variable(Zeller,"Gender"))

test_that("RCM throws errors when covariates are aliased", {
    expect_error(suppressWarnings(RCM(Zeller, covariates = c("Diagnosis", "Country",
                                          "Gender","Age", "bogusVariable"),
                   k = 2)))
})

test_that("RCM throws errors when confounders are aliased", {
  expect_error(RCM(Zeller, confounders = c("Diagnosis", "Country",
                                          "Gender","Age", "bogusVariable"),
                   k = 2))
})

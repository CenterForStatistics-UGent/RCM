context("RCM input")

test_that("RCM throws warnings for integer variables", {
 expect_warning(RCM(Zeller, covariates = c("Age","Diagnosis"), k = 1))
})

test_that("RCM throws errors for wrong input type", {
  expect_error(RCM("Zeller", covariates = c("Diagnosis", "Country", "Gender"), k = 1))
})

test_that("RCM throws errors when only one covariate with one level supplied", {
  expect_error(RCM(Zeller, covariates = "Age"), k = 1)
})

context("RCM input")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(100)],
 prune_samples(sample_names(Zeller)[seq_len(50)], Zeller))

test_that("RCM throws warnings for integer variables", {
 expect_warning(RCM(tmpPhy, covariates = c("Age","Diagnosis"), k = 1))
})

test_that("RCM throws errors for wrong input type", {
<<<<<<< HEAD
  expect_error(RCM(otu_table(tmpPhy), covariates = c("Diagnosis", "Country", "Gender"), k = 1))
=======
  expect_error(RCM("Zeller", covariates = c("Diagnosis", "Country", "Gender"), k = 1))
})

test_that("RCM throws errors when only one covariate with one level supplied", {
  expect_error(RCM(Zeller, covariates = "Age"), k = 1)
>>>>>>> AddS4
})

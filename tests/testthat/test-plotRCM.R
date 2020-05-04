context("RCM plot")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(150)],
                    prune_samples(sample_names(Zeller)[seq_len(100)], Zeller))
unconstrRcm = RCM(tmpPhy, k = 2, allowMissingness = TRUE, confounders = "Age")
constrRcm = RCM(tmpPhy, k = 2, allowMissingness = TRUE,
    covariates = c("Diagnosis", "Country", "Gender", "BMI"),
    confounders = "Age")
#
test_that("Inappropriate plotting commands throw errors", {
  expect_error(plot(unconstrRcm, samColour = "Influence"))
  expect_error(plot(constrRcm, inflVar = "Diagnosis", samColour = "Influence"))
})

test_that("Correct plotting commands yield ggplot2 objects", {
  expect_s3_class(plot(unconstrRcm, samColour = "Influence", inflVar = "psi"),
                  "ggplot")
  expect_s3_class(plot(constrRcm, inflVar = "DiagnosisNormal", samColour = "Influence",
                       samShape = "Diagnosis"),
                  "ggplot")
})

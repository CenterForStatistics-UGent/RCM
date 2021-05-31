context("residual plot")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(150)],
                    prune_samples(sample_names(Zeller)[seq_len(100)], Zeller))
unconstrRcm = RCM(tmpPhy, k = 2, allowMissingness = TRUE, confounders = "Age")
constrRcm = RCM(tmpPhy, k = 2, allowMissingness = TRUE,
                covariates = c("Diagnosis", "Country", "Gender", "BMI"),
                confounders = "Age")
test_that("Residual plotting works, also for higher dimensions", {
  expect_silent(residualPlot(constrRcm))
  expect_silent(residualPlot(constrRcm, Dim = 2))
})
test_that("Residual plotting returns errors for uncosntrained ordination", {
  expect_error(residualPlot(unconstrRcm))
})


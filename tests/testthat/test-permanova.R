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

test_that("permanova throws warning for low number of permutations", {
  expect_warning(permanova(rcmFit, "Diagnosis", nPerm = 20))
})

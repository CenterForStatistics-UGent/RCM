context("RCM output")

tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(100)],
                    prune_samples(sample_names(Zeller)[seq_len(50)], Zeller))

test_that("RCM returns element of class RCM", {
  expect_is(RCM(tmpPhy, k = 1), "RCM")
})

test_that("RCM returns phyloseq object", {
  expect_is(RCM(tmpPhy, k = 1)$physeq, "phyloseq")
})

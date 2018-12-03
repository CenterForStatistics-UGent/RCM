context("RCM output")

test_that("RCM returns element of class RCM", {
  expect_is(RCM(Zeller, k = 1), "RCM")
})

test_that("RCM returns phyloseq object", {
  expect_is(RCM(Zeller, k = 1)$physeq, "phyloseq")
})

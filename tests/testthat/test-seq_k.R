context("Index vector of lagrange multipliers of dimension k")

test_that("Number of lagrange multipliers are correct", {
  expect_equal(seq_k(1), seq_len(2))
  expect_equal(seq_k(2), 3:5)
  expect_equal(seq_k(3), 6:9)
})

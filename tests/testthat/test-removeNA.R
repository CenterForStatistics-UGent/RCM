context("Correct for missingness")

n = 50; p = 100
X = matrix(rnbinom(n*p, size = 1, mu = 2), n,p)
X[sample(seq_len(n*p), n)] = NA
mu = outer(rowSums(X, na.rm = TRUE), colSums(X, na.rm = TRUE))/
  sum(X, na.rm = TRUE)

test_that("correctXMissingness() removes all NAs", {
  expect_true(anyNA(RCM:::correctXMissingness(X, mu, FALSE)))
  expect_false(anyNA(RCM:::correctXMissingness(X, mu, TRUE)))
})

n = 20; p = 1000
X = matrix(rnbinom(n*p, size = 1, mu = 2), n,p)
test_that("RCM() deals with high-dimensions", {
  expect_silent(RCM(X))
})

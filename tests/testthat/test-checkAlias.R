context("Alias structure")

# A dataframe with alias structure
df = data.frame(foo = rnorm(10), baa = rep(c(TRUE, FALSE), each = 5),
                foo2 = factor(rep(c("male", "female"), each = 5)))

test_that("Alias structure is discovered", {
  expect_error(checkAlias(df, names(df)))
})
test_that("Alias structure is not falsely discovered", {
  expect_silent(checkAlias(df, c("foo", "baa")))
})

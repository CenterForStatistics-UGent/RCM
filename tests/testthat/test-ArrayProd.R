x = array(c(1,2,5,3,4,9,5,2,5,4,5,6,4,1,1,2,1,3,5,1,4,5,6,9), dim = c(2,3,4))
y = c(1,1,2,4)

context("Array product")

test_that("Fast array product yields correct outcome", {
  expect_equal(
  RCM:::arrayprod(x = x, y = y),
apply(x, c(1,2), function(z, y){sum(z*y)}, y = y))
}
)

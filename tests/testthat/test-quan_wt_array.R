test_that("quan_wt_array matches apply over third dimension (array input)", {
  set.seed(3)
  a <- array(rnorm(3 * 4 * 5), c(3, 4, 5))
  expect_equal(quan_wt_array(a), apply(a, 1:2, stats::median))

  probs <- c(0.25, 0.75)
  out <- quan_wt_array(a, probs = probs)
  expect_equal(dim(out), c(3, 4, 2))
  expect_equal(out[, , 1], apply(a, 1:2, stats::quantile, probs = 0.25, names = FALSE))
  expect_equal(out[, , 2], apply(a, 1:2, stats::quantile, probs = 0.75, names = FALSE))
})

test_that("quan_wt_array accepts lists of matrices and weights", {
  set.seed(4)
  l <- lapply(1:5, function(i) matrix(rnorm(6), 2, 3))
  a <- array(unlist(l), c(2, 3, 5))
  expect_equal(quan_wt_array(l), quan_wt_array(a))

  wt <- array(1, c(2, 3, 5))
  expect_equal(quan_wt_array(a, wt = wt), quan_wt_array(a))
  expect_equal(quan_wt_array(l, wt = lapply(1:5, function(i) matrix(1, 2, 3))),
               quan_wt_array(a))
})

test_that("quan_wt_array validates inputs", {
  bad <- list(matrix(1, 2, 2), matrix(1, 3, 2))
  expect_error(quan_wt_array(bad), "does not match")
  expect_error(quan_wt_array(list(1:3, 1:3)), "not a matrix")
  a <- array(1, c(2, 2, 3))
  expect_error(quan_wt_array(a, wt = array(1, c(2, 2, 2))), "do not match")
})

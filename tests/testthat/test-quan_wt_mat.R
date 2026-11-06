test_that("quan_wt_mat_col and quan_wt_mat_row match base quantile", {
  set.seed(1)
  m <- matrix(rnorm(40), 8, 5)
  expect_equal(quan_wt_mat_col(m), apply(m, 2, stats::median))
  expect_equal(quan_wt_mat_row(m), apply(m, 1, stats::median))

  probs <- c(0.1, 0.5, 0.9)
  col <- quan_wt_mat_col(m, probs)
  row <- quan_wt_mat_row(m, probs)
  expect_equal(dim(col), c(5, 3))
  expect_equal(dim(row), c(8, 3))
  expect_equal(col, t(apply(m, 2, stats::quantile, probs = probs, names = FALSE)))
  expect_equal(row, t(apply(m, 1, stats::quantile, probs = probs, names = FALSE)))
})

test_that("quan_wt_mat_* with equal weights match unweighted, and bad weight dims error", {
  set.seed(2)
  m <- matrix(rnorm(24), 6, 4)
  w <- matrix(1, 6, 4)
  expect_equal(quan_wt_mat_col(m, wt = w), quan_wt_mat_col(m))
  expect_equal(quan_wt_mat_row(m, wt = w), quan_wt_mat_row(m))
  expect_error(quan_wt_mat_col(m, wt = matrix(1, 5, 4)))
  expect_error(quan_wt_mat_row(m, wt = matrix(1, 6, 3)))
})

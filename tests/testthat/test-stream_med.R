test_that("streaming median matches stats::median", {
  clearHeaps()
  on.exit(clearHeaps())
  expect_true(is.na(getMedian()))

  set.seed(5)
  x <- rnorm(101)
  addVector(x)
  expect_equal(getMedian(), stats::median(x))

  y <- rnorm(20)
  addVector(y)
  expect_equal(getMedian(), stats::median(c(x, y)))
})

test_that("streaming median ignores NA and clearHeaps resets", {
  clearHeaps()
  on.exit(clearHeaps())
  addVector(c(1, NA, 3, 2, NA))
  expect_equal(getMedian(), 2)
  expect_equal(length(getMaxHeap()) + length(getMinHeap()), 3)
  clearHeaps()
  expect_equal(length(getMaxHeap()), 0)
  expect_equal(length(getMinHeap()), 0)
  expect_true(is.na(getMedian()))
})

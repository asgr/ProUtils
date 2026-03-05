## Tests for quan_wt (weighted quantiles), covering types 1-9.

# Helper: equal-weight vector of length n
eq_wt <- function(n) rep(1, n)

# ── Equal-weight collapse to base R ──────────────────────────────────────────

test_that("quan_wt with equal weights matches stats::quantile for types 1-9", {
  set.seed(42)
  x <- rnorm(20)
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  wt    <- eq_wt(length(x))

  for (t in 1:9) {
    expected <- stats::quantile(x, probs = probs, type = t, names = FALSE)
    got      <- quan_wt(x, probs = probs, wt = wt, type = t)
    expect_equal(got, expected, tolerance = 1e-10,
                 label = paste0("type=", t))
  }
})

test_that("quan_wt equal-weight type 1/2/3 matches base R on simple integer sequence", {
  # Use 1:4 so CDF steps land exactly at 0.25, 0.5, 0.75, 1.0 for clean testing
  x  <- 1:4
  wt <- eq_wt(4)
  # Probe values including the exact CDF step-boundaries (for type 2 averaging)
  probs <- c(0, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.9, 1)

  for (t in 1:3) {
    expected <- stats::quantile(x, probs = probs, type = t, names = FALSE)
    got      <- quan_wt(x, probs = probs, wt = wt, type = t)
    expect_equal(got, expected, tolerance = 1e-10,
                 label = paste0("type=", t, " on 1:4"))
  }
})

# ── Weighted case: monotonicity ───────────────────────────────────────────────

test_that("quan_wt types 1-9 produce monotone results with unequal weights", {
  set.seed(7)
  x  <- rnorm(30)
  wt <- runif(30, 0.5, 2)   # positive, unequal weights
  probs <- seq(0, 1, by = 0.1)

  for (t in 1:9) {
    q <- quan_wt(x, probs = probs, wt = wt, type = t)
    expect_true(all(diff(q) >= 0), label = paste0("monotone for type=", t))
  }
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("quan_wt returns NA for empty input", {
  for (t in 1:9) {
    expect_true(is.na(quan_wt(numeric(0), probs = 0.5, wt = numeric(0), type = t)),
                label = paste0("type=", t))
  }
})

test_that("quan_wt returns the single value for length-1 input", {
  for (t in 1:9) {
    expect_equal(quan_wt(42, probs = c(0, 0.5, 1), wt = 1, type = t),
                 c(42, 42, 42),
                 label = paste0("type=", t))
  }
})

test_that("quan_wt returns NA when all weights are zero (weighted types 1-3)", {
  x  <- c(1, 2, 3)
  wt <- c(0, 0, 0)
  for (t in 1:3) {
    expect_true(is.na(quan_wt(x, probs = 0.5, wt = wt, type = t)),
                label = paste0("type=", t))
  }
})

test_that("quan_wt ignores NA values in x or wt", {
  x_na  <- c(1, NA, 3, 4)
  wt_na <- c(1, 1, NA, 1)
  # After filtering: x = c(1, 4), wt = c(1, 1) → median = 2.5
  expect_equal(quan_wt(x_na, probs = 0.5, wt = wt_na, type = 7), 2.5)
})

test_that("quan_wt rejects negative weights", {
  expect_error(quan_wt(1:3, probs = 0.5, wt = c(1, -1, 1)),
               regexp = "non-negative")
})

test_that("quan_wt rejects type outside 1-9", {
  expect_error(quan_wt(1:5, probs = 0.5, wt = eq_wt(5), type = 0),
               regexp = "type")
  expect_error(quan_wt(1:5, probs = 0.5, wt = eq_wt(5), type = 10),
               regexp = "type")
})

# ── No-weight path still supports types 1-3 ──────────────────────────────────

test_that("quan_wt without weights (NULL) supports types 1-3 via stats::quantile", {
  set.seed(9)
  x <- rnorm(20)
  probs <- c(0.1, 0.5, 0.9)
  for (t in 1:3) {
    expected <- stats::quantile(x, probs = probs, type = t, names = FALSE)
    got      <- quan_wt(x, probs = probs, type = t)
    expect_equal(got, expected, tolerance = 1e-10,
                 label = paste0("NULL wt, type=", t))
  }
})

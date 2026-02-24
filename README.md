
# ProUtils

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/ProUtils/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

**ProUtils** is an R package providing utility functions for efficient data processing and statistical operations, including weighted quantiles, running medians, and matrix-based computations. It leverages **Rcpp** for high-performance C++ implementations.

---

## ✨ Features

- **Weighted Quantiles**  
  Compute quantiles for vectors and matrices with optional weights:
  - `quan_wt()` – weighted quantiles for a numeric vector.
  - `quan_wt_mat_col()` – column-wise weighted quantiles for a matrix.
  - `quan_wt_mat_row()` – row-wise weighted quantiles for a matrix.

- **Running Median Tracker**  
  Maintain a running median using two heaps:
  - `addVector()` – add values incrementally.
  - `getMedian()` – retrieve the current median.
  - `clearHeaps()` – reset the state.
  - `getMinHeap()` / `getMaxHeap()` – inspect heap contents.

- **High Performance**  
  All core functions are implemented in C++ via **Rcpp** for speed.

---

## 📦 Installation

You can install the development version from GitHub using **devtools**:

```r
# Install devtools if needed
install.packages("remotes")

# Install ProUtils from GitHub
remotes::install_github("asgr/ProUtils")


library(ProUtils)
```

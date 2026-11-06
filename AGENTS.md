# AGENTS.md

ProUtils: low-level utility functions shared across ProTools packages. R package with Rcpp (C++) code.

## Layout
- `/src/quan_wt.cpp`: weighted quantiles (types 1-9): `quan_wt`, `quan_wt_mat_col`, `quan_wt_mat_row`.
- `/src/stream_med.cpp`: streaming median via global heaps (`addVector`, `getMedian`, `clearHeaps`, `getMinHeap`, `getMaxHeap`). State is global; tests must call `clearHeaps()`.
- `/R/quan_wt_array.R`: `quan_wt_array` (array/list of matrices -> quantiles over 3rd dim).
- `/R/RcppExports.R`, `/src/RcppExports.cpp`: generated; after changing `// [[Rcpp::export]]` signatures run `Rscript -e 'Rcpp::compileAttributes()'`. Do not edit by hand.
- `/NAMESPACE` uses `exportPattern`, so every alphabetically-named R object is exported.
- `/tests/testthat/`: testthat 3rd edition tests (`test-*.R`).

## Commands
- Install: `R CMD INSTALL .`
- Test: `Rscript -e 'devtools::test()'` (or `R CMD check --as-cran --no-manual .`)
- Docs: `Rscript -e 'roxygen2::roxygenise()'` only if roxygen is introduced; man pages in `/man` are otherwise hand-maintained.

## Conventions
- Equal-weight results must agree with `stats::quantile` for each type; keep tests for this.
- Update `Version`/`Date` in `/DESCRIPTION` for user-visible changes.
- CI: `/.github/workflows/main.yml` (R-CMD-check on multiple OSes) and `/.github/workflows/coverage.yml` (covr).

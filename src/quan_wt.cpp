#include <Rcpp.h>
// #ifdef _OPENMP
// #include <omp.h>
// #endif

using namespace Rcpp;

NumericVector _quan_wt(const NumericVector &x, const NumericVector probs = NumericVector::create(0.5), const Nullable<NumericVector> wt = R_NilValue, const int type = 7) {
  int n = x.size();
  int nprobs = probs.size();

  // Special return case where x is empty/NULL/missing:
  if (n == 0) return NumericVector(nprobs, NA_REAL);

  if (type < 1 || type > 9) stop("type must be between 1 and 9 (inclusive)");

  // Validate q values
  for (int i = 0; i < nprobs; i++) {
    if (probs[i] < 0.0 || probs[i] > 1.0) stop("All prob values must be between 0 and 1");
  }

  NumericVector wt_vec;
  if (wt.isNull()) {
    //wt_vec = NumericVector(n, 1.0); // equal weights if none provided
    Rcpp::Environment stats_env = Rcpp::Environment::namespace_env("stats");
    if (type == 7 && probs.size() == 1 && std::abs(probs[0] - 0.5) < 1e-12) {
      Rcpp::Function my_median = stats_env["median"];
      return my_median(x, true); // positional argument for na.rm
    } else {
      Rcpp::Function my_quantile = stats_env["quantile"];
      return my_quantile(x,
                         Named("probs") = probs,
                         Named("na.rm") = true,
                         Named("names") = false,
                         Named("type") = type);
    }
  } else {
    wt_vec = NumericVector(wt);
    if (wt_vec.size() != n) stop("x and wt must have the same length");
    for (int i = 0; i < n; i++) {
      if (!NumericVector::is_na(wt_vec[i]) && wt_vec[i] < 0.0) stop("All wt values must be non-negative");
    }
  }

  // Prepare data
  std::vector<std::pair<double, double>> data;
  data.reserve(n);

  // Filter out NA
  for (int i = 0; i < n; i++) {
    if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(wt_vec[i])) {
      data.push_back(std::make_pair(x[i], wt_vec[i]));
    }
  }

  // Special return case where everything is NA:
  if (data.empty()) return NumericVector(nprobs, NA_REAL);
  // Special return case where remaining length is 1 (so answer should always be that, no matter the quantile)
  if (data.size() == 1) return NumericVector(nprobs, data[0].first);

  // Otherwise we have more values to deal with...
  // Sort by x
  std::sort(data.begin(), data.end(),
            [](const std::pair<double,double>& a, const std::pair<double,double>& b) {
              return a.first < b.first;
            });

  // Types 1-3: step-function (empirical CDF) based
  if (type <= 3) {
    // Build right-continuous CDF: cdf_right[i] = sum(w[0..i]) / total_weight
    double total_weight = 0.0;
    for (size_t i = 0; i < data.size(); i++) total_weight += data[i].second;

    // Special return case where there is no weight at all:
    if (total_weight == 0.0) return NumericVector(nprobs, NA_REAL);

    std::vector<double> cdf_right(data.size());
    double cumw = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
      cumw += data[i].second;
      cdf_right[i] = cumw / total_weight;
    }

    // For type 3, precompute midpoints: mid[i] = (cdf_left[i] + cdf_right[i]) / 2
    std::vector<double> mid;
    if (type == 3) {
      mid.resize(data.size());
      for (size_t i = 0; i < data.size(); i++) {
        double cdf_left = (i > 0) ? cdf_right[i - 1] : 0.0;
        mid[i] = (cdf_left + cdf_right[i]) / 2.0;
      }
    }

    NumericVector result(nprobs);
    for (int k = 0; k < nprobs; k++) {
      double prob = probs[k];
      if (prob <= 0.0) { result[k] = data.front().first; continue; }
      if (prob >= 1.0) { result[k] = data.back().first;  continue; }

      if (type == 1) {
        // Smallest i such that cdf_right[i] >= prob
        auto it = std::lower_bound(cdf_right.begin(), cdf_right.end(), prob);
        size_t i = static_cast<size_t>(it - cdf_right.begin());
        result[k] = data[std::min(i, data.size() - 1)].first;

      } else if (type == 2) {
        // Like type 1, but average adjacent values when prob falls exactly on a CDF step
        auto it = std::lower_bound(cdf_right.begin(), cdf_right.end(), prob);
        size_t i = static_cast<size_t>(it - cdf_right.begin());
        if (i >= data.size()) i = data.size() - 1;
        if (cdf_right[i] == prob && i + 1 < data.size()) {
          result[k] = (data[i].first + data[i + 1].first) / 2.0;
        } else {
          result[k] = data[i].first;
        }

      } else {
        // type == 3: nearest order statistic, with R-compatible tie-breaking.
        //
        // The region for x[i] (0-indexed) is (mid[i], mid[i+1]) strictly.
        // At an exact midpoint mid[j]:
        //   1-indexed j+1 is odd  → assign x[j]   (0-indexed)
        //   1-indexed j+1 is even → assign x[j-1] (0-indexed)
        // For p < mid[0]: assign x[0]; for p > mid[n-1]: assign x[n-1].

        // Find first mid[j] >= prob
        auto it = std::lower_bound(mid.begin(), mid.end(), prob);
        size_t j = static_cast<size_t>(it - mid.begin());

        if (j >= data.size()) {
          result[k] = data.back().first;
        } else if (mid[j] == prob) {
          // Exact midpoint: apply round-half-to-even (by 1-indexed position)
          if ((j + 1) % 2 == 1) {  // 1-indexed j+1 is odd
            result[k] = data[j].first;
          } else {                  // 1-indexed j+1 is even
            result[k] = (j > 0) ? data[j - 1].first : data[0].first;
          }
        } else if (j == 0) {
          // prob < mid[0]: below first midpoint
          result[k] = data[0].first;
        } else {
          // prob in (mid[j-1], mid[j]) strictly → assign x[j-1]
          result[k] = data[j - 1].first;
        }
      }
    }
    return result;
  }

  // Types 4-9: trapezoid/midpoint grid with linear interpolation
  // Adjust quantile position based on type
  double start_mod;
  double end_mod;
  switch (type) {
    case 4: start_mod = 1; end_mod = 0; break;
    case 5: start_mod = 0.5; end_mod = 0.5; break;
    case 6: start_mod = 1; end_mod = 1; break;
    case 7: start_mod = 0; end_mod = 0; break;
    case 8: start_mod = 2.0/3.0; end_mod = 2.0/3.0; break;
    case 9: start_mod = 0.625; end_mod = 0.625; break;
    default: stop("type must be between 1 and 9 (inclusive)"); break;
  }

  // Compute total weight

  // Compute cumulative probabilities using midpoint rule
  std::vector<double> cumProb(data.size());
  double tempWeight = data[0].second * start_mod; //start with a half bin
  cumProb[0] = tempWeight;

  for (size_t i = 1; i < data.size(); i++) {
    tempWeight += (data[i - 1].second + data[i].second) / 2.0;
    cumProb[i] = tempWeight;
  }

  // So we count the half bin at the high end correctly
  double totalWeight = cumProb[data.size() - 1] + data[data.size() - 1].second * end_mod;

  // Special return case where there is no weight at all:
  if (totalWeight == 0.0) return NumericVector(nprobs, NA_REAL);

  for (size_t i = 0; i < data.size(); i++) {
    cumProb[i] /= totalWeight;
  }

  // Compute quantiles for our non-trivial case:
  NumericVector result(nprobs);
  for (int k = 0; k < nprobs; k++) {
    double prob = probs[k];
    if (prob <= 0.0) {
      result[k] = data.front().first;
      continue;
    }
    if (prob >= 1.0) {
      result[k] = data.back().first;
      continue;
    }

    // Use binary search for O(log n) lookup instead of linear scan
    auto it = std::lower_bound(cumProb.begin(), cumProb.end(), prob);
    if (it == cumProb.end()) {
      result[k] = data.back().first;
    } else {
      size_t i = static_cast<size_t>(it - cumProb.begin());
      if (i == 0) {
        result[k] = data[i].first;
      } else {
        double prevVal = data[i - 1].first;
        double currVal = data[i].first;
        double prevProb = cumProb[i - 1];
        double currProb = cumProb[i];
        if (currProb == prevProb) {
          result[k] = (prevVal + currVal) / 2.0;
        } else {
          double fraction = (prob - prevProb) / (currProb - prevProb);
          result[k] = prevVal + fraction * (currVal - prevVal);
        }
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector quan_wt(NumericVector x, NumericVector probs = NumericVector::create(0.5),
                      Nullable<NumericVector> wt = R_NilValue, int type = 7){
  // Rcpp::Environment collapse_env = Rcpp::Environment::namespace_env("collapse");
  // Rcpp::Function my_fquantile = collapse_env["fquantile"];
  // return my_fquantile(x,
  //                      Named("probs") = probs,
  //                      Named("w") = wt,
  //                      Named("na.rm") = true,
  //                      Named("names") = false,
  //                      Named("type") = type);
 return _quan_wt(x, probs, wt, type);
}

// [[Rcpp::export]]
SEXP quan_wt_mat_col(NumericMatrix mat, NumericVector probs = NumericVector::create(0.5),
                              Nullable<NumericMatrix> wt = R_NilValue, int type = 7) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nprobs = probs.size();

  // Rcpp::Environment collapse_env = Rcpp::Environment::namespace_env("collapse");
  // Rcpp::Function my_fquantile = collapse_env["fquantile"];

  // Result matrix: rows = columns of mat, columns = quantiles
  NumericMatrix result(ncol, nprobs);

  if (wt.isNull()) {
    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int j = 0; j < ncol; j++) {
      NumericVector colQuant = _quan_wt(mat(_, j), probs, R_NilValue, type);
      // NumericVector colQuant = my_fquantile(mat(_, j),
      //              Named("probs") = probs,
      //              Named("na.rm") = true,
      //              Named("names") = false,
      //              Named("type") = type);
      for (int k = 0; k < nprobs; k++) {
        result(j, k) = colQuant[k];
      }
    }
  } else {
    NumericMatrix wt_temp = as<NumericMatrix>(wt);
    if (wt_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (wt_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int j = 0; j < ncol; j++) {
      NumericVector wt_vec = wt_temp(_, j);
      NumericVector colQuant = _quan_wt(mat(_, j), probs, wt_vec, type);
      // NumericVector colQuant = my_fquantile(mat(_, j),
      //                                       Named("probs") = probs,
      //                                       Named("w") = wt_vec,
      //                                       Named("na.rm") = true,
      //                                       Named("names") = false,
      //                                       Named("type") = type);
      for (int k = 0; k < nprobs; k++) {
        result(j, k) = colQuant[k];
      }
    }
  }

  // Convert to vector if only one quantile requested
  if (nprobs == 1) {
    NumericVector vecResult(ncol);
    for (int i = 0; i < ncol; i++) {
      vecResult[i] = result(i, 0);
      }
    return vecResult;
  }else{
    return result;
  }
}

// [[Rcpp::export]]
SEXP quan_wt_mat_row(NumericMatrix mat, NumericVector probs = NumericVector::create(0.5),
                              Nullable<NumericMatrix> wt = R_NilValue, int type = 7) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nprobs = probs.size();

  // Rcpp::Environment collapse_env = Rcpp::Environment::namespace_env("collapse");
  // Rcpp::Function my_fquantile = collapse_env["fquantile"];

  // Result matrix: rows = rows of mat, columns = quantiles
  NumericMatrix result(nrow, nprobs);

  if (wt.isNull()) {
    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int i = 0; i < nrow; i++) {
      NumericVector rowQuant = _quan_wt(mat(i, _), probs, R_NilValue, type);
      // NumericVector rowQuant = my_fquantile(mat(i, _),
      //                                       Named("probs") = probs,
      //                                       Named("na.rm") = true,
      //                                       Named("names") = false,
      //                                       Named("type") = type);
      for (int k = 0; k < nprobs; k++) {
        result(i, k) = rowQuant[k];
      }
    }
  } else {
    NumericMatrix wt_temp = as<NumericMatrix>(wt);
    if (wt_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (wt_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int i = 0; i < nrow; i++) {
      NumericVector wt_vec = wt_temp(i, _);
      NumericVector rowQuant = _quan_wt(mat(i, _), probs, wt_vec, type);
      // NumericVector rowQuant = my_fquantile(mat(i, _),
      //                                       Named("probs") = probs,
      //                                       Named("w") = wt_vec,
      //                                       Named("na.rm") = true,
      //                                       Named("names") = false,
      //                                       Named("type") = type);
      for (int k = 0; k < nprobs; k++) {
        result(i, k) = rowQuant[k];
      }
    }
  }

  // Convert to vector if only one quantile requested
  if (nprobs == 1) {
    NumericVector vecResult(nrow);
    for (int i = 0; i < nrow; i++) {
      vecResult[i] = result(i, 0);
    }
    return vecResult;
  }else{
    return result;
  }
}

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

  if (type < 4 || type > 9) stop("type must be between 5 and 9 (inclusive)");

  // Validate q values
  for (int i = 0; i < nprobs; i++) {
    if (probs[i] < 0.0 || probs[i] > 1.0) stop("All prob values must be between 0 and 1");
  }

  // Prepare data
  std::vector<std::pair<double, double>> data;
  data.reserve(n);

  NumericVector wt_vec;
  if (wt.isNull()) {
    //wt_vec = NumericVector(n, 1.0); // equal weights if none provided
    Rcpp::Environment stats_env = Rcpp::Environment::namespace_env("stats");
    if (probs.size() == 1 && std::abs(probs[0] - 0.5) < 1e-12) {
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
  }

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

  // Adjust quantile position based on type
  double start_mod;
  double end_mod;
  switch (type) {
    case 4: start_mod = 1; end_mod = 0; break;
    case 5: start_mod = 0.5; end_mod = 0.5; break;
    case 6: start_mod = 1; end_mod = 1; break;
    case 7: start_mod = 0; end_mod = 0; break;
    case 8: start_mod = 0.6666667; end_mod = 0.6666667; break;
    case 9: start_mod = 0.625; end_mod = 0.625; break;
  }

  // Compute total weight

  // Compute cumulative probabilities using midpoint rule
  std::vector<double> cumProb(data.size());
  double tempWeight = data[0].second * start_mod; //start with a half bin
  cumProb[0] = tempWeight;

  for (size_t i = 1; i < data.size(); i++) {
    tempWeight += data[i - 1].second / 2.0 + data[i].second / 2.0;
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

    result[k] = data.back().first; // default
    for (size_t i = 0; i < data.size(); i++) {
      if (cumProb[i] >= prob) {
        if (i == 0) {
          result[k] = data[i].first;
        } else {
          double prevVal = data[i - 1].first;
          double currVal = data[i].first;
          double prevProb = cumProb[i - 1];
          double currProb = cumProb[i];
          double fraction = (prob - prevProb) / (currProb - prevProb);
          result[k] = prevVal + fraction * (currVal - prevVal);
        }
        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector quan_wt(NumericVector x, NumericVector probs = NumericVector::create(0.5),
                      Nullable<NumericVector> wt = R_NilValue, int type = 7){
  return _quan_wt(x, probs, wt, type);
}

// [[Rcpp::export]]
SEXP quan_wt_mat_col(NumericMatrix mat, NumericVector probs = NumericVector::create(0.5),
                              Nullable<NumericMatrix> wt = R_NilValue, int type = 7) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nprobs = probs.size();

  // Result matrix: rows = columns of mat, columns = quantiles
  NumericMatrix result(ncol, nprobs);

  if (wt.isNull()) {
    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int j = 0; j < ncol; j++) {
      NumericVector colQuant = _quan_wt(mat(_, j), probs);
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

  // Result matrix: rows = rows of mat, columns = quantiles
  NumericMatrix result(nrow, nprobs);

  if (wt.isNull()) {
    // #ifdef _OPENMP
    //     // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
    // #pragma omp parallel for schedule(dynamic, 10) if(nrow > 100) num_threads(nthreads)
    // #endif
    for (int i = 0; i < nrow; i++) {
      NumericVector rowQuant = _quan_wt(mat(i, _), probs);
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

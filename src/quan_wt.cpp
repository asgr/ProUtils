#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector quan_wt(NumericVector x, NumericVector probs = NumericVector::create(0.5),
                      Nullable<NumericVector> wt = R_NilValue, int type = 7) {
  int n = x.size();
  int nprobs = probs.size();
  
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
    wt_vec = NumericVector(n, 1.0); // equal weights if none provided
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
  if (data.empty()) return NumericVector(nprobs, NA_REAL);

  // Sort by x
  std::sort(data.begin(), data.end(),
            [](const std::pair<double,double>& a, const std::pair<double,double>& b) {
              return a.first < b.first;
            });

  // Adjust plotting position based on type
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

  // so we count the half bin at the high end correctly
  double totalWeight = cumProb[data.size() - 1] + data[data.size() - 1].second * end_mod;

  if (totalWeight == 0.0) return NumericVector(nprobs, NA_REAL);

  for (size_t i = 0; i < data.size(); i++) {
    cumProb[i] /= totalWeight;
  }

  // Compute quantiles
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
SEXP quan_wt_mat_col(NumericMatrix mat, NumericVector probs = NumericVector::create(0.5),
                              Nullable<NumericMatrix> wt = R_NilValue, int type = 7) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nprobs = probs.size();

  // Result matrix: rows = columns of mat, columns = quantiles
  NumericMatrix result(ncol, nprobs);

  if (wt.isNull()) {
    for (int j = 0; j < ncol; j++) {
      NumericVector colQuant = quan_wt(mat(_, j), probs);
      for (int k = 0; k < nprobs; k++) {
        result(j, k) = colQuant[k];
      }
    }
  } else {
    NumericMatrix wt_temp = as<NumericMatrix>(wt);
    if (wt_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (wt_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    for (int j = 0; j < ncol; j++) {
      NumericVector wt_vec = wt_temp(_, j);
      NumericVector colQuant = quan_wt(mat(_, j), probs, wt_vec, type);
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
    for (int i = 0; i < nrow; i++) {
      NumericVector rowQuant = quan_wt(mat(i, _), probs);
      for (int k = 0; k < nprobs; k++) {
        result(i, k) = rowQuant[k];
      }
    }
  } else {
    NumericMatrix wt_temp = as<NumericMatrix>(wt);
    if (wt_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (wt_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    for (int i = 0; i < nrow; i++) {
      NumericVector wt_vec = wt_temp(i, _);
      NumericVector rowQuant = quan_wt(mat(i, _), probs, wt_vec, type);
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


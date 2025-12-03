#include <Rcpp.h>
using namespace Rcpp;


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector quan_wt(NumericVector x, NumericVector q = NumericVector::create(0.5),
                      Nullable<NumericVector> w = R_NilValue) {
  int n = x.size();
  int nq = q.size();

  // Validate q values
  for (int i = 0; i < nq; i++) {
    if (q[i] < 0.0 || q[i] > 1.0) stop("All q values must be between 0 and 1");
  }

  // Prepare data
  std::vector<std::pair<double, double>> data;
  data.reserve(n);

  NumericVector wvec;
  if (w.isNull()) {
    wvec = NumericVector(n, 1.0); // equal weights if none provided
  } else {
    wvec = NumericVector(w);
    if (wvec.size() != n) stop("x and w must have the same length");
  }

  // Filter out NA
  for (int i = 0; i < n; i++) {
    if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(wvec[i])) {
      data.push_back(std::make_pair(x[i], wvec[i]));
    }
  }
  if (data.empty()) return NumericVector(nq, NA_REAL);

  // Sort by x
  std::sort(data.begin(), data.end(),
            [](const std::pair<double,double>& a, const std::pair<double,double>& b) {
              return a.first < b.first;
            });

  // Compute total weight


  // Compute cumulative probabilities using midpoint rule
  std::vector<double> cumProb(data.size());
  double tempWeight = data[0].second / 2.0; //start with a half bin
  cumProb[0] = tempWeight;

  for (size_t i = 1; i < data.size(); i++) {
    tempWeight += data[i - 1].second / 2.0 + data[i].second / 2.0;
    cumProb[i] = tempWeight;
  }

  // so we count the half bin at the high end correctly
  double totalWeight = cumProb[data.size() - 1] + data[data.size() - 1].second / 2.0;

  if (totalWeight == 0.0) return NumericVector(nq, NA_REAL);

  for (size_t i = 0; i < data.size(); i++) {
    //Rcout << cumProb[i] << " ";
    cumProb[i] /= totalWeight;
    //Rcout << cumProb[i] << " ";
  }
  //Rcout << std::endl;

  // Compute quantiles
  NumericVector result(nq);
  for (int k = 0; k < nq; k++) {
    double qq = q[k];
    if (qq <= 0.0) {
      result[k] = data.front().first;
      continue;
    }
    if (qq >= 1.0) {
      result[k] = data.back().first;
      continue;
    }

    result[k] = data.back().first; // default
    for (size_t i = 0; i < data.size(); i++) {
      if (cumProb[i] >= qq) {
        if (i == 0) {
          result[k] = data[i].first;
        } else {
          double prevVal = data[i - 1].first;
          double currVal = data[i].first;
          double prevProb = cumProb[i - 1];
          double currProb = cumProb[i];
          double fraction = (qq - prevProb) / (currProb - prevProb);
          result[k] = prevVal + fraction * (currVal - prevVal);
        }
        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
SEXP quan_wt_mat_col(NumericMatrix mat, NumericVector q = NumericVector::create(0.5),
                              Nullable<NumericMatrix> w = R_NilValue) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nq = q.size();

  // Result matrix: rows = columns of mat, columns = quantiles
  NumericMatrix result(ncol, nq);

  if (w.isNull()) {
    for (int j = 0; j < ncol; j++) {
      NumericVector colQuant = quan_wt(mat(_, j), q);
      for (int k = 0; k < nq; k++) {
        result(j, k) = colQuant[k];
      }
    }
  } else {
    NumericMatrix w_temp = as<NumericMatrix>(w);
    if (w_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (w_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    for (int j = 0; j < ncol; j++) {
      NumericVector wvec = w_temp(_, j);
      NumericVector colQuant = quan_wt(mat(_, j), q, wvec);
      for (int k = 0; k < nq; k++) {
        result(j, k) = colQuant[k];
      }
    }
  }

  // Convert to vector if only one quantile requested
  if (nq == 1) {
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
SEXP quan_wt_mat_row(NumericMatrix mat, NumericVector q = NumericVector::create(0.5),
                              Nullable<NumericMatrix> w = R_NilValue) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int nq = q.size();

  // Result matrix: rows = rows of mat, columns = quantiles
  NumericMatrix result(nrow, nq);

  if (w.isNull()) {
    for (int i = 0; i < nrow; i++) {
      NumericVector rowQuant = quan_wt(mat(i, _), q);
      for (int k = 0; k < nq; k++) {
        result(i, k) = rowQuant[k];
      }
    }
  } else {
    NumericMatrix w_temp = as<NumericMatrix>(w);
    if (w_temp.nrow() != nrow) stop("Row dims of weights must match mat");
    if (w_temp.ncol() != ncol) stop("Column dims of weights must match mat");

    for (int i = 0; i < nrow; i++) {
      NumericVector wvec = w_temp(i, _);
      NumericVector rowQuant = quan_wt(mat(i, _), q, wvec);
      for (int k = 0; k < nq; k++) {
        result(i, k) = rowQuant[k];
      }
    }
  }

  // Convert to vector if only one quantile requested
  if (nq == 1) {
    NumericVector vecResult(nrow);
    for (int i = 0; i < nrow; i++) {
      vecResult[i] = result(i, 0);
    }
    return vecResult;
  }else{
    return result;
  }
}


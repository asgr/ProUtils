#include <Rcpp.h>
#include <queue>
#include <vector>

using namespace Rcpp;

// Define two heaps:
// maxHeap for lower half (largest element on top)
// minHeap for upper half (smallest element on top)
std::priority_queue<double> maxHeap;
std::priority_queue<double, std::vector<double>, std::greater<double>> minHeap;

void addNumber(double num) {
  // Do nothing if num is NA
  if (NumericVector::is_na(num)) return;

  // Insert into appropriate heap
  if (maxHeap.empty() || num <= maxHeap.top()) {
    maxHeap.push(num);
  } else {
    minHeap.push(num);
  }

  // Balance heaps so size difference <= 1
  if (maxHeap.size() > minHeap.size() + 1) {
    minHeap.push(maxHeap.top());
    maxHeap.pop();
  } else if (minHeap.size() > maxHeap.size()) {
    maxHeap.push(minHeap.top());
    minHeap.pop();
  }
}

// [[Rcpp::export]]
void addVector(NumericVector vec){
  for (int i = 0; i < vec.size(); i++) {
    addNumber(vec[i]);
  }
}

// [[Rcpp::export]]
double getMedian() {
  if (maxHeap.empty() && minHeap.empty()) {
    return NA_REAL; // No elements yet
  }
  if (maxHeap.size() == minHeap.size()) {
    return (maxHeap.top() + minHeap.top()) / 2.0;
  } else {
    return maxHeap.top(); // maxHeap will have one extra element
  }
}

// [[Rcpp::export]]
void clearHeaps() {
  maxHeap = std::priority_queue<double>();
  minHeap = std::priority_queue<double, std::vector<double>, std::greater<double>>();
}


// [[Rcpp::export]]
NumericVector getMinHeap() {
  // Copy elements from minHeap into a vector
  std::vector<double> temp;
  auto copyHeap = minHeap; // Make a copy so we don't destroy the original

  while (!copyHeap.empty()) {
    temp.push_back(copyHeap.top());
    copyHeap.pop();
  }

  return NumericVector(temp.begin(), temp.end());
}

// [[Rcpp::export]]
NumericVector getMaxHeap() {
  // Copy elements from minHeap into a vector
  std::vector<double> temp;
  auto copyHeap = maxHeap; // Make a copy so we don't destroy the original

  while (!copyHeap.empty()) {
    temp.push_back(copyHeap.top());
    copyHeap.pop();
  }

  return NumericVector(temp.begin(), temp.end());
}

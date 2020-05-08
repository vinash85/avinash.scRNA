#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


// [[Rcpp::export]]
Eigen::SparseMatrix<double> myLogNorm(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(2*data.outerSize(), display_progress);
  // Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      // std::cout<< it.value()<< std::endl;
      it.valueRef() = log1p(double(it.value()) * scale_factor);
    }
  }
  // this does not works because this iterates over only nonzero elements.
  // Eigen::VectorXd colMeans = data.transpose() * Eigen::VectorXd::Constant(data.rows(), 1.0/data.cols());
  // std::cout<< colMeans << std::endl;
  // for (int k=0; k < data.outerSize(); ++k){
  //   p.increment();
  //   for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
  //     it.valueRef() = it.value() - colMeans[k];
  //   }
  // }
  return data;
}

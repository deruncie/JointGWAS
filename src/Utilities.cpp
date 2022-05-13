#include "JointGWAS_types.h"

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export()]]
MatrixXd matrix_multiply_toDense(SEXP X_, SEXP Y_){
  if(Rf_isNull(X_)) return(as<Map<MatrixXd> >(Y_));
  if(Rf_isMatrix(X_)) {
    Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
  else {
    MSpMat X = as<MSpMat>(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      setNbThreads(0);
      if(Eigen::nbThreads( ) > 1) {
        SparseMatrix<double,RowMajor> Xr = X;  // Convert to RowMajor so it can be parallelized
        return(Xr*Y);
      } else{
        return(X*Y);
      }
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
}

// [[Rcpp::export()]]
MatrixXd partial_matrix_multiply_toDense(SEXP X_, SEXP Y_,VectorXi j){
  if(Rf_isNull(X_)) return(as<Map<MatrixXd> >(Y_));
  if(Rf_isMatrix(X_)) {
    Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      MatrixXd Z = MatrixXd::Zero(X.rows(),Y.cols());
      for(int i=0;i<X.cols();i++){
        if(j[i]) {
          Z += X.col(i) * Y.row(i);
        }
      }
      return(Z);
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      SparseMatrix<double,RowMajor> Yr = Y;  // Convert to RowMajor so it can be parallelized
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      MatrixXd Z = MatrixXd::Zero(X.rows(),Y.cols());
      for(int i=0;i<X.cols();i++){
        if(j[i]) {
          Z += X.col(i) * Yr.row(i);
        }
      }
      return(Z);
    }
  }
  else {
    MSpMat X = as<MSpMat>(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      setNbThreads(0);
      if(Eigen::nbThreads( ) > 1) {
        MatrixXd Z = MatrixXd::Zero(X.rows(),Y.cols());
        for(int i=0;i<X.cols();i++){
          if(j[i]) {
            Z += X.col(i) * Y.row(i);
          }
        }
        return(Z);
      } else{
        MatrixXd Z = MatrixXd::Zero(X.rows(),Y.cols());
        for(int i=0;i<X.cols();i++){
          if(j[i]) {
            Z += X.col(i) * Y.row(i);
          }
        }
        return(Z);
      }
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      SparseMatrix<double,RowMajor> Yr = Y;  // Convert to RowMajor so it can be parallelized
      MatrixXd Z = MatrixXd::Zero(X.rows(),Y.cols());
      for(int i=0;i<X.cols();i++){
        if(j[i]) {
          Z += X.col(i) * Yr.row(i);
        }
      }
      return(Z);
    }
  }
}



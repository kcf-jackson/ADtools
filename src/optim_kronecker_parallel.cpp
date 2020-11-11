#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix lcircledast2(NumericMatrix x, NumericMatrix y) {
  int num_x_blk_out = x.nrow();
  int num_x_blk_in = x.ncol();
  int size_x_blk = y.nrow() / num_x_blk_in;
  int size_y_blk = y.ncol();
  NumericMatrix retval(num_x_blk_out * size_x_blk, size_y_blk);

  for (int i = 0; i < num_x_blk_out; i++) {
    NumericMatrix::Sub sub = retval( Range(i*size_x_blk, (i+1)*size_x_blk-1) , _ );
    for (int j = 0; j < num_x_blk_in; j++) {
      for (int h = 0; h < sub.nrow(); h++) {
        for (int k = 0; k < sub.ncol(); k++) {
          sub(h, k) += x(i, j) * y(j * size_x_blk + h, k);
        }
      }
    }
  }
  return retval;
}


// [[Rcpp::export]]
NumericMatrix lboxdot2(NumericMatrix x, NumericMatrix y) {
  int size_x_blk_in = x.ncol();
  int size_x_blk_out = x.nrow();
  int size_y_blk = y.ncol();
  int num_x_blk = y.nrow() / size_x_blk_in;
  NumericMatrix retval(size_x_blk_out * num_x_blk, size_y_blk);

  for (int i = 0; i < num_x_blk; i++) {
    for (int h = 0; h < size_x_blk_out; h++) {
      for (int k = 0; k < retval.ncol(); k++) {
        for (int j = 0; j < size_x_blk_in; j++) {
          retval(size_x_blk_out*i + h, k) += x(h, j) * y(size_x_blk_in * i + j, k);
        }
      }
    }
  }
  return retval;
}

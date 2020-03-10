#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat lcircledast(arma::mat x, arma::mat y) {
  int num_x_blk_out = x.n_rows;
  int num_x_blk_in = x.n_cols;
  int size_x_blk = y.n_rows / num_x_blk_in;
  int size_y_blk = y.n_cols;
  arma::mat retval = arma::zeros(num_x_blk_out * size_x_blk, size_y_blk);

  for (int i = 0; i < num_x_blk_out; i++) {
    arma::mat tempval = arma::zeros(size_x_blk, size_y_blk);
    for (int j = 0; j < num_x_blk_in; j++) {
      tempval += x(i, j) * y(arma::span(j*size_x_blk, (j+1)*size_x_blk-1), arma::span::all);
    }
    retval(arma::span(i*size_x_blk, (i+1)*size_x_blk-1), arma::span::all) = tempval;
  }
  return retval;
}


// [[Rcpp::export]]
arma::mat lboxdot(arma::mat x, arma::mat y) {
  int size_x_blk_in = x.n_cols;
  int size_x_blk_out = x.n_rows;
  int size_y_blk = y.n_cols;
  int num_x_blk = y.n_rows / size_x_blk_in;

  arma::mat retval = arma::zeros(size_x_blk_out * num_x_blk, size_y_blk);
  for (int i = 0; i < num_x_blk; i++) {
    retval(arma::span(size_x_blk_out*i, size_x_blk_out*(i+1)-1), arma::span::all) =
      x * y(arma::span(size_x_blk_in*i, size_x_blk_in*(i+1)-1), arma::span::all);
  }
  return retval;
}


// [[Rcpp::export]]
arma::mat rcircledast(arma::mat x, arma::mat y) {
  int num_y_blk_in = y.n_rows;
  int num_y_blk_out = y.n_cols;
  int size_x_blk = x.n_rows;
  int size_y_blk = x.n_cols / num_y_blk_in;
  arma::mat retval = arma::zeros(size_x_blk, size_y_blk * num_y_blk_out);

  for (int j = 0; j < num_y_blk_out; j++) {
    arma::mat tempval = arma::zeros(size_x_blk, size_y_blk);
    for (int i = 0; i < num_y_blk_in; i++) {
      tempval += x(arma::span::all, arma::span(i*size_y_blk, (i+1)*size_y_blk-1)) * y(i, j);
    }
    retval(arma::span::all, arma::span(j*size_y_blk, (j+1)*size_y_blk-1)) = tempval;
  }
  return retval;
}


// [[Rcpp::export]]
arma::mat rboxdot(arma::mat x, arma::mat y) {
  int size_x_blk = x.n_rows;
  int size_y_blk_in = y.n_rows;
  int size_y_blk_out = y.n_cols;
  int num_y_blk = x.n_cols / y.n_rows;
  arma::mat retval = arma::zeros(size_x_blk, num_y_blk * size_y_blk_out);

  for (int i = 0; i < num_y_blk; i++) {
    retval(arma::span::all, arma::span(size_y_blk_out*i, size_y_blk_out*(i+1)-1)) =
      x(arma::span::all, arma::span(size_y_blk_in*i, size_y_blk_in*(i+1)-1)) * y;
  }
  return retval;
}

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct lboxdotWorker : public Worker {
  const RMatrix<double> x;
  const RMatrix<double> y;
  RMatrix<double> rmat;

  int size_x_blk_in = x.ncol();
  int size_x_blk_out = x.nrow();
  int size_y_blk = y.ncol();
  int num_x_blk = y.nrow() / size_x_blk_in;

  lboxdotWorker(const NumericMatrix x, NumericMatrix y, NumericMatrix rmat)
  : x(x), y(y), rmat(rmat) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (int h = 0; h < size_x_blk_out; h++) {
        for (int k = 0; k < rmat.ncol(); k++) {
          for (int j = 0; j < size_x_blk_in; j++) {
            rmat(size_x_blk_out*i + h, k) += x(h, j) * y(size_x_blk_in * i + j, k);
          }
        }
      }
    }
  }
};


// [[Rcpp::export]]
NumericMatrix lboxdot(NumericMatrix x, NumericMatrix y) {
  int size_x_blk_in = x.ncol();
  int size_x_blk_out = x.nrow();
  int size_y_blk = y.ncol();
  int num_x_blk = y.nrow() / size_x_blk_in;
  NumericMatrix retval(size_x_blk_out * num_x_blk, size_y_blk);

  lboxdotWorker Loop(x, y, retval);
  parallelFor(0, num_x_blk, Loop);
  return retval;
}

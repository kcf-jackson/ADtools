// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct lcircledastWorker : public Worker {
  const RMatrix<double> mat;
  const RMatrix<double> mat2;
  RMatrix<double> rmat;

  int num_x_blk_out = mat.nrow();
  int num_x_blk_in = mat.ncol();
  int size_x_blk = mat2.nrow() / num_x_blk_in;
  int size_y_blk = mat2.ncol();

  lcircledastWorker(const NumericMatrix mat, NumericMatrix mat2, NumericMatrix rmat)
    : mat(mat), mat2(mat2), rmat(rmat) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < num_x_blk_in; j++) {
        for (int h = 0; h < size_x_blk; h++) {
          for (int k = 0; k < rmat.ncol(); k++) {
            rmat(i * size_x_blk + h, k) += mat(i, j) * mat2(j * size_x_blk + h, k);
          }
        }
      }
    }
  }
};


// [[Rcpp::export]]
NumericMatrix lcircledast(NumericMatrix x, NumericMatrix y) {
  int num_x_blk_out = x.nrow();
  int num_x_blk_in = x.ncol();
  int size_x_blk = y.nrow() / num_x_blk_in;
  int size_y_blk = y.ncol();
  NumericMatrix retval(num_x_blk_out * size_x_blk, size_y_blk);

  lcircledastWorker Loop(x, y, retval);
  parallelFor(0, num_x_blk_out, Loop);
  return retval;
}

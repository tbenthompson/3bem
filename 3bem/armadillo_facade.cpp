#include "armadillo_facade.h"
#include "dense_operator.h"
#include <armadillo>

namespace tbem {

arma::mat arma_mat(const Operator& op)
{
    return arma::mat(&op[0], op.n_rows(), op.n_cols());
}
Operator arma_invert(const Operator& op) 
{
    arma::mat inv_A = arma::inv(arma_mat(op));
    std::vector<double> out(inv_A.begin(), inv_A.end());
    return Operator(op.n_rows(), op.n_cols(), out);
}

double arma_cond(const Operator& op) {
    auto A = arma_mat(op);
    return arma::cond(A);
}

} // end namespace tbem

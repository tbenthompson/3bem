#include "armadillo_interface.h"
#include <armadillo>

std::vector<double> arma_invert(const std::vector<double>& vec_mat) 
{
    auto entries = vec_mat.size();
    auto rows = std::sqrt(entries);
    auto cols = rows;
    arma::mat A(&vec_mat[0], rows, cols);
    arma::mat inv_A = arma::inv(A);
    std::vector<double> out(inv_A.begin(), inv_A.end());
    return out;
}

std::vector<double> arma_invert_block(size_t n_comp_rows, size_t n_comp_cols,
    const std::vector<std::vector<double>>& block_mat) 
{

}

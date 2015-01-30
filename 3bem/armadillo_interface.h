#ifndef __kljewtlkha_ARMADILLO_H
#define __kljewtlkha_ARMADILLO_H

#include <vector>
#include <cstdlib>

std::vector<double> arma_invert(const std::vector<double>& vec_mat);

std::vector<double> arma_invert_block(size_t n_comp_rows, size_t n_comp_cols,
    const std::vector<std::vector<double>>& block_mat);

#endif

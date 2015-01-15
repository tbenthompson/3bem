#ifndef __123123123789798798_OPERATOR_H
#define __123123123789798798_OPERATOR_H

#include <cassert>
#include <vector>
#include "numbers.h"
//TODO: Try to remove this dependency on vec_ops
#include "vec_ops.h"

namespace tbem {

struct MatrixOperator 
{
    //TODO: rename to n_rows
    const size_t rows;
    //TODO: rename to n_cols
    const size_t cols;
    const size_t n_comp_rows;
    const size_t n_comp_cols;
    const std::vector<std::vector<double>> data;
};


MatrixOperator reshape_to_operator(const size_t rows, const size_t cols,
    const std::vector<double>& A) 
{
    return {rows, cols, 1, 1, {A}};
}

template <size_t dim>
MatrixOperator reshape_to_operator(const size_t rows, const size_t cols, 
    std::vector<Vec<Vec<double,dim>,dim>> A) 
{
    std::vector<std::vector<double>> out(dim * dim, std::vector<double>(A.size()));
    for (size_t i = 0; i < A.size(); i++) {
        for (int d1 = 0; d1 < dim; d1++) {
            for (int d2 = 0; d2 < dim; d2++) {
                out[d1 * dim + d2][i] = A[i][d1][d2];
            }
        }
    }
    return {rows, cols, dim, dim, out};
}

std::vector<std::vector<double>> apply_operator(const MatrixOperator& A,
    const std::vector<std::vector<double>>& x) 
{
    std::vector<std::vector<double>> res(A.n_comp_rows, std::vector<double>(A.rows, 0.0));
    for (size_t d1 = 0; d1 < A.n_comp_rows; d1++) {
        for (size_t d2 = 0; d2 < A.n_comp_cols; d2++) {
            size_t comp_idx = d1 * A.n_comp_cols + d2;
            assert(A.rows * x[d2].size() == A.data[comp_idx].size());
#pragma omp parallel for
            for (size_t i = 0; i < A.rows; i++) {
                for (size_t j = 0; j < A.cols; j++) {
                    size_t matrix_idx = i * A.cols + j;
                    res[d1][i] += A.data[comp_idx][matrix_idx] * x[d2][j];
                }
            }
        }
    }
    return res;
}

std::vector<double> apply_operator(const MatrixOperator& A, const std::vector<double>& x) {
    assert(A.n_comp_cols == 1 && A.n_comp_rows == 1);
    return apply_operator(A, std::vector<std::vector<double>>{x})[0];
}

} // end namespace tbem

#endif

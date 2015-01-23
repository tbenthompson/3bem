#include "operator.h"

namespace tbem {

MatrixOperator reshape_to_operator(const size_t n_rows, const size_t n_cols,
    const std::vector<double>& A) 
{
    return {n_rows, n_cols, 1, 1, {A}};
}

template <size_t dim>
MatrixOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
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
    return {n_rows, n_cols, dim, dim, out};
}

template 
MatrixOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    std::vector<Vec<Vec<double,2>,2>> A);
template 
MatrixOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    std::vector<Vec<Vec<double,3>,3>> A); 

std::vector<std::vector<double>> apply_operator(const MatrixOperator& A,
    const std::vector<std::vector<double>>& x) 
{
    assert(A.n_comp_rows * x.size() == A.data.size());
    assert(A.n_comp_cols == x.size());

    std::vector<std::vector<double>> res(A.n_comp_rows,
        std::vector<double>(A.n_rows, 0.0));
    for (size_t d1 = 0; d1 < A.n_comp_rows; d1++) {
        for (size_t d2 = 0; d2 < A.n_comp_cols; d2++) {
            size_t comp_idx = d1 * A.n_comp_cols + d2;
            assert(A.n_rows * x[d2].size() == A.data[comp_idx].size());
#pragma omp parallel for
            for (size_t i = 0; i < A.n_rows; i++) {
                for (size_t j = 0; j < A.n_cols; j++) {
                    size_t matrix_idx = i * A.n_cols + j;
                    res[d1][i] += A.data[comp_idx][matrix_idx] * x[d2][j];
                }
            }
        }
    }
    return res;
}

std::vector<double> apply_operator(const MatrixOperator& A, const std::vector<double>& x) 
{
    assert(A.n_comp_cols == 1 && A.n_comp_rows == 1);
    return apply_operator(A, std::vector<std::vector<double>>{x})[0];
}

} //end namespace tbem

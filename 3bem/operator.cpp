#include "operator.h"

namespace tbem {

BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols,
    const std::vector<double>& A) 
{
    return {1, 1, {{n_rows, n_cols, A}}};
}

template <size_t dim>
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    std::vector<Vec<Vec<double,dim>,dim>> A) 
{
    std::vector<std::vector<double>> out(dim * dim, std::vector<double>(A.size()));
    std::vector<Operator> ops;
    for (int d1 = 0; d1 < dim; d1++) {
        for (int d2 = 0; d2 < dim; d2++) {
            std::vector<double> data(n_rows * n_cols);
            for (size_t i = 0; i < A.size(); i++) {
                data[i] = A[i][d1][d2];
            }
            ops.push_back({n_rows, n_cols, data});
        }
    }
    return {dim, dim, ops};
}

template 
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    std::vector<Vec<Vec<double,2>,2>> A);
template 
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    std::vector<Vec<Vec<double,3>,3>> A); 

BlockFunction apply_operator(const BlockOperator& A, const BlockFunction& x)
{
    assert(A.n_comp_rows * x.size() == A.data.size());
    assert(A.n_comp_cols == x.size());

    BlockFunction res(A.n_comp_rows, Function(A.ops[0].n_rows, 0.0));
    for (size_t d1 = 0; d1 < A.n_comp_rows; d1++) {
        for (size_t d2 = 0; d2 < A.n_comp_cols; d2++) {
            size_t comp_idx = d1 * A.n_comp_cols + d2;

            const auto& op = A.ops[comp_idx];
            assert(op.n_rows * x[d2].size() == op.size());
            assert(op.n_cols == x[d2].size());

#pragma omp parallel for
            for (size_t i = 0; i < op.n_rows; i++) {
                for (size_t j = 0; j < op.n_cols; j++) {
                    size_t matrix_idx = i * op.n_cols + j;
                    res[d1][i] += op.data[matrix_idx] * x[d2][j];
                }
            }
        }
    }
    return res;
}

Function apply_operator(const BlockOperator& A, const Function& x) 
{
    assert(A.n_comp_cols == 1 && A.n_comp_rows == 1);
    return apply_operator(A, BlockFunction{x})[0];
}

BlockOperator combine_components(const BlockOperator& block_op) {
    size_t total_cols = block_op.n_comp_cols * block_op.ops[0].n_cols;
    size_t total_rows = block_op.n_comp_rows * block_op.ops[0].n_rows;
    int n_elements = total_cols * total_rows;

    std::vector<double> out(n_elements);

    for (size_t d1 = 0; d1 < block_op.n_comp_rows; d1++) {
        for (size_t row_idx = 0; row_idx < block_op.ops[0].n_rows; row_idx++) {
            size_t out_row_idx = d1 * block_op.ops[0].n_rows + row_idx;

            for (size_t d2 = 0; d2 < block_op.n_comp_cols; d2++) {
                size_t in_comp = d1 * block_op.n_comp_cols + d2;

                const auto& op = block_op.ops[in_comp];
                for (size_t col_idx = 0; col_idx < op.n_cols; col_idx++) {

                    size_t out_col_idx = d2 * op.n_cols + col_idx;
                    size_t in_element = row_idx * op.n_cols + col_idx;
                    size_t out_element = out_row_idx * total_cols + out_col_idx;
                    out[out_element] = op.data[in_element];

                }
            }
        }
    }

    return {1, 1, {{total_rows, total_cols, out}}};
}


} //end namespace tbem

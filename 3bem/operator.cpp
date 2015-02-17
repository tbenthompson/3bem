#include "operator.h"

namespace tbem {

Operator make_operator(size_t n_rows, size_t n_cols, const std::vector<double>& data) 
{
    return {n_rows, n_cols, std::make_shared<std::vector<double>>(data)};
}

BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols,
    const std::vector<double>& A) 
{
    return {1, 1, {make_operator(n_rows, n_cols, A)}};
}

template <size_t dim>
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    const std::vector<Vec<Vec<double,dim>,dim>>& A) 
{
    std::vector<Operator> ops;
    for (int d1 = 0; d1 < dim; d1++) {
        for (int d2 = 0; d2 < dim; d2++) {
            auto data = std::make_shared<std::vector<double>>(n_rows * n_cols);
            for (size_t i = 0; i < A.size(); i++) {
                (*data)[i] = A[i][d1][d2];
            }
            ops.push_back({n_rows, n_cols, data});
        }
    }
    return {dim, dim, ops};
}

template 
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    const std::vector<Vec<Vec<double,2>,2>>& A);
template 
BlockOperator reshape_to_operator(const size_t n_rows, const size_t n_cols, 
    const std::vector<Vec<Vec<double,3>,3>>& A); 

BlockFunction apply_operator(const BlockOperator& A, const BlockFunction& x)
{
    assert(A.n_comp_rows * x.size() == A.ops.size());
    assert(A.n_comp_cols == x.size());

    BlockFunction res(A.n_comp_rows, Function(A.ops[0].n_rows, 0.0));
    for (size_t d1 = 0; d1 < A.n_comp_rows; d1++) {
        for (size_t d2 = 0; d2 < A.n_comp_cols; d2++) {
            size_t comp_idx = d1 * A.n_comp_cols + d2;

            const auto& op = A.ops[comp_idx];
            assert(op.n_rows * x[d2].size() == op.data.size());
            assert(op.n_cols == x[d2].size());

#pragma omp parallel for
            for (size_t i = 0; i < op.n_rows; i++) {
                for (size_t j = 0; j < op.n_cols; j++) {
                    size_t matrix_idx = i * op.n_cols + j;
                    res[d1][i] += (*op.data)[matrix_idx] * x[d2][j];
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

size_t total_cols(const BlockOperator& block_op) {
    size_t n_cols = 0;
    for (size_t d2 = 0; d2 < block_op.n_comp_cols; d2++) {
        n_cols += block_op.ops[d2].n_cols;
    }
    return n_cols;
}

size_t total_rows(const BlockOperator& block_op) {
    size_t n_rows = 0;
    for (size_t d1 = 0; d1 < block_op.n_comp_rows; d1++) {
        n_rows += block_op.ops[d1 * block_op.n_comp_cols].n_rows;
    }
    return n_rows;
}

BlockOperator combine_components(const BlockOperator& block_op) {
    auto n_cols = total_cols(block_op);
    auto n_rows = total_rows(block_op);
    auto n_elements = n_cols * n_rows;

    auto out = std::make_shared<std::vector<double>>(n_elements);

    size_t n_rows_so_far = 0;
    for (size_t d1 = 0; d1 < block_op.n_comp_rows; d1++) {
        auto n_this_comp_rows = block_op.ops[d1 * block_op.n_comp_cols].n_rows;
        for (size_t row_idx = 0; row_idx < n_this_comp_rows; row_idx++) {
            auto out_row_idx = n_rows_so_far + row_idx;

            size_t n_cols_so_far = 0;
            for (size_t d2 = 0; d2 < block_op.n_comp_cols; d2++) {
                auto in_comp = d1 * block_op.n_comp_cols + d2;

                const auto& op = block_op.ops[in_comp];
                for (size_t col_idx = 0; col_idx < op.n_cols; col_idx++) {
                    auto out_col_idx = n_cols_so_far + col_idx;

                    auto in_element = row_idx * op.n_cols + col_idx;
                    auto out_element = out_row_idx * n_cols + out_col_idx;
                    (*out)[out_element] = (*op.data)[in_element];
                }
                n_cols_so_far += op.n_cols;
            }
        }
        n_rows_so_far += n_this_comp_rows;
    }

    return {1, 1, {{n_rows, n_cols, out}}};
}


} //end namespace tbem
